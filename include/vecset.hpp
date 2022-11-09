//
// Created by max on 09.11.22.
//
#ifndef PDS_VECSET_HPP
#define PDS_VECSET_HPP

#include <vector>
#include <memory>
#include <type_traits>
#include <cassert>

namespace pds {

template<std::integral T, std::integral Timestamp = bool>
requires std::is_trivially_destructible_v<Timestamp> && std::is_constructible_v<T, size_t>
class VecSet {
    static constexpr const Timestamp MIN_TIMESTAMP = std::numeric_limits<Timestamp>::min();
    static constexpr const Timestamp INITIAL_TIMESTAMP = std::numeric_limits<Timestamp>::min() + 1;
    static constexpr const Timestamp MAX_TIMESTAMP = std::numeric_limits<Timestamp>::max();

    struct Iterator {
        const VecSet* base;
        size_t i;

        using value_type = T;
        using reference = void;
        using pointer = std::unique_ptr<T>;
        using difference_type = ssize_t;
        using iterator_category = std::forward_iterator_tag;

        Iterator() = default;
        Iterator(const Iterator&) = default;
        Iterator(Iterator&&) = default;
        Iterator& operator=(const Iterator&) = default;
        Iterator& operator=(Iterator&&) = default;

        Iterator(const VecSet* base, size_t i) : base(base), i(i) {
            while (i < base->capacity() && !base->contains(i)) { ++i; }
        }

        auto operator<=>(const Iterator&) const noexcept = default;

        value_type operator*() const {
            return value_type{i};
        }

        pointer operator->() const {
            return std::make_unique<value_type>(i);
        }

        Iterator& operator++() {
            do {
                ++i;
            } while (i < base->capacity() && !base->contains(i));
            return *this;
        }

        Iterator& operator++(int) {
            auto old = *this;
            ++this;
            return old;
        }
    };

    inline bool validKey(const T &x) const { return x >= 0 && static_cast<size_t>(x) < m_present.size(); }

    void makeValidKey(const T& k) {
        assert(k >= 0);
        reserve(static_cast<size_t>(k) + 1);
    }
public:
    using iterator = Iterator;
    using const_iterator = iterator;

    VecSet(size_t capacity) : m_size {0}, m_present(capacity, MIN_TIMESTAMP), m_timestamp(INITIAL_TIMESTAMP) { }
    VecSet() : VecSet(0) { }

    VecSet(const VecSet &other) = default;

    VecSet(VecSet &&other) = default;

    VecSet &operator=(const VecSet &) = default;

    VecSet &operator=(VecSet &&) = default;

    inline bool empty() const { return m_size == 0; }
    inline size_t size() const { return m_size; }

    void reserve(size_t maxCapacity) {
        if (capacity() < maxCapacity) {
            m_present.resize(maxCapacity, MIN_TIMESTAMP);
        }
    }

    void shrink_to(size_t maxCapacity) {
        if (capacity() > maxCapacity) {
            size_t size = 0;
            for (size_t i = 0; i < maxCapacity; ++i) {
                if (contains(i)) {
                    ++size;
                    m_present[i] = INITIAL_TIMESTAMP;
                } else {
                    m_present[i] = MIN_TIMESTAMP;
                }
            }
            m_size = size;
            m_timestamp = INITIAL_TIMESTAMP;
            m_present.resize(maxCapacity, MIN_TIMESTAMP);
            m_present.shrink_to_fit();
        }
    }

    const_iterator begin() const {
        if (empty()) return end();
        else return {this, 0};
    }

    const_iterator end() const {
        return {this, capacity()};
    }

    inline size_t capacity() const noexcept { return m_present.size(); }

    void clear() {
        if (m_timestamp == MAX_TIMESTAMP) {
            m_present.clear();
            m_size = 0;
        } else {
            m_timestamp = m_timestamp + 1;
            m_size = 0;
        }
    }

    void insert(const T &x) {
        makeValidKey(x);
        if (!contains(x)) {
            size_t key = static_cast<size_t>(x);
            m_present[key] = m_timestamp;
            ++m_size;
        }
    }

    template<class... Args>
    void emplace(Args&&... args) {
        T key(std::forward<Args>(args)...);
        insert(key);
    }

    inline void erase(const T &x) {
        if (validKey(x) && m_present[x]) {
            m_present[x] = false;
            --m_size;
        }
    }

    inline void swap(VecSet& other) {
        std::swap(m_present, other.m_present);
        std::swap(m_timestamp, other.m_timestamp);
        std::swap(m_size, other.m_size);
    }

    inline bool contains(const T &x) const { return validKey(x) && m_present[x] == m_timestamp; }

    inline size_t count(const T &x) const { return contains(x); }

    // TODO find

private:
    size_t m_size;
    std::vector <Timestamp> m_present;
    Timestamp m_timestamp;
};
}

#endif //PDS_VECSET_HPP
