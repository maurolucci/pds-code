//
// Created by max on 08.11.22.
//
#ifndef PDS_VECMAP_HPP
#define PDS_VECMAP_HPP

#include <vector>
#include <memory>
#include <type_traits>
#include <cassert>

namespace pds {
template<std::integral T, class Bool = bool>
class VecSet {
public:
    VecSet() : m_size{0}, m_content{} {};

    VecSet(const VecSet &other) = default;

    VecSet(VecSet &&other) = default;

    VecSet &operator=(const VecSet &) = default;

    VecSet &operator=(VecSet &&) = default;

    void clear() {
        m_content.clear();
    }

    void insert(const T &x) {
        assert(x > 0);
        if (m_content.size() <= x) {
            m_content.resize(x + 1, false);
        }
        if (!m_content[x]) {
            m_content[x] = true;
            ++m_size;
        }
    }

    inline void erase(const T &x) {
        if (validKey(x) && m_content[x]) {
            m_content[x] = false;
            --m_size;
        }
    }

    inline bool validKey(const T &x) const { return x >= 0 && static_cast<size_t>(x) < m_content.size(); }

    inline bool contains(const T &x) const { validKey(x) && m_content[x]; }

    inline size_t count(const T &x) const { return contains(x); }

private:
    std::vector <Bool> m_content;
    size_t m_size;
};

template<std::integral K, class V, std::integral Timestamp = bool, class Allocator = std::allocator <V>>
class VecMap {
    union Entry {
        struct {} empty;
        V value;

        Entry() {};

        /// only safe to be called with empty Entries!
        Entry(const Entry &) {};

        /// only safe to be called with empty Entries!
        Entry &operator=(const Entry &) { return *this; }

        ~Entry() {};

        void clear() { std::destroy_at(&value); }

        template<class... T>
        void updateValue(T... args) {
            clear();
            emplaceValue(args...);
        }

        template<class... T>
        void emplaceValue(T &&... args) {
            std::construct_at(&value, std::forward<T>(args)...);
        }

    };

    template<class Map>
    struct IteratorPair {
        size_t first;
        std::conditional_t<std::is_const_v < Map>, const V&, V&> second;
    };

    template<class Map>
    struct Iterator {
        Map *base;
        size_t i;

        using value_type = IteratorPair<Map>;
        using reference = value_type;
        using pointer = std::unique_ptr<IteratorPair<Map>>;
        using difference_type = ssize_t;
        using iterator_category = std::forward_iterator_tag;

        auto operator<=>(const Iterator &other) const = default;

        Iterator() : base(nullptr), i(0) { }
        Iterator(const Iterator<VecMap> &other) requires std::is_const_v<Map>: base(other.base), i(other.i) {}
        Iterator(Iterator&&) = default;
        Iterator(const Iterator&) = default;
        Iterator& operator=(const Iterator&) = default;
        Iterator& operator=(Iterator&&) = default;

        Iterator(Map *base, size_t i) : base(base), i(i) {
            while (i < base->m_present.size() && !base->m_present[i]) {
                ++i;
            }
        }

        IteratorPair<Map> operator*() requires (!std::is_const_v < Map > ) {
            return {i, const_cast<VecMap *>(base)->m_content[i]};
        }

        IteratorPair<const Map> operator*() const {
            return {i, base->m_present[i]};
        }

        pointer operator->() requires (!std::is_const_v < Map > ) {
            return std::make_unique<IteratorPair<Map>>(IteratorPair<Map>{i, base->m_present[i]});
        }

        pointer operator->() const {
            return std::make_unique<IteratorPair<Map>>(IteratorPair<Map>{i, base->m_present[i]});
        }

        Iterator &operator++() {
            do {
                ++i;
            } while (i < base->capacity() && !base->contains(i));
            return *this;
        }

        Iterator operator++(int) {
            auto old = *this;
            ++*this;
            return old;
        }
    };

public:
    using iterator = Iterator<VecMap>;
    using const_iterator = Iterator<const VecMap>;

private:
    inline bool validKey(const K &x) const { return x >= 0 && static_cast<size_t>(x) < capacity(); }

    inline void makeValidKey(const K &s) {
        if (s < 0) {
            throw std::out_of_range("negative key");
        }
        reserve(s + 1);
    }

    void clearEntryUnchecked(size_t k) {
        std::destroy_at(&m_content[k]);
        if (m_present[k] == m_timestamp) {
            --m_size;
        }
        m_present[k] = MIN_TIMESTAMP;
    }

    template<class... T>
    iterator emplaceNewEntryUnchecked(size_t k, T &&... args) {
        std::construct_at<V>(&m_content[k], std::forward<T>(args)...);
        m_present[k] = m_timestamp;
        ++m_size;
        return {this, k};
    }

    bool clearEntry(size_t k) {
        if (m_present[k] != MIN_TIMESTAMP) {
            clearEntryUnchecked(k);
            return true;
        } else {
            return false;
        }
    }

    template<class... T>
    std::pair<iterator, bool> updateEntry(const K &k, T &&... args) {
        assert(k >= 0);
        makeValidKey(k);
        if (m_present[k] != MIN_TIMESTAMP) {
            m_content[k] = V(std::forward<T>(args)...);
            if (m_present[k] != m_timestamp) {
                m_present[k] = m_timestamp;
                ++m_size;
            }
            return {{this, k}, false};
        } else {
            emplaceNewEntryUnchecked(static_cast<size_t>(k), std::forward<T>(args)...);
            return {{this, k}, true};
        }
    }

    static constexpr const Timestamp MIN_TIMESTAMP = std::numeric_limits<Timestamp>::min();
    static constexpr const Timestamp INITIAL_TIMESTAMP = std::numeric_limits<Timestamp>::min() + 1;
    static constexpr const Timestamp MAX_TIMESTAMP = std::numeric_limits<Timestamp>::max();

    void moveFront(V* newContent, size_t count) {
        for (size_t i = 0; i < count; ++i) {
            if (m_present[i] == m_timestamp) {
                std::construct_at(&newContent[i], std::move(m_content[i]));
                std::destroy_at(&m_content[i]);
                m_present[i] = INITIAL_TIMESTAMP;
            } else if (m_present[i] != MIN_TIMESTAMP) {
                std::destroy_at(&m_content[i]);
            }
        }
        m_timestamp = INITIAL_TIMESTAMP;
    }

    // must be called before moveFront! (moveFront invalidates timestamps)
    void destroyTail(size_t start) {
        for (size_t i = start; i < capacity(); ++i) {
            clearEntry(i);
        }
    }


public:

    explicit VecMap(size_t capacity) :
            m_alloc{},
            m_timestamp(INITIAL_TIMESTAMP),
            m_capacity(0),
            m_present(0),
            m_content(nullptr),
            m_size(0)
    {
        reserve(capacity);
        assert(this->capacity() >= size());
    }

    VecMap() : VecMap(size_t{0}) {}

    virtual ~VecMap() {
        destroyTail(0);
        m_alloc.deallocate(m_content, capacity());
    }

    VecMap(const VecMap &other) : VecMap(other.capacity()) {
        for (size_t i = 0; i < other.capacity(); ++i) {
            if (other.m_present[i]) {
                emplaceNewEntryUnchecked(i, other.m_content[i]);
            }
        }
    }

    VecMap(VecMap &&other) noexcept: VecMap() {
        swap(other);
    }

    VecMap& operator=(const VecMap& other) {
        assert(capacity() >= size());
        assert(other.capacity() >= other.size());
        reserve(other.capacity());
        for (size_t i = 0; i < other.capacity(); ++i) {
            if (other.m_present[i]) {
                emplaceNewEntryUnchecked(i, other.m_content[i]);
            }
        }
        return *this;
    }

    VecMap& operator=(VecMap&& other) {
        swap(other);
        return *this;
    }

    iterator begin() noexcept {
        return {this, 0};
    }

    iterator end() noexcept {
        return {this, m_present.size()};
    }

    const_iterator begin() const noexcept {
        return {this, 0};
    }

    const_iterator end() const noexcept {
        return {this, m_present.size()};
    }

    friend auto begin(VecMap& map) { return map.begin(); }
    friend auto begin(const VecMap& map) { return map.begin(); }
    friend auto end(VecMap& map) { return map.end(); }
    friend auto end(const VecMap& map) { return map.end(); }

    void reserve(size_t minCapacity) {
        if (capacity() < static_cast<size_t>(minCapacity)) {
            size_t oldCapacity = capacity();
            m_present.resize(minCapacity, MIN_TIMESTAMP);
            assert(m_present.capacity() >= minCapacity);
            m_capacity = m_present.capacity();
            m_present.resize(m_capacity, MIN_TIMESTAMP);
            if (m_capacity != oldCapacity) {
                V *newContent = m_alloc.allocate(m_capacity);
                moveFront(newContent, oldCapacity);
                m_alloc.deallocate(m_content, oldCapacity);
                m_content = newContent;
            }
        }
        assert(capacity() >= size());
    }

    void shrink_to(size_t newCapacity) {
        assert(newCapacity > 0);
        if (capacity() > newCapacity) {
            destroyTail(newCapacity);
            size_t oldCapacity = capacity();
            m_present.resize(newCapacity, MIN_TIMESTAMP);
            m_present.shrink_to_fit();
            m_capacity = m_present.capacity();
            if (m_capacity != oldCapacity) {
                V* newContent = m_alloc.allocate(capacity());
                moveFront(newContent, m_capacity);
                m_alloc.deallocate(m_content, oldCapacity);
                m_content = newContent;
            }
        }
        assert(capacity() >= size());
    }

    inline size_t size() const noexcept { return m_size; }

    inline size_t capacity() const noexcept { return m_capacity; }

    inline bool empty() const noexcept { return m_size == 0; }

    inline bool contains(const K &k) const noexcept {
        return validKey(k) && m_present[k] == m_timestamp;
    }

    size_t count(const K &k) const noexcept {
        return contains(k);
    }

    iterator find(const K &key) noexcept {
        if (contains(key)) return {this, key};
        else return end();
    }

    const_iterator find(const K &key) const noexcept {
        if (contains(key)) return {this, key};
        else return end();
    }

    void clear() {
        if (empty()) return;
        if (std::is_same_v<bool, Timestamp>
                || m_timestamp == std::numeric_limits<Timestamp>::max())
        {
            assert(m_size == size_t(std::ranges::count(m_present, m_timestamp)));
            destroyTail(0);
            m_timestamp = INITIAL_TIMESTAMP;
            std::fill(m_present.begin(), m_present.end(), MIN_TIMESTAMP);
            assert(m_size == 0);
        } else {
            m_timestamp = m_timestamp + 1;
            m_size = 0;
        }
    }

    std::pair<iterator, bool> insert(const std::pair <K, V> &content) requires std::is_copy_constructible_v<V> {
        return emplace(content.first, content.second);
    }

    template<class M>
    std::pair<iterator, bool> insert_or_assign(const K &key, M &&value) requires std::is_constructible_v<V, M> {
        return updateEntry(key, value);
    }

    template<class ...T>
    std::pair<iterator, bool> emplace(const K &k, T &&... args) requires std::is_constructible_v<V, T...> {
        if (!contains(k)) {
            if (validKey(k)) {
                clearEntry(k);
            }else {
                makeValidKey(k);
            }
            return {emplaceNewEntryUnchecked(k, std::forward<T>(args)...), true};
        } else {
            return {{this, static_cast<size_t>(k)}, false};
        }
    }

    template<class... T>
    std::pair<iterator, bool> try_emplace(const K &k, T &&... args) requires std::is_constructible_v<V, T...> {
        return emplace(k, std::forward<T>(args)...);
    }

    size_t erase(const K &k) {
        return clearEntry(k);
    }

    void erase(const_iterator it) {
        assert(it.base == this);
        assert(contains(it.i));
        clearEntryUnchecked(it.i);
    }

    void swap(VecMap &other) noexcept(std::is_nothrow_swappable_v<Timestamp> && std::is_nothrow_swappable_v<decltype(m_present)>) {
        std::swap(m_alloc, other.m_alloc);
        std::swap(m_capacity, other.m_capacity);
        std::swap(m_timestamp, other.m_timestamp);
        std::swap(m_present, other.m_present);
        std::swap(m_content, other.m_content);
        std::swap(m_size, other.m_size);
    }

    const V &at(const K &k) const {
        if (!contains(k)) throw std::out_of_range("invalid key");
        return m_content[k];
    }

    V &at(const K &k) {
        if (!contains(k)) throw std::out_of_range("invalid key");
        return m_content[k];
    }

    const V &operator[](const K &k) const requires std::is_default_constructible_v<V> {
        if (!contains(k))
            emplace(k);
        return m_content[k];
    }

    V &operator[](const K &k) requires std::is_default_constructible_v<V> {
        if (!contains(k))
            emplace(k);
        return m_content[k];
    }

private:
    Allocator m_alloc;
    Timestamp m_timestamp;
    size_t m_capacity;
    std::vector <Timestamp> m_present;
    V *m_content;
    size_t m_size;
};

} // namespace pds

#endif //PDS_VECMAP_HPP
