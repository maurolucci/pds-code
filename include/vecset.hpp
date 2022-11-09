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

template<std::integral T, class Timestamp = bool>
class VecSet {
    static constexpr const Timestamp MIN_TIMESTAMP = std::numeric_limits<Timestamp>::min();
    static constexpr const Timestamp INITIAL_TIMESTAMP = std::numeric_limits<Timestamp>::min() + 1;
    static constexpr const Timestamp MAX_TIMESTAMP = std::numeric_limits<Timestamp>::max();
public:
    VecSet() : m_size{0}, m_present{} {};

    VecSet(const VecSet &other) = default;

    VecSet(VecSet &&other) = default;

    VecSet &operator=(const VecSet &) = default;

    VecSet &operator=(VecSet &&) = default;

    void clear() {
    }

    void insert(const T &x) {
        assert(x > 0);
        if (m_present.size() <= x) {
            m_present.resize(x + 1, false);
        }
        if (!m_present[x]) {
            m_present[x] = true;
            ++m_size;
        }
    }

    inline void erase(const T &x) {
        if (validKey(x) && m_present[x]) {
            m_present[x] = false;
            --m_size;
        }
    }

    inline bool validKey(const T &x) const { return x >= 0 && static_cast<size_t>(x) < m_present.size(); }

    inline bool contains(const T &x) const { validKey(x) && m_present[x]; }

    inline size_t count(const T &x) const { return contains(x); }

private:
    std::vector <Timestamp> m_present;
    Timestamp m_timestamp;
    size_t m_size;
};
}

#endif //PDS_VECSET_HPP
