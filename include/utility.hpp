

#ifndef PDS_UTILITY_HPP
#define PDS_UTILITY_HPP

#include <range/v3/view/subrange.hpp>

namespace pds {

template<class IteratorPair>
auto range_pair(IteratorPair pair) {
    return ranges::subrange(pair.first, pair.second);
}

}

#endif //PDS_UTILITY_HPP