//
// Created by max on 02.08.22.
//

#ifndef PDS_MAP_HPP
#define PDS_MAP_HPP

#include <unordered_map>
#include <unordered_set>

#include <ankerl/unordered_dense.h>

namespace pds {
template<class K, class V>
using map = ankerl::unordered_dense::map<K, V>;

template<class K>
using set = ankerl::unordered_dense::set<K>;

}

#endif //PDS_MAP_HPP
