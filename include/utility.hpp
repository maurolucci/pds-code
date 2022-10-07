//
// Created by max on 19.08.22.
//

#ifndef PDS_UTILITY_HPP
#define PDS_UTILITY_HPP

#include <algorithm>
namespace pds {

template<class... T>
void unused(T&&...) { }

template<typename CharT>
inline void ltrim(std::basic_string<CharT>& s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](const CharT c) { return !std::isspace(c); }));
}

template<typename CharT>
inline void rtrim(std::basic_string<CharT>& s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), [](const CharT c) { return !std::isspace(c); }).base(), s.end());
}

template<typename CharT>
inline void trim(std::basic_string<CharT>& s) {
    ltrim(s);
    rtrim(s);
}
}

#endif //PDS_UTILITY_HPP