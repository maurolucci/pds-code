//
// Created by max on 07.11.22.
//
#include <boost/test/unit_test.hpp>
#include <map>
#include "vecmap.hpp"
#include <fmt/format.h>
#include <range/v3/all.hpp>

struct DeStruct {
    static int& instances() {
        static int instances = 0;
        return instances;
    }
    static int& total() {
        static int total = 0;
        return total;
    }

    int value;
    int index;

    DeStruct() : DeStruct(-1) { }

    DeStruct(int value) : value(value), index(++total()) {
        ++instances();
    }

    DeStruct(const DeStruct& other) : DeStruct(other.value) { }

    virtual ~DeStruct() {
        --instances();
    }

    friend auto operator <=>(const DeStruct& a, const DeStruct& b);
};

BOOST_AUTO_TEST_CASE(test_DeStruct) {
    DeStruct::total() = 0;
    DeStruct::instances() = 0;
    DeStruct* reserved = (DeStruct*)malloc(sizeof(DeStruct));
    {
        DeStruct blub;
        BOOST_TEST(blub.index == 1);
        BOOST_TEST(DeStruct::instances() == 1);
        BOOST_TEST(DeStruct::total() == 1);
        DeStruct bla(blub);
        BOOST_TEST(bla.index == 2);
        BOOST_TEST(DeStruct::instances() == 2);
        BOOST_TEST(DeStruct::total() == 2);
        std::construct_at(reserved, 3);
        BOOST_TEST(reserved->index == 3);
        BOOST_TEST(DeStruct::instances() == 3);
        BOOST_TEST(DeStruct::total() == 3);
    }
    BOOST_TEST(DeStruct::instances() == 1);
    BOOST_TEST(DeStruct::total() == 3);
    std::destroy_at(reserved);
    BOOST_TEST(DeStruct::instances() == 0);
    BOOST_TEST(DeStruct::total() == 3);
    free(reserved);
}

template<class... T> struct PrintType;

BOOST_AUTO_TEST_CASE(traits) {
    using Map = pds::VecMap<size_t, DeStruct>;
    Map testMap;
    auto it = testMap.begin();
    using It = decltype(it);
    BOOST_TEST(std::input_or_output_iterator<It>);
    BOOST_TEST(std::weakly_incrementable<It>);
    //BOOST_TEST((requires (It m) { typename std::iter_difference_t<Map>; };));
    BOOST_TEST(std::movable<It>);
    BOOST_TEST(std::is_object_v<It>);
    BOOST_TEST(std::move_constructible<It>);
    BOOST_TEST((std::assignable_from<It&, It>));
    BOOST_TEST(std::swappable<It>);
    BOOST_TEST((std::sentinel_for<It, It>));
    auto begin = std::ranges::begin(testMap);
    auto end = std::ranges::end(testMap);
    BOOST_TEST((begin == end));
}

BOOST_AUTO_TEST_CASE(insert_delete) {
    DeStruct::total() = 0;
    DeStruct::instances() = 0;
    pds::VecMap<size_t, DeStruct> testMap;
    constexpr const size_t NUM_INSERT = 100;
    for (size_t i = 0; i < NUM_INSERT; ++i) {
        testMap.insert({i, int(i)});
        BOOST_TEST(testMap.contains(i), "expected " << i << " in map");
        BOOST_TEST(DeStruct::instances() == i + 1);
        BOOST_TEST(testMap[i].value == i);
        BOOST_TEST(testMap.size() == i + 1);
    }
    BOOST_TEST(DeStruct::instances() == NUM_INSERT);
    BOOST_TEST(testMap.size() == NUM_INSERT);
    for (size_t i = 0; i < NUM_INSERT; ++i) {
        BOOST_TEST(testMap.contains(i), "expected " << i << " in map");
        BOOST_TEST(testMap[i].value == i);
    }
    BOOST_TEST(DeStruct::instances() == NUM_INSERT);
    for (auto kv: testMap) {
        auto k = kv.first;
        auto v = kv.second;
        BOOST_TEST(testMap.contains(k));
        BOOST_TEST(testMap[k].value == v.value);
    }
    BOOST_TEST(DeStruct::instances() == NUM_INSERT);
    for (const auto& kv: static_cast<const decltype(testMap)&>(testMap)) {
        auto k = kv.first;
        auto& v = kv.second;
        BOOST_TEST(testMap.contains(k));
        BOOST_TEST(testMap[k].value == v.value);
    }
    BOOST_TEST(DeStruct::instances() == NUM_INSERT);
    for (auto it = begin(testMap); it != end(testMap); ++it) { }
    BOOST_TEST(DeStruct::instances() == NUM_INSERT);
    const auto& constMap = testMap;
    BOOST_TEST(DeStruct::instances() == NUM_INSERT);
    BOOST_TEST(std::distance(testMap.begin(), testMap.end()) == testMap.size());
    BOOST_TEST(DeStruct::instances() == NUM_INSERT);
    BOOST_TEST(std::distance(constMap.begin(), constMap.end()) == testMap.size());
    BOOST_TEST(DeStruct::instances() == NUM_INSERT);
    //PrintType<std::ranges::iterator_t<decltype(testMap)>, std::ranges::range_difference_t<decltype(testMap)>> { };

    //std::ranges::distance(testMap);
    //ranges::distance(testMap) == testMap.size();
    //for (auto [i, kv]: ranges::enumerate(testMap)) {
    //    auto k = kv.first;
    //    auto v = kv.second;
    //    BOOST_TEST(testMap[k] == v);
    //    BOOST_TEST(testMap.at(k) == v);
    //    BOOST_TEST(k == i);
    //}
    for (size_t i = 0; i < NUM_INSERT; ++i) {
        BOOST_TEST(testMap.contains(i), "expected " << i << " in map");
        BOOST_TEST(DeStruct::instances() == NUM_INSERT - i);
        testMap.erase(i);
        BOOST_TEST(!testMap.contains(i), "expected " << i << " to be deleted");
        BOOST_TEST(DeStruct::instances() == NUM_INSERT - i - 1);
        BOOST_TEST(testMap.size() == NUM_INSERT - i - 1);
    }
    BOOST_TEST(testMap.size() == 0);
    BOOST_TEST(DeStruct::instances() == 0);
}

void test_insert_delete(auto& testMap) {
    constexpr const size_t NUM_INSERT = 100;
    for (size_t i = 0; i < NUM_INSERT; ++i) {
        testMap.insert({i, int(i)});
        BOOST_TEST(testMap.contains(i), "expected " << i << " in map");
        BOOST_TEST(testMap[i].value == i);
        BOOST_TEST(testMap.size() == i + 1);
    }
    BOOST_TEST(testMap.size() == NUM_INSERT);
    for (size_t i = 0; i < NUM_INSERT; ++i) {
        BOOST_TEST(testMap.contains(i), "expected " << i << " in map");
        BOOST_TEST(testMap[i].value == i);
    }
    for (auto kv: testMap) {
        auto k = kv.first;
        auto v = kv.second;
        BOOST_TEST(testMap.contains(k));
        BOOST_TEST(testMap[k].value == v.value);
    }
    for (const auto& kv: static_cast<const decltype(testMap)&>(testMap)) {
        auto k = kv.first;
        auto& v = kv.second;
        BOOST_TEST(testMap.contains(k));
        BOOST_TEST(testMap[k].value == v.value);
    }
    for (auto it = begin(testMap); it != end(testMap); ++it) { }
    const auto& constMap = testMap;
    BOOST_TEST(std::distance(testMap.begin(), testMap.end()) == testMap.size());
    BOOST_TEST(std::distance(constMap.begin(), constMap.end()) == testMap.size());
    //PrintType<std::ranges::iterator_t<decltype(testMap)>, std::ranges::range_difference_t<decltype(testMap)>> { };

    //std::ranges::distance(testMap);
    //ranges::distance(testMap) == testMap.size();
    //for (auto [i, kv]: ranges::enumerate(testMap)) {
    //    auto k = kv.first;
    //    auto v = kv.second;
    //    BOOST_TEST(testMap[k] == v);
    //    BOOST_TEST(testMap.at(k) == v);
    //    BOOST_TEST(k == i);
    //}
    for (size_t i = 0; i < NUM_INSERT; ++i) {
        BOOST_TEST(testMap.contains(i), "expected " << i << " in map");
        testMap.erase(i);
        BOOST_TEST(!testMap.contains(i), "expected " << i << " to be deleted");
        BOOST_TEST(testMap.size() == NUM_INSERT - i - 1);
    }
    BOOST_TEST(testMap.size() == 0);
}

BOOST_AUTO_TEST_CASE(insert_delete_cleared) {
    BOOST_TEST_CONTEXT("Timestamp=bool"){
        DeStruct::total() = 0;
        DeStruct::instances() = 0;
        pds::VecMap<size_t, DeStruct, bool> testMap;
        test_insert_delete(testMap);
        testMap.clear();
        test_insert_delete(testMap);
    }
    BOOST_TEST(DeStruct::instances() == 0);
    BOOST_TEST_CONTEXT("Timestamp=uint8_t"){
        DeStruct::total() = 0;
        DeStruct::instances() = 0;
        pds::VecMap<size_t, DeStruct, std::uint8_t> testMap;
        test_insert_delete(testMap);
    }
    BOOST_TEST(DeStruct::instances() == 0);
    BOOST_TEST_CONTEXT("Timestamp=uint8_t, cleared map") {
        DeStruct::total() = 0;
        DeStruct::instances() = 0;
        pds::VecMap<size_t, DeStruct, std::uint8_t> testMap;
        test_insert_delete(testMap);
        for (size_t i = 0; i < 300; ++i) {
            testMap.emplace(0, 0);
            testMap.clear();
        }
        test_insert_delete(testMap);
    }
    BOOST_TEST(DeStruct::instances() == 0);
    BOOST_TEST_CONTEXT("Timestamp=char, cleared map") {
        DeStruct::total() = 0;
        DeStruct::instances() = 0;
        pds::VecMap<size_t, DeStruct, char> testMap;
        test_insert_delete(testMap);
        for (size_t i = 0; i < 300; ++i) {
            testMap.emplace(0, 0);
            testMap.clear();
        }
        test_insert_delete(testMap);
    }
    BOOST_TEST(DeStruct::instances() == 0);
}

void test_insert_delete2(auto& testMap) {
    constexpr const size_t NUM_INSERT = 100;
    for (size_t i = 0; i <= NUM_INSERT; ++i) {
        if (i < NUM_INSERT) {
            auto it = testMap.emplace(i, i);
            BOOST_TEST(testMap.size() == ((i > 0) ? 2 : 1));
            BOOST_TEST(it.second);
            BOOST_TEST(it.first->second.value == i);
        }
        if (i > 0) {
            BOOST_TEST(testMap.erase(i - 1));
            BOOST_TEST(!testMap.erase(i - 1));
        }
    }
}
BOOST_AUTO_TEST_CASE(insert_delete2) {
    BOOST_TEST_CONTEXT("Timestamp=bool"){
        DeStruct::total() = 0;
        DeStruct::instances() = 0;
        pds::VecMap<size_t, DeStruct, bool> testMap;
        test_insert_delete2(testMap);
    }
    BOOST_TEST(DeStruct::instances() == 0);
    BOOST_TEST_CONTEXT("Timestamp=uint8_t"){
        DeStruct::total() = 0;
        DeStruct::instances() = 0;
        pds::VecMap<size_t, DeStruct, std::uint8_t> testMap;
        test_insert_delete2(testMap);
    }
    BOOST_TEST(DeStruct::instances() == 0);
    BOOST_TEST_CONTEXT("Timestamp=uint8_t, cleared map") {
        DeStruct::total() = 0;
        DeStruct::instances() = 0;
        pds::VecMap<size_t, DeStruct, std::uint8_t> testMap;
        test_insert_delete2(testMap);
        for (size_t i = 0; i < 255; ++i) {
            testMap.emplace(0, 0);
            testMap.clear();
        }
        test_insert_delete2(testMap);
    }
    BOOST_TEST(DeStruct::instances() == 0);
    BOOST_TEST_CONTEXT("Timestamp=char, cleared map") {
        DeStruct::total() = 0;
        DeStruct::instances() = 0;
        pds::VecMap<size_t, DeStruct, char> testMap;
        test_insert_delete2(testMap);
        for (size_t i = 0; i < 255; ++i) {
            testMap.emplace(0, 0);
            testMap.clear();
        }
        test_insert_delete2(testMap);
    }
    BOOST_TEST(DeStruct::instances() == 0);
}

BOOST_AUTO_TEST_CASE(copy_move) {
    DeStruct::total() = 0;
    DeStruct::instances() = 0;
    constexpr size_t NUM_INSERT = 100;
    pds::VecMap<size_t, DeStruct> testMap(NUM_INSERT);
    BOOST_TEST(testMap.capacity() >= NUM_INSERT);
    for (size_t i = 0; i < NUM_INSERT; ++i) {
        auto it = testMap.try_emplace(2 * i, 2 * i);
        BOOST_TEST(testMap.size() == i + 1);
        BOOST_TEST(it.second);
        BOOST_TEST(it.first->second.value == 2 * i);
        BOOST_TEST(testMap.contains(2 * i));
    }
    BOOST_TEST(testMap.size() == NUM_INSERT);
    BOOST_TEST(DeStruct::instances() == testMap.size());
    auto testMap2 = testMap;
    BOOST_TEST(DeStruct::instances() == (2 * testMap.size()));
    BOOST_TEST(testMap.capacity() >= (NUM_INSERT * 2) - 1);
    BOOST_TEST(testMap2.size() == testMap.size());
    for (auto kv: testMap) {
        auto k = kv.first;
        auto& v = kv.second;
        BOOST_TEST(testMap2.contains(k));
        BOOST_TEST(testMap2.at(k).value == v.value);
    }
    for (auto kv: testMap2) {
        auto k = kv.first;
        auto& v = kv.second;
        BOOST_TEST(testMap.contains(k));
        BOOST_TEST(testMap.at(k).value == v.value);
    }
    decltype(testMap) testMap5 = testMap2;
    {
        decltype(testMap) testMap3;
        size_t totalCreated = DeStruct::total();
        testMap3 = std::move(testMap2);
        BOOST_TEST(testMap2.size() == 0);
        BOOST_TEST(testMap3.size() == testMap.size());
        BOOST_TEST(DeStruct::instances() == (3 * testMap.size()));
        BOOST_TEST(totalCreated == DeStruct::total());
        for (auto kv: testMap) {
            auto k = kv.first;
            auto &v = kv.second;
            BOOST_TEST(testMap3.contains(k));
            BOOST_TEST(testMap3.at(k).value == v.value);
        }
        std::swap(testMap2, testMap3);
        BOOST_TEST(testMap3.size() == 0);
        BOOST_TEST(testMap2.size() == testMap.size());
        BOOST_TEST(DeStruct::instances() == (3 * testMap.size()));
        BOOST_TEST(totalCreated == DeStruct::total());

        auto testMap4(std::move(testMap2));
        BOOST_TEST(testMap2.size() == 0);
        BOOST_TEST(testMap4.size() == testMap.size());
        BOOST_TEST(DeStruct::instances() == (3 * testMap.size()));
        BOOST_TEST(totalCreated == DeStruct::total());
        for (auto kv: testMap) {
            auto k = kv.first;
            auto &v = kv.second;
            BOOST_TEST(testMap4.contains(k));
            BOOST_TEST(testMap4.at(k).value == v.value);
            BOOST_TEST(testMap4.at(k).value == k);
        }

        testMap3 = testMap4;
        BOOST_TEST(DeStruct::instances() == (4 * testMap.size()));
        for (auto kv: testMap) {
            kv.second = 1337;
        }
        for (auto kv: testMap3) {
            auto k = kv.first;
            auto &v = kv.second;
            BOOST_TEST(testMap3.contains(k));
            BOOST_TEST(testMap3.at(k).value == v.value);
            BOOST_TEST(testMap3.at(k).value == k);
            BOOST_TEST(testMap.contains(k));
            BOOST_TEST(testMap.at(k).value == 1337);
        }
        BOOST_TEST(DeStruct::instances() == (4 * testMap.size()));
    } // destructs testMap3 and testMap4
    BOOST_TEST(DeStruct::instances() == (2 * testMap.size()));
    for (size_t i = 0; i < NUM_INSERT; ++i) {
        BOOST_TEST(testMap.contains(2 * i));
        BOOST_TEST(testMap[2*i].value == 1337);
        BOOST_TEST(testMap5.contains(2*i));
        BOOST_TEST(testMap5[2*i].value == (2*i));
    }
}

BOOST_AUTO_TEST_CASE(find_contains) {
    constexpr size_t NUM_INSERT = 100;
    pds::VecMap<size_t, DeStruct> testMap(NUM_INSERT);
    for (size_t i = 0; i < NUM_INSERT; ++i) {
        testMap.insert_or_assign(2*i, 2*i);
        auto it = testMap.find(2*i);
        auto it2 = testMap.find(2*i + 1);
        BOOST_TEST((it != testMap.end()));
        BOOST_TEST((it2 == testMap.end()));
        BOOST_TEST(it->second.value == 2*i);
        BOOST_TEST(testMap.count(2*i) == 1);
        BOOST_TEST(testMap.count(2*i + 1) == 0);
    }
}

template<class Timestamp = bool>
void test_clear() {
    DeStruct::total() = 0;
    DeStruct::instances() = 0;
    constexpr size_t NUM_INSERT = 10;
    constexpr size_t NUM_CLEAR = ssize_t(std::numeric_limits<Timestamp>::max()) - ssize_t(std::numeric_limits<Timestamp>::min());
    pds::VecMap<size_t, DeStruct, Timestamp> testMap(NUM_INSERT);
    BOOST_TEST(testMap.empty());
    for (size_t i = 0; i < NUM_INSERT; ++i) {
        testMap[i] = i;
    }
    BOOST_TEST(DeStruct::instances() == NUM_INSERT);
    BOOST_TEST(testMap.size() == NUM_INSERT);
    BOOST_TEST(!testMap.empty());
    BOOST_TEST(testMap.capacity() >= NUM_INSERT);
    for (size_t i = 1; i < NUM_CLEAR; ++i) {
        testMap.clear();
        BOOST_TEST(DeStruct::instances() <= NUM_INSERT);
        BOOST_TEST(testMap.size() == 0);
        BOOST_TEST(testMap.empty());
        testMap.emplace(i % NUM_INSERT, 1337);
        BOOST_TEST(DeStruct::instances() <= NUM_INSERT);
        BOOST_TEST(testMap.capacity() >= NUM_INSERT);
        BOOST_TEST(testMap.size() == 1);
        BOOST_TEST(!testMap.empty());
    }
    BOOST_TEST(DeStruct::instances() != 0);
    BOOST_TEST(DeStruct::instances() <= NUM_INSERT);
    testMap.clear();
    BOOST_TEST(testMap.size() == 0);
    BOOST_TEST(testMap.empty());
    BOOST_TEST(DeStruct::instances() == 0);
}

BOOST_AUTO_TEST_CASE(clear) {
    BOOST_TEST_CONTEXT("bool timestamp") {
        test_clear<bool>();
    }
    BOOST_TEST(DeStruct::instances() == 0);
    BOOST_TEST_CONTEXT("uint8_t timestamp") {
        test_clear<uint8_t>();
    }
    BOOST_TEST(DeStruct::instances() == 0);
    BOOST_TEST_CONTEXT("char timestamp") {
        test_clear<char>();
    }
    BOOST_TEST(DeStruct::instances() == 0);
}
