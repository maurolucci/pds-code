//
// Created by max on 09.11.22.
//
#include <boost/test/unit_test.hpp>
#include "vecset.hpp"
#include <fmt/format.h>
#include <range/v3/all.hpp>

BOOST_AUTO_TEST_CASE(test_insert_remove_continous) {
    pds::VecSet<size_t> testSet;
    constexpr size_t NUM_INSERT = 100;
    for (size_t i = 0; i < NUM_INSERT; ++i) {
        BOOST_TEST(!testSet.contains(i));
        BOOST_TEST(testSet.count(i) == 0);
        testSet.insert(i);
        BOOST_TEST(testSet.contains(i));
        BOOST_TEST(testSet.count(i) == 1);
    }
    BOOST_TEST(testSet.size() == NUM_INSERT);
    BOOST_TEST(!testSet.empty());
    for (size_t i = 0; i < NUM_INSERT; ++i) {
        BOOST_TEST(testSet.contains(i));
        BOOST_TEST(testSet.count(i) == 1);
        testSet.erase(i);
        BOOST_TEST(!testSet.contains(i));
        BOOST_TEST(testSet.count(i) == 0);
    }
    BOOST_TEST(testSet.size() == 0);
    BOOST_TEST(testSet.empty());
    for (size_t i = 0; i < NUM_INSERT; ++i) {
        BOOST_TEST(!testSet.contains(i));
        testSet.erase(i);
        BOOST_TEST(testSet.size() == 0);
    }
}

void testIteration(auto& set) {
    BOOST_TEST(set.empty());
    BOOST_TEST((set.begin() == set.end()));
    constexpr size_t NUM_INSERT = 200;
    for (size_t i = 0; i < NUM_INSERT; ++i) {
        set.insert(2 * i);
        BOOST_TEST(set.contains(2*i));
        BOOST_TEST(set.size() == i + 1);
    }
    BOOST_TEST(set.size() == NUM_INSERT);
    auto it = set.begin();
    BOOST_TEST(*it == 0);
    size_t previous = 0;
    for (++it; it != set.end(); ++it) {
        BOOST_TEST(*it == previous + 2);
        previous = *it;
    }
    BOOST_TEST(previous == NUM_INSERT * 2 - 2);
    bool first = true;
    previous = -1;
    for (auto x: set) {
        if (!first) {
            BOOST_TEST(x == previous + 2);
        }
        previous = x;
    }
    BOOST_TEST(previous == NUM_INSERT * 2 - 2);
}

BOOST_AUTO_TEST_CASE(iteration) {
    BOOST_TEST_CONTEXT("timestamp=bool") {
        pds::VecSet<size_t, bool> testMap;
        testIteration(testMap);
        testMap.clear();
        testIteration(testMap);
        testMap.clear();
        testIteration(testMap);
    }
    BOOST_TEST_CONTEXT("timestamp=uint8_t") {
        pds::VecSet<size_t, std::uint8_t> testMap;
        BOOST_TEST_CONTEXT("fresh set") {
            testIteration(testMap);
        }
        testMap.clear();
        BOOST_TEST_CONTEXT("cleared") {
            testIteration(testMap);
        }
        testMap.clear();
        for (size_t i = 0; i < 254; ++i) {
            testMap.insert(i % 113);
            testMap.clear();
        }
        BOOST_TEST_CONTEXT("cleared until reset") {
            testIteration(testMap);
        }
    }
    BOOST_TEST_CONTEXT("timestamp=char") {
        pds::VecSet<size_t, signed char> testMap;
        BOOST_TEST_CONTEXT("fresh set") {
            testIteration(testMap);
        }
        testMap.clear();
        BOOST_TEST_CONTEXT("cleared") {
            testIteration(testMap);
        }
        testMap.clear();
        for (size_t i = 0; i < 254; ++i) {
            testMap.insert(i % 111);
            testMap.clear();
        }
        BOOST_TEST_CONTEXT("cleared until reset") {
            testIteration(testMap);
        }
    }
}
