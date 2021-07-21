#include <boost/test/unit_test.hpp>
#include "../src/Lenstra.h"
#include "../src/WeierstrassModel.h"
#include "../src/EdwardsModel.h"


struct TestFixture {
    TestFixture() {
        options = std::make_shared<Options>();
        options->weierstrass = true;
        *options->composite_number = NTL::ZZ(1000730021);
        weierstrass_model = std::make_shared<WeierstrassModel>(options);
        edwards_model = std::make_shared<EdwardsModel>(options);
    }
    std::shared_ptr<Options> options;
    std::shared_ptr<AbstractModel> weierstrass_model;
    std::shared_ptr<AbstractModel> edwards_model;
};

BOOST_FIXTURE_TEST_SUITE(lenstra_test, TestFixture)

BOOST_AUTO_TEST_CASE(test_weierstrass) {
    Lenstra test(options, weierstrass_model);
    auto result = test.factorize();
    BOOST_TEST((result == 100003 || result == 10007));
}

BOOST_AUTO_TEST_CASE(test_edwards) {
    Lenstra test(options, edwards_model);
    auto result = test.factorize();
    BOOST_TEST((result == 100003 || result == 10007));
}

BOOST_AUTO_TEST_SUITE_END()