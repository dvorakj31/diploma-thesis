#include <iostream>
#include <NTL/ZZ.h>
#include <memory>

#include <boost/program_options.hpp>

#include "Lenstra.h"
#include "Options.h"
#include "EdwardsModel.h"
#include "WeierstrassModel.h"

namespace po = boost::program_options;

int main(int argc, char **argv) {
    std::shared_ptr<Options> options = std::make_shared<Options>();

    po::options_description desc("OPTIONS");
    po::variables_map vm;
    desc.add_options()
            ("help,h", "produce help message")
            ("weierstrass_model,w", po::bool_switch(&options->weierstrass), "set Weierstrass model")
            ("edwards_model,e", po::bool_switch(&options->edwards), "set Edwards model")
            ("timer,t", po::bool_switch(&options->timer), "time measurement")
            ("parallel,p", po::bool_switch(&options->parallel), "start parallel")
            ("bound,b", po::value<NTL::ZZ>(options->bound.get()), "Maximal bound for iterations (Default square root of composite number)")
            ("composite-number,n", po::value<NTL::ZZ>(options->composite_number.get())->required(), "Positive integer bigger than 1 to factorize");
    try {
        po::store(po::parse_command_line(argc, argv, desc), vm);
        if (vm.count("help")) {
            std::cout << argv[0] << " [OPTIONS] --composite-number/-n COMPOSITE NUMBER\n";
            std::cout << desc << "\n";
            return 1;
        }
        po::notify(vm);
    } catch (std::exception &exception) {
        std::cout << exception.what() << "\n";
        return 1;
    }

    if (options->weierstrass && options->edwards) {
        std::cerr << "Only one model can be specified!\n";
        return 2;
    }

    options->weierstrass = !options->edwards;

    if (*options->composite_number < 2) {
        std::cerr << "Composite number must be positive integer bigger than 1!\n";
        return 3;
    }

    std::cout << "Factorizing number: " << *options->composite_number << '\n';
    std::cout << "Using model: " << (options->edwards ? "Edwards" : "Weierstrass") << '\n';
    std::cout << "Using timer: " << (options->timer ? "yes" : "no") << '\n';
    double start_time = 0.0, end_time;
    std::shared_ptr<AbstractModel> model;
    if (options->weierstrass) {
        model = std::make_shared<WeierstrassModel>(options);
    } else {
        model = std::make_shared<EdwardsModel>(options);
    }
    Lenstra ecm(options, model);
    if (options->timer) {
        start_time = NTL::GetTime();
    }
    NTL::ZZ factor;
    if (!options->parallel) {
        factor = ecm.factorize();
    } else {
        factor = ecm.factorize_parallel(argc, argv);
    }

    end_time = NTL::GetTime();
    if (options->timer && !options->parallel) {
        std::cout << "time = " << end_time - start_time << " s\n";
    }

    if (factor != 0)
        std::cout << "Factor = " << factor << "\n";

    return 0;
}
