#ifndef DIP_OPTIONS_H
#define DIP_OPTIONS_H

#include <NTL/ZZ.h>

struct Options {
    /// This struct stores options from command-line
    std::shared_ptr<NTL::ZZ> composite_number = std::make_shared<NTL::ZZ>();
    bool edwards = false;
    bool weierstrass = false;
    std::shared_ptr<NTL::ZZ> bound = std::make_shared<NTL::ZZ>(0);
    bool timer = false;
    bool parallel = false;
};

#endif //DIP_OPTIONS_H
