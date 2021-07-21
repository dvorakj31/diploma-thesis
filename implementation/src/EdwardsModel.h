#ifndef DIP_EDWARDSMODEL_H
#define DIP_EDWARDSMODEL_H

#include "AbstractModel.h"

#include <set>
#include <sstream>
#include <boost/serialization/serialization.hpp>

class EdwardsModel final : public AbstractModel {
    public:

    struct EllipticCurve {
        /// This struct represents elliptic curve in Edwards form x^2 + y^2 = 1 + dx^2y^2
        friend class boost::serialization::access;

        NTL::ZZ d;

        std::shared_ptr<NTL::ZZ> modulus = nullptr;

        /// Auxiliary function for std::set, which requires this operator
        bool operator<(const EllipticCurve &rhs) const {
            return d < rhs.d;
        }
        /// This function is used for serialization operation
        template<class Archive>
        void save(Archive &ar, const unsigned int version) const {
            std::stringstream buffer;

            buffer << d;
            ar & buffer.str();
        }
        /// This function is used for loading serialized data
        template<class Archive>
        void load(Archive &ar, const unsigned int version) {
            std::string value;

            ar & value;
            NTL::conv(d, value.c_str());
        }
        BOOST_SERIALIZATION_SPLIT_MEMBER()
    };

    explicit EdwardsModel(std::shared_ptr<Options> options) : AbstractModel(std::move(options), {
                                                                        NTL::conv<NTL::ZZ>(0),
                                                                        NTL::conv<NTL::ZZ>(1),
                                                                        NTL::conv<NTL::ZZ>(1)})
    {}
    /// This function implements adding two points P and Q on Edwards curve
    [[nodiscard]] ProjectivePoint add_points(const ProjectivePoint &P, const ProjectivePoint &Q) const override;

    /// This function implements doubling point on Edwards curve
    [[nodiscard]] ProjectivePoint double_point(const ProjectivePoint &P) const override;

    /// This function implements generation of new elliptic curve and returns point on this curve
    ProjectivePoint generate_elliptic_curve() override;

    /// This function tries to find inversion of Z point. It can return divisor of modulus or another value (1 or modulus)
    [[nodiscard]] NTL::ZZ try_get_factor(const ProjectivePoint &point) const noexcept override;

    /// This function gets new elliptic curve
    [[nodiscard]] EllipticCurve get_elliptic_curve() const noexcept;

    /// This function sets new elliptic curve
    void set_elliptic_curve(const EllipticCurve &curve) noexcept;

private:

        EllipticCurve _ecc;

        std::set<EllipticCurve> _duplicates;
};


#endif //DIP_EDWARDSMODEL_H
