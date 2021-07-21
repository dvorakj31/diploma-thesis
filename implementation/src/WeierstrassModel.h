#ifndef DIP_WEIERSTRASSMODEL_H
#define DIP_WEIERSTRASSMODEL_H

#include <set>
#include <utility>
#include <boost/serialization/string.hpp>

#include "AbstractModel.h"

class WeierstrassModel final : public AbstractModel {
public:
    struct EllipticCurve {
        friend class boost::serialization::access;

        // Elliptic curve in Weierstrass form y^2 = x^3 + ax + b
        NTL::ZZ a;
        NTL::ZZ b;

        std::shared_ptr<NTL::ZZ> modulus = nullptr;

        template<class Archive>
        void save(Archive &ar, const unsigned int version) const {
            std::stringstream buffer;

            buffer << a;
            ar & buffer.str();
            buffer.str("");

            buffer << b;
            ar & buffer.str();
        }

        template<class Archive>
        void load(Archive &ar, const unsigned int version) {
            std::string value;

            ar & value;
            NTL::conv(a, value.c_str());

            value.clear();
            ar & value;
            NTL::conv(b, value.c_str());
        }
        BOOST_SERIALIZATION_SPLIT_MEMBER()

        bool operator<(const EllipticCurve &rhs) const {
            return a < rhs.a || (a == rhs.a && b < rhs.b);
        }
    };

    explicit WeierstrassModel(std::shared_ptr<Options> options) : AbstractModel(std::move(options), {
                                                                NTL::conv<NTL::ZZ>(0),
                                                                NTL::conv<NTL::ZZ>(1),
                                                                NTL::conv<NTL::ZZ>(0)})
            {}

    WeierstrassModel(const WeierstrassModel &rhs) : AbstractModel(rhs._options, rhs.INFINITY_POINT), _ecc(rhs._ecc) {}

    [[nodiscard]] ProjectivePoint add_points(const ProjectivePoint &P, const ProjectivePoint &Q) const override;

    [[nodiscard]] ProjectivePoint double_point(const ProjectivePoint &P) const override;

    ProjectivePoint generate_elliptic_curve() override;

    [[nodiscard]] NTL::ZZ try_get_factor(const ProjectivePoint &point) const noexcept override;

    [[nodiscard]] EllipticCurve get_elliptic_curve() const noexcept;

    void set_elliptic_curve(const EllipticCurve &curve) noexcept;

private:

    EllipticCurve _ecc;
    std::set<EllipticCurve> _duplicates;

    bool _is_nonsingular(const ProjectivePoint &point);
};


#endif //DIP_WEIERSTRASSMODEL_H
