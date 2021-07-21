#ifndef DIP_ABSTRACTMODEL_H
#define DIP_ABSTRACTMODEL_H

#include <memory>
#include <utility>
#include <NTL/ZZ.h>
#include <sstream>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/split_member.hpp>


#include "Options.h"

struct ProjectivePoint {
    friend class boost::serialization::access;
    /// This struct represents point in projective coordinate system P = (X : Y: Z)
    NTL::ZZ x{0};
    NTL::ZZ y{1};
    NTL::ZZ z{0};

    bool operator==(const ProjectivePoint &rhs) const {
        return x == rhs.x && y == rhs.y && z == rhs.z;
    }
    /// Help function for getting point coordinates
    [[nodiscard]] std::string get_str() const {
        std::ostringstream oss;
        oss << "(" << x << ", " << y << ", " << z << ")";
        return oss.str();
    }
    /// This function is used for serialization operation
    template<class Archive>
    void save(Archive &ar, const unsigned int version) const {
        std::ostringstream buffer;

        buffer << x;
        ar & buffer.str();
        buffer.str("");

        buffer << y;
        ar & buffer.str();
    }
    /// This function is used for loading serialized data
    template<class Archive>
    void load(Archive &ar, const unsigned int version) {
        std::string buffer;

        ar & buffer;
        NTL::conv(x, buffer.c_str());
        buffer.clear();

        ar & buffer;
        NTL::conv(y, buffer.c_str());
    }
    BOOST_SERIALIZATION_SPLIT_MEMBER()

};

class AbstractModel {
    /// Class represents abstract elliptic curve model
public:
    explicit AbstractModel(std::shared_ptr<Options> options, ProjectivePoint point) : _options(std::move(options)), INFINITY_POINT(std::move(point)) {
    }
    virtual ~AbstractModel() = default;
    /// This function represents adding two points P and Q
    [[nodiscard]] virtual ProjectivePoint add_points(const ProjectivePoint &P, const ProjectivePoint &Q) const = 0;
    /// This function represents doubling point P
    [[nodiscard]] virtual ProjectivePoint double_point(const ProjectivePoint &P) const = 0;
    /// This function multiplies point P with multiplier k.
    [[nodiscard]] virtual ProjectivePoint mul_points(const NTL::ZZ &k, const ProjectivePoint &P) const {
        // Double-and-add algorithm
        ProjectivePoint N = P, Q = INFINITY_POINT;
        auto multiplier = k;
        while (multiplier != 0) {
            if ((multiplier & 1) == 1) {
                Q = add_points(Q, N);
            }
            N = double_point(N);
            if (is_infinity_point(N)) {
                break;
            }
            multiplier >>= 1;
        }
        return Q;
    }
    /// This function generates new elliptic curve and returns point on this curve
    virtual ProjectivePoint generate_elliptic_curve() = 0;
    [[nodiscard]] virtual bool is_infinity_point(const ProjectivePoint &point) const noexcept {
        return point == INFINITY_POINT;
    }
    /// Abstract method for converting projective system to affine system
    [[nodiscard]] virtual NTL::ZZ try_get_factor(const ProjectivePoint &point) const noexcept = 0;

protected:
    /// Options from command-line
    std::shared_ptr<Options> _options;
    /// Stores value of infinity (neutral) point on elliptic curve
    ProjectivePoint INFINITY_POINT;
};


#endif //DIP_ABSTRACTMODEL_H
