#include "EdwardsModel.h"

ProjectivePoint EdwardsModel::add_points(const ProjectivePoint &P, const ProjectivePoint &Q) const {
    /// This function uses formulas for calculation point addition on Edwards curve as described in thesis
    if (P == INFINITY_POINT) {
        return Q;
    }

    if (Q == INFINITY_POINT) {
        return P;
    }

    if (P == Q) {
        return double_point(P);
    }

    NTL::ZZ A = P.z * Q.z;
    NTL::ZZ B = A * A;
    NTL::ZZ C = P.x * Q.x;
    NTL::ZZ D = P.y * Q.y;
    NTL::ZZ E = C * D;
    NTL::ZZ F = B - E;
    NTL::ZZ G = B + E;

    return {
            (A * F * ((P.x + P.y) * (Q.x + Q.y) - C - D)) % *_ecc.modulus,
            (A * G * (D - C)) % *_ecc.modulus,
            (F * G) % *_ecc.modulus
    };
}

ProjectivePoint EdwardsModel::double_point(const ProjectivePoint &P) const {
    /// This function uses formulas for calculation point addition on Edwards curve as described in thesis
    if (P == INFINITY_POINT) {
        return P;
    }

    NTL::ZZ B = (P.x + P.y) * (P.x + P.y);
    NTL::ZZ C = P.x * P.x;
    NTL::ZZ D = P.y * P.y;
    NTL::ZZ F = C + D;
    NTL::ZZ H = P.z * P.z;
    NTL::ZZ J = F - (H << 1);

    return {
            ((B - C - D) * J) % *_ecc.modulus,
            (F * (C - D)) % *_ecc.modulus,
            (F * J) % *_ecc.modulus,
    };
}

ProjectivePoint EdwardsModel::generate_elliptic_curve() {
    ProjectivePoint ret;
    ret.z = 1;
    _ecc.d = 1;
    _ecc.modulus = _options->composite_number;
    /// While d = 1 then generate new value of d.
    while (_ecc.d < 2) {
        ret.x = NTL::RandomBnd(*_ecc.modulus);
        ret.y = NTL::RandomBnd(*_ecc.modulus);
        auto square_x = NTL::PowerMod(ret.x, 2, *_ecc.modulus);
        auto square_y = NTL::PowerMod(ret.y, 2, *_ecc.modulus);
        auto mult = NTL::MulMod(square_x, square_y, *_ecc.modulus);
        if (NTL::GCD(mult, *_ecc.modulus) == 1) {
            auto inv = NTL::InvMod(mult, *_ecc.modulus);
            _ecc.d = NTL::MulMod(square_x + square_y - 1, inv, *_ecc.modulus);
        }
        /// Duplicate check
        if (_duplicates.find(_ecc) != _duplicates.end()) {
            _ecc.d = 1;
        }
    }

    _duplicates.insert(_ecc);

    return ret;
}

NTL::ZZ EdwardsModel::try_get_factor(const ProjectivePoint &point) const noexcept {
    /// Tries to convert (X : Y : Z) to (X/Z : Y/Z)
    return NTL::GCD(point.z, *_ecc.modulus);
}

EdwardsModel::EllipticCurve EdwardsModel::get_elliptic_curve() const noexcept {
    return _ecc;
}

void EdwardsModel::set_elliptic_curve(const EdwardsModel::EllipticCurve &curve) noexcept {
    _ecc = curve;
    /// Modulus is set only to composite number value from command-line
    _ecc.modulus = _options->composite_number;
}
