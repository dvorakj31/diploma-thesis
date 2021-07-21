#include "WeierstrassModel.h"

ProjectivePoint
WeierstrassModel::add_points(const ProjectivePoint &P, const ProjectivePoint &Q) const {
    /// This function uses formulas for calculation point addition on Weierstrass curve as described in thesis
    // if 0 + Q, return Q
    if (P == INFINITY_POINT) {
        return Q;
    }
    // if P + 0, return P
    if (Q == INFINITY_POINT) {
        return P;
    }

    NTL::ZZ A = Q.y * P.z;
    NTL::ZZ B = P.y * Q.z;
    NTL::ZZ C = Q.x * P.z;
    NTL::ZZ D = P.x * Q.z;
    NTL::ZZ E = A - B;
    NTL::ZZ F = C - D;
    NTL::ZZ G = F * F;
    NTL::ZZ H = G * F;
    NTL::ZZ I = P.z * Q.z;
    NTL::ZZ J = E * E * I - H - (G << 1) * D;
    return {
        (F * J) % *_ecc.modulus,
        (E * (G * D - J) - H * B) % *_ecc.modulus,
        (H * I) % *_ecc.modulus
    };
}

ProjectivePoint WeierstrassModel::double_point(const ProjectivePoint &P) const {
    /// This function uses formulas for calculation point doubling on Weierstrass curve as described in thesis
    if (P == INFINITY_POINT) {
        return P;
    }
    NTL::ZZ A = _ecc.a * P.z * P.z + P.x * P.x * 3;
    NTL::ZZ B = P.y * P.z;
    NTL::ZZ C = P.x * P.y * B;
    NTL::ZZ D = A * A - (C << 3);
    return {(B * D << 1) % *_ecc.modulus,
            (A * ((C << 2) - D) - P.y * P.y * (B << 3) * B) % *_ecc.modulus,
            ((B * B * B) << 3) % *_ecc.modulus
    };
}

ProjectivePoint WeierstrassModel::generate_elliptic_curve() {
    ProjectivePoint p;
    p.z = 1;
    if (!_ecc.modulus) {
        _ecc.modulus = _options->composite_number;
    }
    /// While duplicates or elliptic curve is singular try to generate new points and get compute elliptic curve parameters
    while (true) {
        p.x = NTL::RandomBnd(*_options->composite_number);
        p.y = NTL::RandomBnd(*_options->composite_number);
        _ecc.a = NTL::RandomBnd(*_options->composite_number);
        _ecc.b = (p.y * p.y - p.x * p.x * p.x - _ecc.a * p.x) % *_options->composite_number;
        if (_duplicates.find(_ecc) != _duplicates.end()) {
            continue;
        }
        if (_is_nonsingular(p)) {
            break;
        }
    }
    _duplicates.insert(_ecc);
    return p;
}

bool WeierstrassModel::_is_nonsingular(const ProjectivePoint &point) {
    /// Check if gcd(4A^3 + 27B^2, COMPOSITE_NUMBER) = 1
    return NTL::GCD(((_ecc.a * _ecc.a * _ecc.a) << 2) + _ecc.b * _ecc.b * 27, *_ecc.modulus) == 1;
}

NTL::ZZ WeierstrassModel::try_get_factor(const ProjectivePoint &point) const noexcept {
    /// Computes conversion from projective coordinates to affine
    return NTL::GCD(point.z, *_options->composite_number);
}


WeierstrassModel::EllipticCurve WeierstrassModel::get_elliptic_curve() const noexcept {
    return _ecc;
}

void WeierstrassModel::set_elliptic_curve(const WeierstrassModel::EllipticCurve &curve) noexcept {
    _ecc = curve;
    _ecc.modulus = _options->composite_number;
}
