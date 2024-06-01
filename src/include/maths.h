#pragma once

#include <limits>
#include <cmath>
#include <complex>
#include <array>

namespace scurvy::impl {
    using namespace std::complex_literals;

    constexpr double RELTOL = 1e-9;
    constexpr double ABSTOL = 1e-9;
    constexpr double RELTOL_DIST = 1e-6;
    constexpr double ABSTOL_DIST = 1e-6;
    constexpr auto NAN_D = std::numeric_limits<double>::quiet_NaN();
    constexpr auto NAN_CD = std::numeric_limits<std::complex<double>>::quiet_NaN();

    inline bool is_close(double x, double target, double reltol = RELTOL, double abstol = ABSTOL) {
        auto err = std::abs(x - target);

        if(target == 0) {
            return err <= abstol;
        }

        return err <= abstol || err / std::abs(target) <= reltol;
    }

    inline bool near_zero(double x) {
        return is_close(x, 0);
    }

    template<typename T>
    T diff_of_products(T a, T b, T c, T d) {
        auto w = d * c;
        auto e = std::fma(-d, c, w);
        auto f = std::fma(a, b, -w);
        return f + e;
    }
}
