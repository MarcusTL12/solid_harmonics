#include <cmath>
#include <cstdint>
#include <vector>
#include <stdio.h>

struct gtoDef {
    int ncrt;
    int nsph;

    std::vector<std::vector<int>> cart_to_sph_indices;
    std::vector<std::vector<double>> crt_to_sph_coeffs;
};

uint64_t binomial(uint64_t n, uint64_t k) {
    if (k > n) {
        return 0;
    }

    if (k > (n >> 1)) {
        k = n - k;
    }

    if (k == 0) {
        return 1;
    }

    uint64_t x = n - k + 1;

    for (uint64_t nn = x + 1, rr = 2; rr <= k; rr++, nn++) {
        x = (x * nn) / rr;
    }

    return x;
}

double pow_by_squaring(double b, int e) {
    if (e < 0) {
        e = -e;
        b = 1.0 / b;
    }

    double p = b;
    double acc = 1.0;

    while (e != 0) {
        if (e & 1 != 0) {
            acc *= p;
        }
        p *= p;
        e >>= 1;
    }

    return acc;
}

double c_coeff(int l, int m, int t, int u, int v) {
    uint64_t c = binomial(l, t) * binomial(l - t, std::abs(m) + t) *
                 binomial(t, u) *
                 binomial(std::abs(m), 2 * v + (m < 0 ? 1 : 0));

    return static_cast<double>(c) * pow_by_squaring(4.0, -t);
}

double factorial(uint64_t n) {
    double acc = 1.0;

    for (uint64_t i = 2; i <= n; i++) {
        acc *= static_cast<double>(i);
    }

    return acc;
}

double ns_coeff(int l, int m) {
    double lf = factorial(l);
    double f1 = factorial(l + m) / lf;
    double f2 = factorial(l - m) / lf;
    double f3 = pow_by_squaring(2.0, -(2 * std::abs(m) + (m == 0 ? 1 : 0) - 1));

    return std::sqrt(f1 * f2 * f3);
}

int crt_to_ind(int lx, int ly, int lz) {
    return (lz * (2 * lx + 2 * ly + lz + 3) + 2 * ly) / 2;
}

gtoDef initialize_gto_coeffs(int l) {
    return {};
}

int main() {
    double x = c_coeff(0, 0, 0, 0, 0);

    printf("%.15f\n", x);

    return 0;
}
