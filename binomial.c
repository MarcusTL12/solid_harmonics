#include <stdint.h>

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
