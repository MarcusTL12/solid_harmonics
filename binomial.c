#include <stdint.h>
#include <stdio.h>

uint64_t binomial(uint64_t n, uint64_t k) {
    if (k > n) {
        return 0;
    }

    if (k > (n >> 1)) {
        k = n - k;
    }

    if (k == 0 || k == n) {
        return 1;
    }

    uint64_t x = n - k + 1, nn = x + 1, rr = 2;

    while (rr <= k) {
        x = (x * nn) / rr;
        rr += 1;
        nn += 1;
    }

    return x;
}
