#ifndef UNIT_TESTS_M1_TEST_PRNG_H_
#define UNIT_TESTS_M1_TEST_PRNG_H_

#include <stdint.h>

typedef struct {
  uint64_t state;
} m1_test_rng;

/* SplitMix64-v1: state advance, finalizer, and high-bit conversion are part of
 * the deterministic replay contract. Unsigned overflow is defined by C99 for
 * uint64_t. */
static inline uint64_t m1_test_rng_next_u64(m1_test_rng *restrict rng) {
  uint64_t z = (rng->state += UINT64_C(0x9e3779b97f4a7c15));
  z = (z ^ (z >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
  z = (z ^ (z >> 27)) * UINT64_C(0x94d049bb133111eb);
  return z ^ (z >> 31);
}

static inline double m1_test_rng_unit(m1_test_rng *restrict rng) {
  return (double)(m1_test_rng_next_u64(rng) >> 11) * 0x1.0p-53;
}

static inline double
m1_test_rng_between(m1_test_rng *restrict rng, const double lower, const double upper) {
  return lower + (upper - lower) * m1_test_rng_unit(rng);
}

#endif
