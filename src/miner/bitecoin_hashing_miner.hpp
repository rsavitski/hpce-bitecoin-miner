#ifndef bitecoin_hashing_miner_hpp
#define bitecoin_hashing_miner_hpp

// Provies a basic (non-cryptographic) hash function
#include "contrib/fnv.hpp"

// Basic operations for dealing with multi-word integers.
// The individual words are often called "limbs", so a 256 bit
// integer contains 8 x 32-bit limbs.
#include "wide_int.h"

#include <cstdint>
#include <stdexcept>
#include <algorithm>

#include "bitecoin_protocol.hpp"
#include "bitecoin_hashing.hpp"

typedef __uint128_t uint128_t;

namespace bitecoin {

// http://locklessinc.com/articles/256bit_arithmetic/

uint128_t fastermul65(uint64_t x1, uint64_t x2, uint64_t y1, uint64_t y2,
                      uint64_t *carry = nullptr) {
  uint64_t c = 0;
  uint64_t t1 = x1 - x2;
  uint64_t t2 = y2 - y1;
  uint128_t result = (uint128_t)t1 * t2;

  if ((x2 > x1) && t2) {
    result -= (uint128_t)t2 << 64;
    c = ~c;
  }

  if ((y1 > y2) && t1) {
    result -= (uint128_t)t1 << 64;
    c = ~c;
  }

  if (carry) *carry = c;
  return result;
}

void fasterMul128(uint32_t *result, const uint32_t *x, const uint32_t *y) {
  uint64_t x1 = uint64_t(x[0]) | (uint64_t(x[1]) << 32);
  uint64_t x2 = uint64_t(x[2]) | (uint64_t(x[3]) << 32);
  uint64_t x3 = 0, x4 = 0;
  uint64_t y1 = uint64_t(y[0]) | (uint64_t(y[1]) << 32);
  uint64_t y2 = uint64_t(y[2]) | (uint64_t(y[3]) << 32);
  uint64_t y3 = 0, y4 = 0;

  uint128_t xy11 = (uint128_t) x1 * y1;
  uint128_t xy22 = (uint128_t) x2 * y2;
  uint128_t xy33 = (uint128_t) x3 * y3;
  uint64_t xy14 = x1 * y4 + x4 * y1;
  uint64_t xy23 = (x3 - x2)*(y2 - y3);

  uint64_t carry;

  uint128_t t1 = xy11;
  uint128_t t2 = fastermul65(x1, x2, y1, y2, &carry);
  uint128_t t3 = fastermul65(x1, x3, y1, y3) + xy11 + xy22 + xy33 + (t2 >> 64);
  uint64_t t4 = xy14 + xy23 + xy22 + xy33 + carry;

  t2 = (uint64_t) t2;
  t2 += t1 >> 64;
  t3 += t2 >> 64;

  t2 = (uint64_t) t2;
  t2 += xy11;
  t3 += t2 >> 64;

  t2 = (uint64_t) t2;
  t2 += xy22;
  t3 += t2 >> 64;

  result[0] = uint32_t(t1);
  result[1] = uint32_t(t1 >> 32);
  result[2] = uint32_t(t2);
  result[3] = uint32_t(t2 >> 32);

  uint128_t high = t3 + ((uint128_t) t4 << 64);
  result[4] = uint32_t(high);
  result[5] = uint32_t(high >> 32);
  result[6] = uint32_t(high >> 64);
  result[7] = uint32_t(high >> 96);
}

class fnvIterative {
  static const uint64_t FNV_64_PRIME = 0x100000001b3ULL;
  uint64_t m_offset;

  fnvIterative(const uint64_t init = INIT) : m_offset(init) {}
  fnvIterative(fnvIterative const &);
  void operator=(fnvIterative const &);

public:
  static fnvIterative &getInstance() {
    static fnvIterative instance;
    return instance;
  }
  static const uint64_t INIT = 0xcbf29ce484222325ULL;

  void reset() { offset(INIT); }
  void offset(uint64_t init = INIT) { m_offset = init; }

  uint64_t operator()(const std::string &_buf) {
    return operator()(_buf.c_str(), _buf.length());
  }

  uint64_t operator()(const char *_buf, size_t _len) {
    const unsigned char *bp =
        reinterpret_cast<const unsigned char *>(_buf); /* start of buffer */
    const unsigned char *be = bp + _len; /* beyond end of buffer */

    uint64_t hval = m_offset;

    /*
     * FNV-1a hash each octet of the buffer
     */
    while (bp < be) {

      /* xor the bottom with the current octet */
      hval ^= (uint64_t) * bp++;

/* multiply by the 64 bit FNV magic prime mod 2^64 */

#if defined(NO_FNV_GCC_OPTIMIZATION)
      hval *= FNV_64_PRIME;
#else
      hval += (hval << 1) + (hval << 4) + (hval << 5) + (hval << 7) +
              (hval << 8) + (hval << 40);
#endif
    }

    return m_offset = hval;
  }
};

// Given the various round parameters, this calculates the hash for a particular
// index value.
// Multiple hashes of different indices will be combined to produce the overall
// result.
bigint_t PoolHashMiner(const Packet_ServerBeginRound *pParams, uint32_t index,
                       uint64_t chainHash) {
  // assert(NLIMBS == 4 * 2);

  // The value x is 8 words long (8*32 bits in total)
  // We build (MSB to LSB) as  [ chainHash ; roundSalt ; roundId ; index ]
  bigint_t x;
  x.limbs[0] = index;
  x.limbs[2] = (uint32_t)(pParams->roundId & 0xFFFFFFFFULL);
  x.limbs[3] = (uint32_t)(pParams->roundId & 0xFFFFFFFFULL);
  x.limbs[4] = (uint32_t)(pParams->roundSalt & 0xFFFFFFFFULL);
  x.limbs[5] = (uint32_t)(pParams->roundSalt & 0xFFFFFFFFULL);
  x.limbs[6] = (uint32_t)(chainHash & 0xFFFFFFFFULL);
  x.limbs[7] = (uint32_t)(chainHash & 0xFFFFFFFFULL);

  // Now step forward by the number specified by the server
  for (unsigned j = 0; j < pParams->hashSteps; j++) {
    bigint_t tmp;
    // tmp=lo(x)*c;
    fasterMul128(tmp.limbs, x.limbs, pParams -> c);
    // [carry,lo(x)] = lo(tmp)+hi(x)
    uint32_t carry = wide_add(4, x.limbs, tmp.limbs, x.limbs + 4);
    // hi(x) = hi(tmp) + carry
    wide_add(4, x.limbs + 4, tmp.limbs + 4, carry);
  }
  return x;
}
};

#endif
