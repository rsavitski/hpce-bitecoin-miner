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
#include <gmp.h>

#include "bitecoin_protocol.hpp"
#include "bitecoin_hashing.hpp"

namespace bitecoin {

// wrapper for 128 bits integers
struct mpWrapper {
  mp_limb_t *limbs;
  const size_t words;
  mpWrapper(size_t size = 128) : words(size / mp_bits_per_limb) {
    assert(size % mp_bits_per_limb == 0);
    limbs = new mp_limb_t[words]();
  }

  ~mpWrapper() { delete[] limbs; }

  void importLimbs(const uint32_t *in) {
    size_t offset = 0;
    unsigned jEnd = mp_bits_per_limb / 32;
    for (unsigned i = 0; i < words; ++i) {
      mp_limb_t word = 0;
      for (unsigned j = 0; j < jEnd; ++j) {
        word += mp_limb_t(*(in + offset)) << (j * 32);
        offset++;
      }
      limbs[i] = word;
    }
  }

  void exportLimbs(uint32_t *out) {
    size_t offset = 0;
    unsigned jEnd = mp_bits_per_limb / 32;
    for (unsigned i = 0; i < words; ++i) {
      for (unsigned j = 0; j < jEnd; ++j) {
        out[offset] = uint32_t(limbs[i] >> j * 32);
        offset++;
      }
    }
  }

  void reset() {
    for (unsigned i = 0; i < words; ++i) {
      limbs[i] = 0;
    }
  }
};

void multiply128(mp_limb_t *results, mp_limb_t *a, mp_limb_t *b) {
  mpn_mul_n(results, a, b, 128 / mp_bits_per_limb);
}

mp_limb_t add128(mp_limb_t *results, mp_limb_t *a, mp_limb_t *b) {
  return mpn_add_n(results, a, b, 128 / mp_bits_per_limb);
}

// mp_limb_t add128(mpWrapper &results, mpWrapper &a, mp_limb_t b) {
//   mp_limb_t *_b = new mp_limb_t[128/mp_bits_per_limb]();
//   *_b = b;
//   mp_limb_t carry =  mpn_add_n(results.limbs, a.limbs, _b,
// 128/mp_bits_per_limb);
//   delete[] _b;
//   return carry;
// }

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

  mpWrapper x1, x2, carry, c;
  x1.importLimbs(x.limbs);
  x2.importLimbs(x.limbs + 4);
  c.importLimbs(pParams->c);
  mpWrapper temp(256);
  // Now step forward by the number specified by the server
  for (unsigned j = 0; j < pParams->hashSteps; j++) {
    temp.reset();
    // tmp=lo(x)*c;
    multiply128(temp.limbs, x1.limbs, c.limbs);
    // [carry,lo(x)] = lo(tmp)+hi(x)
    carry.limbs[0] = add128(x1.limbs, temp.limbs, x2.limbs);
    // hi(x) = hi(tmp) + carry
    add128(x2.limbs, temp.limbs + (128/mp_bits_per_limb), carry.limbs);
  }
  x1.exportLimbs(x.limbs);
  x2.exportLimbs(x.limbs + 4);
  return x;
}
};

#endif
