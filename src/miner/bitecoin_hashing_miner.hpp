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

namespace bitecoin {

// This provides a primitive randomness step. It is not cryptographic quality,
// but suffices for these purposes. There is a constant c that comes from the
// server at the beginning of the round that gets used here.
void PoolHashMinerStep(bigint_t &x, const Packet_ServerBeginRound *pParams) {
  assert(NLIMBS == 4 * 2);

  bigint_t tmp;
  // tmp=lo(x)*c;
  wide_mul(4, tmp.limbs + 4, tmp.limbs, x.limbs, pParams->c);
  // [carry,lo(x)] = lo(tmp)+hi(x)
  uint32_t carry = wide_add(4, x.limbs, tmp.limbs, x.limbs + 4);
  // hi(x) = hi(tmp) + carry
  wide_add(4, x.limbs + 4, tmp.limbs + 4, carry);

  // overall:  tmp=lo(x)*c; x=tmp>hi(x)
}

// Given the various round parameters, this calculates the hash for a particular
// index value.
// Multiple hashes of different indices will be combined to produce the overall
// result.
bigint_t PoolHashMiner(const Packet_ServerBeginRound *pParams, uint32_t index) {
  assert(NLIMBS == 4 * 2);

  // Incorporate the existing block chain data - in a real system this is the
  // list of transactions we are signing. This is the FNV hash:
  // http://en.wikipedia.org/wiki/Fowler%E2%80%93Noll%E2%80%93Vo_hash_function
  hash::fnv<64> hasher;
  uint64_t chainHash =
      hasher((const char *)&pParams->chainData[0], pParams->chainData.size());

  // The value x is 8 words long (8*32 bits in total)
  // We build (MSB to LSB) as  [ chainHash ; roundSalt ; roundId ; index ]
  bigint_t x;
  wide_zero(8, x.limbs);
  wide_add(8, x.limbs, x.limbs, index); // chosen index goes in at two low limbs
  wide_add(6, x.limbs + 2, x.limbs + 2,
           pParams->roundId); // Round goes in at limbs 3 and 2
  wide_add(4, x.limbs + 4, x.limbs + 4,
           pParams->roundSalt); // Salt goes in at limbs 5 and 4
  wide_add(2, x.limbs + 6, x.limbs + 6,
           chainHash); // chainHash at limbs 7 and 6

  // Now step forward by the number specified by the server
  for (unsigned j = 0; j < pParams->hashSteps; j++) {
    PoolHashMinerStep(x, pParams);
  }
  return x;
}

// This is the complete hash reference function. Given the current round
// parameters,
// and the solution, which is a vector of indices, it calculates the proof. The
// goodness
// of the solution is determined by the numerical value of the proof.
bigint_t HashMiner(const Packet_ServerBeginRound *pParams,
                       unsigned nIndices, const uint32_t *pIndices) {
  if (nIndices > pParams->maxIndices)
    throw std::invalid_argument(
        "HashMiner - Too many indices for parameter set.");

  bigint_t acc;
  wide_zero(8, acc.limbs);

  for (unsigned i = 0; i < nIndices; i++) {
    if (i > 0) {
      if (pIndices[i - 1] >= pIndices[i])
        throw std::invalid_argument("HashMiner - Indices are not in "
                                    "monotonically increasing order.");
    }

    // Calculate the hash for this specific point
    bigint_t point = PoolHashMiner(pParams, pIndices[i]);

    // Combine the hashes of the points together using xor
    wide_xor(8, acc.limbs, acc.limbs, point.limbs);
  }

  return acc;
}
};

#endif
