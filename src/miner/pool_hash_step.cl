typedef struct {
  uint limbs[8];
} bigint_t;

void wide_zero(uint n, uint *res) {
  for (uint i = 0; i < n; i++) {
    res[i] = 0;
  }
}

uint wide_add(uint n, uint *res, const uint *a, const uint *b) {
  ulong carry = 0;
  for (uint i = 0; i < n; i++) {
    ulong tmp = convert_ulong(a[i]) + b[i] + carry;
    res[i] = convert_uint(tmp & 0xFFFFFFFFUL);
    carry = tmp >> 32;
  }
  return carry;
}

uint wide_add_carry(uint n, uint *res, const uint *a, uint b) {
  ulong carry = b;
  for (uint i = 0; i < n; i++) {
    ulong tmp = a[i] + carry;
    res[i] = convert_uint(tmp & 0xFFFFFFFFUL);
    carry = tmp >> 32;
  }
  return carry;
}

/*! Multiply two n-limb numbers to produce a 2n-limb result
        \note All the integers must be distinct, the output cannot overlap the
   input */
void wide_mul(uint n, uint *res_hi, uint *res_lo, const uint *a,
              __global const uint *b) {
  ulong carry = 0, acc = 0;
  for (uint i = 0; i < n; i++) {
    for (uint j = 0; j <= i; j++) {
      ulong tmp = convert_ulong(a[j]) * b[i - j];
      acc += tmp;
      if (acc < tmp)
        carry++;
      // fprintf(stderr, " (%d,%d)", j,i-j);
    }
    res_lo[i] = convert_uint(acc & 0xFFFFFFFFul);
    // fprintf(stderr, "\n  %d : %u\n", i, res_lo[i]);
    acc = (carry << 32) | (acc >> 32);
    carry = carry >> 32;
  }

  for (uint i = 1; i < n; i++) {
    for (uint j = i; j < n; j++) {
      ulong tmp = convert_ulong(a[j]) * b[n - j + i - 1];
      acc += tmp;
      if (acc < tmp)
        carry++;
    }
    res_hi[i - 1] = convert_uint(acc & 0xFFFFFFFFul);
    acc = (carry << 32) | (acc >> 32);
    carry = carry >> 32;
  }
  res_hi[n - 1] = acc;
}

// Kenrl to hash 2 pairs top half
__kernel void poolhash_pair_tophalf(__global const uint *indices,
                                    __global ulong *word1,
                                    __global ulong *word2, __global uint *c,
                                    uint roundId, uint roundSalt,
                                    uint chainHash, uint hashSteps,
                                    uint offset) {
  uint i = get_global_id(0);
  uint index = indices[i];
  bigint_t x1;
  x1.limbs[0] = index;
  x1.limbs[1] = 0;
  x1.limbs[2] = roundId;
  x1.limbs[3] = roundId;
  x1.limbs[4] = roundSalt;
  x1.limbs[5] = roundSalt;
  x1.limbs[6] = chainHash;
  x1.limbs[7] = chainHash;

  for (unsigned j = 0; j < hashSteps; j++) {
    bigint_t tmp;
    wide_zero(8, tmp.limbs);
    // tmp=lo(x)*c;
    wide_mul(4, tmp.limbs + 4, tmp.limbs, x1.limbs, c);
    // [carry,lo(x)] = lo(tmp)+hi(x)
    uint carry = wide_add(4, x1.limbs, tmp.limbs, x1.limbs + 4);
    // hi(x) = hi(tmp) + carry
    wide_add_carry(4, x1.limbs + 4, tmp.limbs + 4, carry);
  }

  bigint_t x2;
  x2.limbs[0] = index + offset;
  x2.limbs[1] = 0;
  x2.limbs[2] = roundId;
  x2.limbs[3] = roundId;
  x2.limbs[4] = roundSalt;
  x2.limbs[5] = roundSalt;
  x2.limbs[6] = chainHash;
  x2.limbs[7] = chainHash;

  for (unsigned j = 0; j < hashSteps; j++) {
    bigint_t tmp;
    wide_zero(8, tmp.limbs);
    // tmp=lo(x)*c;
    wide_mul(4, tmp.limbs + 4, tmp.limbs, x2.limbs, c);
    // [carry,lo(x)] = lo(tmp)+hi(x)
    uint carry = wide_add(4, x2.limbs, tmp.limbs, x2.limbs + 4);
    // hi(x) = hi(tmp) + carry
    wide_add_carry(4, x2.limbs + 4, tmp.limbs + 4, carry);
  }

  uint msw = x1.limbs[7] ^ x2.limbs[7];
  uint msw2 = x1.limbs[6] ^ x2.limbs[6];
  uint msw3 = x1.limbs[5] ^ x2.limbs[5];
  uint msw4 = x1.limbs[4] ^ x2.limbs[4];

  word1[i] = convert_ulong(msw2) | (convert_ulong(msw) << 32);
  word2[i] = convert_ulong(msw4) | (convert_ulong(msw3) << 32);
}
