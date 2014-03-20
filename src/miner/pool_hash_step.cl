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
              __constant uint *b) {
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

// Credits:
// https://devtalk.nvidia.com/default/topic/610914/cuda-programming-and-performance/modular-exponentiation-amp-biginteger/
// Inline PTX assembly

void add_asm(ulong *result, const ulong *a, const ulong *b) {
  asm("{\n\t"
      ".reg .u32 o0, o1, o2, o3, o4; \n\t"
      ".reg .u32 o5, o6, o7, i8, i9; \n\t"
      ".reg .u32 i10, i11, i12, i13; \n\t"
      ".reg .u32 i14, i15, i16, i17; \n\t"
      ".reg .u32 i18, i19, i20, i21; \n\t"
      ".reg .u32 i22, i23;           \n\t"
      "mov.b64         { i8, i9}, %4;\n\t"
      "mov.b64         {i10,i11}, %5;\n\t"
      "mov.b64         {i12,i13}, %6;\n\t"
      "mov.b64         {i14,i15}, %7;\n\t"
      "mov.b64         {i16,i17}, %8;\n\t"
      "mov.b64         {i18,i19}, %9;\n\t"
      "mov.b64         {i20,i21},%10;\n\t"
      "mov.b64         {i22,i23},%11;\n\t"
      "add.cc.u32      o0,  i8, i16; \n\t"
      "addc.cc.u32     o1,  i9, i17; \n\t"
      "addc.cc.u32     o2, i10, i18; \n\t"
      "addc.cc.u32     o3, i11, i19; \n\t"
      "addc.cc.u32     o4, i12, i20; \n\t"
      "addc.cc.u32     o5, i13, i21; \n\t"
      "addc.cc.u32     o6, i14, i22; \n\t"
      "addc.u32        o7, i15, i23; \n\t"
      "mov.b64         %0, {o0,o1};  \n\t"
      "mov.b64         %1, {o2,o3};  \n\t"
      "mov.b64         %2, {o4,o5};  \n\t"
      "mov.b64         %3, {o6,o7};  \n\t"
      "}"
      : "=l"(result[0]), "=l"(result[1]), "=l"(result[2]), "=l"(result[3])
      : "l"(a[0]), "l"(a[1]), "l"(a[2]), "l"(a[3]), "l"(b[0]), "l"(b[1]),
        "l"(b[2]), "l"(b[3]));
}

void mul_asm(ulong *result, const ulong *a, __constant ulong *b) {
  asm("{\n\t"
      ".reg .u32 o0, o1, o2, o3, o4;    \n\t"
      ".reg .u32 o5, o6, o7, i8, i9;    \n\t"
      ".reg .u32 i10, i11, i12, i13;    \n\t"
      ".reg .u32 i14, i15, i16, i17;    \n\t"
      ".reg .u32 i18, i19, i20, i21;    \n\t"
      ".reg .u32 i22, i23;              \n\t"
      "mov.b64         { i8, i9}, %4;   \n\t"
      "mov.b64         {i10,i11}, %5;   \n\t"
      "mov.b64         {i12,i13}, %6;   \n\t"
      "mov.b64         {i14,i15}, %7;   \n\t"
      "mov.b64         {i16,i17}, %8;   \n\t"
      "mov.b64         {i18,i19}, %9;   \n\t"
      "mov.b64         {i20,i21},%10;   \n\t"
      "mov.b64         {i22,i23},%11;   \n\t"
      "mul.lo.u32      o0,  i8, i16;    \n\t"
      "mul.hi.u32      o1,  i8, i16;    \n\t"
      "mad.lo.cc.u32   o1,  i8, i17, o1;\n\t"
      "madc.hi.u32     o2,  i8, i17,  0;\n\t"
      "mad.lo.cc.u32   o1,  i9, i16, o1;\n\t"
      "madc.hi.cc.u32  o2,  i9, i16, o2;\n\t"
      "madc.hi.u32     o3,  i8, i18,  0;\n\t"
      "mad.lo.cc.u32   o2,  i8, i18, o2;\n\t"
      "madc.hi.cc.u32  o3,  i9, i17, o3;\n\t"
      "madc.hi.u32     o4,  i8, i19,  0;\n\t"
      "mad.lo.cc.u32   o2,  i9, i17, o2;\n\t"
      "madc.hi.cc.u32  o3, i10, i16, o3;\n\t"
      "madc.hi.cc.u32  o4,  i9, i18, o4;\n\t"
      "madc.hi.u32     o5,  i8, i20,  0;\n\t"
      "mad.lo.cc.u32   o2, i10, i16, o2;\n\t"
      "madc.lo.cc.u32  o3,  i8, i19, o3;\n\t"
      "madc.hi.cc.u32  o4, i10, i17, o4;\n\t"
      "madc.hi.cc.u32  o5,  i9, i19, o5;\n\t"
      "madc.hi.u32     o6,  i8, i21,  0;\n\t"
      "mad.lo.cc.u32   o3,  i9, i18, o3;\n\t"
      "madc.hi.cc.u32  o4, i11, i16, o4;\n\t"
      "madc.hi.cc.u32  o5, i10, i18, o5;\n\t"
      "madc.hi.cc.u32  o6,  i9, i20, o6;\n\t"
      "madc.hi.u32     o7,  i8, i22,  0;\n\t"
      "mad.lo.cc.u32   o3, i10, i17, o3;\n\t"
      "madc.lo.cc.u32  o4,  i8, i20, o4;\n\t"
      "madc.hi.cc.u32  o5, i11, i17, o5;\n\t"
      "madc.hi.cc.u32  o6, i10, i19, o6;\n\t"
      "madc.hi.u32     o7,  i9, i21, o7;\n\t"
      "mad.lo.cc.u32   o3, i11, i16, o3;\n\t"
      "madc.lo.cc.u32  o4,  i9, i19, o4;\n\t"
      "madc.hi.cc.u32  o5, i12, i16, o5;\n\t"
      "madc.hi.cc.u32  o6, i11, i18, o6;\n\t"
      "madc.hi.u32     o7, i10, i20, o7;\n\t"
      "mad.lo.cc.u32   o4, i10, i18, o4;\n\t"
      "madc.lo.cc.u32  o5,  i8, i21, o5;\n\t"
      "madc.hi.cc.u32  o6, i12, i17, o6;\n\t"
      "madc.hi.u32     o7, i11, i19, o7;\n\t"
      "mad.lo.cc.u32   o4, i11, i17, o4;\n\t"
      "madc.lo.cc.u32  o5,  i9, i20, o5;\n\t"
      "madc.hi.cc.u32  o6, i13, i16, o6;\n\t"
      "madc.hi.u32     o7, i12, i18, o7;\n\t"
      "mad.lo.cc.u32   o4, i12, i16, o4;\n\t"
      "madc.lo.cc.u32  o5, i10, i19, o5;\n\t"
      "madc.lo.cc.u32  o6,  i8, i22, o6;\n\t"
      "madc.hi.u32     o7, i13, i17, o7;\n\t"
      "mad.lo.cc.u32   o5, i11, i18, o5;\n\t"
      "madc.lo.cc.u32  o6,  i9, i21, o6;\n\t"
      "madc.hi.u32     o7, i14, i16, o7;\n\t"
      "mad.lo.cc.u32   o5, i12, i17, o5;\n\t"
      "madc.lo.cc.u32  o6, i10, i20, o6;\n\t"
      "madc.lo.u32     o7,  i8, i23, o7;\n\t"
      "mad.lo.cc.u32   o5, i13, i16, o5;\n\t"
      "madc.lo.cc.u32  o6, i11, i19, o6;\n\t"
      "madc.lo.u32     o7,  i9, i22, o7;\n\t"
      "mad.lo.cc.u32   o6, i12, i18, o6;\n\t"
      "madc.lo.u32     o7, i10, i21, o7;\n\t"
      "mad.lo.cc.u32   o6, i13, i17, o6;\n\t"
      "madc.lo.u32     o7, i11, i20, o7;\n\t"
      "mad.lo.cc.u32   o6, i14, i16, o6;\n\t"
      "madc.lo.u32     o7, i12, i19, o7;\n\t"
      "mad.lo.u32      o7, i13, i18, o7;\n\t"
      "mad.lo.u32      o7, i14, i17, o7;\n\t"
      "mad.lo.u32      o7, i15, i16, o7;\n\t"
      "mov.b64         %0, {o0,o1};     \n\t"
      "mov.b64         %1, {o2,o3};     \n\t"
      "mov.b64         %2, {o4,o5};     \n\t"
      "mov.b64         %3, {o6,o7};     \n\t"
      "}"
      : "=l"(result[0]), "=l"(result[1]), "=l"(result[2]), "=l"(result[3])
      : "l"(a[0]), "l"(a[1]), "l"(a[2]), "l"(a[3]), "l"(b[0]), "l"(b[1]),
        "l"(b[2]), "l"(b[3]));
}
// Kenrl to hash 2 pairs top half
__kernel void poolhash_pair_tophalf(__global const uint *indices,
                                    __global ulong *word1,
                                    __global ulong *word2, __constant uint *c,
                                    const uint roundId, const uint roundSalt,
                                    const uint chainHash, const uint hashSteps,
                                    const uint offset) {
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
