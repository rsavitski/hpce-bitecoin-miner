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

// Credits: http://goo.gl/NtaADC
// Inline PTX assembly

// Not used because too many registers used -- register thrashing!
uint add_asm(uint *result, const uint *a, const uint *b) {
  uint carry;
  asm("{\n\t"
      "add.cc.u32      %0,  %9, %17; \n\t"
      "addc.cc.u32     %1,  %10, %18; \n\t"
      "addc.cc.u32     %2, %11, %19; \n\t"
      "addc.cc.u32     %3, %12, %20; \n\t"
      "addc.cc.u32     %4, %13, %21; \n\t"
      "addc.cc.u32     %5, %14, %22; \n\t"
      "addc.cc.u32     %6, %15, %23; \n\t"
      "addc.cc.u32     %7, %16, %24; \n\t"
      "addc.u32         %8, 0, 0; \n\t"
      "}"
      : "=r"(result[0]), "=r"(result[1]), "=r"(result[2]), "=r"(result[3]),
        "=r"(result[4]), "=r"(result[5]), "=r"(result[6]), "=r"(result[7]),
        "=r"(carry)
      : "r"(a[0]), "r"(a[1]), "r"(a[2]), "r"(a[3]), "r"(a[4]), "r"(a[5]),
        "r"(a[6]), "r"(a[7]), "r"(b[0]), "r"(b[1]), "r"(b[2]), "r"(b[3]),
        "r"(b[4]), "r"(b[5]), "r"(b[6]), "r"(b[7]));
}

// Not used because too many registers used -- register thrashing!
void add_asm_carry(uint *result, const uint *a, const uint b) {
  asm("{\n\t"
      "add.cc.u32      %0,  %8, %16; \n\t"
      "addc.cc.u32     %1,  %9, 0; \n\t"
      "addc.cc.u32     %2, %10, 0; \n\t"
      "addc.cc.u32     %3, %11, 0; \n\t"
      "addc.cc.u32     %4, %12, 0; \n\t"
      "addc.cc.u32     %5, %13, 0; \n\t"
      "addc.cc.u32     %6, %14, 0; \n\t"
      "addc.u32        %7, %15, 0; \n\t"
      "}"
      : "=r"(result[0]), "=r"(result[1]), "=r"(result[2]), "=r"(result[3]),
        "=r"(result[4]), "=r"(result[5]), "=r"(result[6]), "=r"(result[7])
      : "r"(a[0]), "r"(a[1]), "r"(a[2]), "r"(a[3]), "r"(a[4]), "r"(a[5]),
        "r"(a[6]), "r"(a[7]), "r"(b));
}

void mul_asm(uint *result, const uint *a, __constant uint *b) {
  asm("{\n\t"
      "mul.lo.u32      %0,  %8, %16;    \n\t"
      "mul.hi.u32      %1,  %8, %16;    \n\t"
      "mad.lo.cc.u32   %1,  %8, %17, %1;\n\t"
      "madc.hi.u32     %2,  %8, %17,  0;\n\t"
      "mad.lo.cc.u32   %1,  %9, %16, %1;\n\t"
      "madc.hi.cc.u32  %2,  %9, %16, %2;\n\t"
      "madc.hi.u32     %3,  %8, %18,  0;\n\t"
      "mad.lo.cc.u32   %2,  %8, %18, %2;\n\t"
      "madc.hi.cc.u32  %3,  %9, %17, %3;\n\t"
      "madc.hi.u32     %4,  %8, %19,  0;\n\t"
      "mad.lo.cc.u32   %2,  %9, %17, %2;\n\t"
      "madc.hi.cc.u32  %3, %10, %16, %3;\n\t"
      "madc.hi.cc.u32  %4,  %9, %18, %4;\n\t"
      "madc.hi.u32     %5,  %8, %20,  0;\n\t"
      "mad.lo.cc.u32   %2, %10, %16, %2;\n\t"
      "madc.lo.cc.u32  %3,  %8, %19, %3;\n\t"
      "madc.hi.cc.u32  %4, %10, %17, %4;\n\t"
      "madc.hi.cc.u32  %5,  %9, %19, %5;\n\t"
      "madc.hi.u32     %6,  %8, %21,  0;\n\t"
      "mad.lo.cc.u32   %3,  %9, %18, %3;\n\t"
      "madc.hi.cc.u32  %4, %11, %16, %4;\n\t"
      "madc.hi.cc.u32  %5, %10, %18, %5;\n\t"
      "madc.hi.cc.u32  %6,  %9, %20, %6;\n\t"
      "madc.hi.u32     %7,  %8, %23,  0;\n\t"
      "mad.lo.cc.u32   %3, %10, %17, %3;\n\t"
      "madc.lo.cc.u32  %4,  %8, %20, %4;\n\t"
      "madc.hi.cc.u32  %5, %11, %17, %5;\n\t"
      "madc.hi.cc.u32  %6, %10, %19, %6;\n\t"
      "madc.hi.u32     %7,  %9, %21, %7;\n\t"
      "mad.lo.cc.u32   %3, %11, %16, %3;\n\t"
      "madc.lo.cc.u32  %4,  %9, %19, %4;\n\t"
      "madc.hi.cc.u32  %5, %12, %16, %5;\n\t"
      "madc.hi.cc.u32  %6, %11, %18, %6;\n\t"
      "madc.hi.u32     %7, %10, %20, %7;\n\t"
      "mad.lo.cc.u32   %4, %10, %18, %4;\n\t"
      "madc.lo.cc.u32  %5,  %8, %21, %5;\n\t"
      "madc.hi.cc.u32  %6, %12, %17, %6;\n\t"
      "madc.hi.u32     %7, %11, %19, %7;\n\t"
      "mad.lo.cc.u32   %4, %11, %17, %4;\n\t"
      "madc.lo.cc.u32  %5,  %9, %20, %5;\n\t"
      "madc.hi.cc.u32  %6, %13, %16, %6;\n\t"
      "madc.hi.u32     %7, %12, %18, %7;\n\t"
      "mad.lo.cc.u32   %4, %12, %16, %4;\n\t"
      "madc.lo.cc.u32  %5, %10, %19, %5;\n\t"
      "madc.lo.cc.u32  %6,  %8, %23, %6;\n\t"
      "madc.hi.u32     %7, %13, %17, %7;\n\t"
      "mad.lo.cc.u32   %5, %11, %18, %5;\n\t"
      "madc.lo.cc.u32  %6,  %9, %21, %6;\n\t"
      "madc.hi.u32     %7, %14, %16, %7;\n\t"
      "mad.lo.cc.u32   %5, %12, %17, %5;\n\t"
      "madc.lo.cc.u32  %6, %10, %20, %6;\n\t"
      "madc.lo.u32     %7,  %8, %23, %7;\n\t"
      "mad.lo.cc.u32   %5, %13, %16, %5;\n\t"
      "madc.lo.cc.u32  %6, %11, %19, %6;\n\t"
      "madc.lo.u32     %7,  %9, %23, %7;\n\t"
      "mad.lo.cc.u32   %6, %12, %18, %6;\n\t"
      "madc.lo.u32     %7, %10, %21, %7;\n\t"
      "mad.lo.cc.u32   %6, %13, %17, %6;\n\t"
      "madc.lo.u32     %7, %11, %20, %7;\n\t"
      "mad.lo.cc.u32   %6, %14, %16, %6;\n\t"
      "madc.lo.u32     %7, %12, %19, %7;\n\t"
      "mad.lo.u32      %7, %13, %18, %7;\n\t"
      "mad.lo.u32      %7, %14, %17, %7;\n\t"
      "mad.lo.u32      %7, %15, %16, %7;\n\t"
      "}"
      : "=r"(result[0]), "=r"(result[1]), "=r"(result[2]), "=r"(result[3]),
        "=r"(result[4]), "=r"(result[5]), "=r"(result[6]), "=r"(result[7])
      : "r"(a[0]), "r"(a[1]), "r"(a[2]), "r"(a[3]), "r"(a[4]), "r"(a[5]),
        "r"(a[6]), "r"(a[7]), "r"(b[0]), "r"(b[1]), "r"(b[2]), "r"(b[3]),
        "r"(b[4]), "r"(b[5]), "r"(b[6]), "r"(b[7]));
}
// kernel to hash 2 pairs top half
__kernel void poolhash_pair_tophalf(__global const uint *indices,
                                    __global uint *hashes, __constant uint *c,
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
    // mul_asm(tmp.limbs, x1.limbs, c);
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
    // mul_asm(tmp.limbs, x2.limbs, c);
    wide_mul(4, tmp.limbs + 4, tmp.limbs, x2.limbs, c);
    // [carry,lo(x)] = lo(tmp)+hi(x)
    uint carry = wide_add(4, x2.limbs, tmp.limbs, x2.limbs + 4);
    // hi(x) = hi(tmp) + carry
    wide_add_carry(4, x2.limbs + 4, tmp.limbs + 4, carry);
  }

  hashes[i*8] = x1.limbs[0] ^ x2.limbs[0];
  hashes[i*8 + 1] = x1.limbs[1] ^ x2.limbs[1];
  hashes[i*8 + 2] = x1.limbs[2] ^ x2.limbs[2];
  hashes[i*8 + 3] = x1.limbs[3] ^ x2.limbs[3];
  hashes[i*8 + 4] = x1.limbs[4] ^ x2.limbs[4];
  hashes[i*8 + 5] = x1.limbs[5] ^ x2.limbs[5];
  hashes[i*8 + 6] = x1.limbs[6] ^ x2.limbs[6];
  hashes[i*8 + 7] = x1.limbs[7] ^ x2.limbs[7];
}
