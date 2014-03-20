struct bigint_t {
  uint limbs[8];
};

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
              const uint *b) {
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
