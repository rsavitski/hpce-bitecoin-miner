===============================================================================
HPCE Coursework 6: bitecoin miner
Ryan Savitski & Yong Wen Chua
===============================================================================
General Idea:

Generating a pool of points (in parallel as they are independent, as long as all
indices are unique, which can be done in a cull check at a merging stage).
Followed by a search to put the most arithmetically close points together, which
is equivalent to doing the naive N^2 work and checking all pairs since the
arithmetic closeness will guarantee that xors will also have the most leading
zeroes.

The sorted pool is then walked over and each pair is examined (xored) and also
put into an array of metapoints (2+ points xored) for a later pass, if such is
allowed by maxindices. This metapoint array can then be similarly sorted and
iterated over, effectively clearing the next chunk of leading zeroes by
examining the optimal pairings of these new metapoints. This can be repeated
(effectively generating 2^x index metapoints) until maxindices is hit.  Note
that this iterative metapass approach also needs only one generation step of
intial points, where each subsequent metapass does at most as much work as the
previous. In other words, they all do NlogN sort + N serial scan, where N is the
initial point pool size, in fact, since duplicate indices cannot exist in the
final proof, the metapoint folding culls any such duplicates, meaning that each
successive metapoint pass can do less work than the previous. Additionally, the
arrays are more sorted and stable in subsequent passes, making this iterative
scheme both effective (by empirical evidence) and cheaper than scaling the size
of the initial point seeds.

Additionally, the implementation takes advantage of a differential
cryptographic weakness of poolhash to get around 32 leading zero clears for
free. The weakness is such that, given an index i and an index j=i+k (where k
is a constant), the arithmetic difference of their corresponding points after
poolhashing is the same as for any other pair of indices separated by k.
Therefore, if we find one pair of indices whose arithmetic difference (and thus
most xorings) has a lot of leading zeroes, we can then generate pairing of
points (metapoints) that are guaranteed to always have the same arithmetic
difference and thus amount of leading zeroes. Interestingly enough, these
offsets only depend on the prime constant c and numsteps of poolhash. Given
that c does not seem to be varied, it would be possible to construct a cache of
such offsets to save ~140 milliseconds for known rounds. As otherwise the code
needs to run an online attack against the offset (2^16 - 2^17 range as sorting
trick gives us an effective hashrate of N^2 and the 17th power is needed for a
different crypto property not explained here).

With these approaches combined, on laptop grade machines, we can generate
proofs with around 140 leading zeroes (average for 16 index rounds) in a couple
of seconds. For longer rounds (tens of seconds), there are diminishing returns
but up to 190 leading zeroes were cleared.

===============================================================================
OpenCL/TBB:

On a CPU-only implementation, the bottleneck is poolhashing for metapoint
generation. For profiled runs, around 80% of time was spent in wide_mul.

Therefore an OpenCL implementation of metapoint generation was written, using
inline ptx assembly to accelerate multiplication.

With a GPU-accelerated metapoint generation, the new bottleneck of the algorithm
is the NlogN sorting of metapoint arrays. This is accelerated through
parallel_sort from TBB, giving nearly linear further speedup with number of CPU
cores.
