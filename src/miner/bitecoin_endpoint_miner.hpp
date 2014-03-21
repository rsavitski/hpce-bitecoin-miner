#ifndef bitecoin_miner_endpoint_hpp
#define bitecoin_miner_endpoint_hpp

#define TBB

#define __STDC_FORMAT_MACROS
#include <inttypes.h>
#include <cstdio>

#include <cstdarg>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <cassert>

#include <random>
#include <vector>
#include <memory>
#include <map>
#include <algorithm>
#include <set>
#include <utility>
#include <functional>

#define __CL_ENABLE_EXCEPTIONS
#include "CL/cl.hpp"

#ifdef TBB
#include "tbb/parallel_for.h"
#include "tbb/parallel_sort.h"
#endif

#include "bitecoin_protocol.hpp"
#include "bitecoin_endpoint.hpp"
#include "bitecoin_endpoint_client.hpp"
#include "bitecoin_hashing_miner.hpp"
#include "cl_boilerplate.hpp"

namespace bitecoin {

// merge step of two sorted vectors into a third with additional check of
// uniqueness of all values (important for point fold uniqueness checks)
// returns 0 if identical elements were found
int merge_and_check_uniq(std::vector<uint32_t> &res, std::vector<uint32_t> &v1,
                         std::vector<uint32_t> &v2) {
  res.clear();
  auto it1 = v1.begin();
  auto it2 = v2.begin();

  // run of the mill merge loop
  while (it1 != v1.end() && it2 != v2.end()) {
    if (*it1 < *it2) {
      res.push_back(*it1++);
    } else if (*it2 < *it1) {
      res.push_back(*it2++);
    } else
      return 0; // found identical element, will skip pair
  }
  while (it1 != v1.end())
    res.push_back(*it1++);
  while (it2 != v2.end())
    res.push_back(*it2++);

  return 1;
};

//////////////////////////////////////////////////////////

// glue to finaliseBid
struct timing_data{
  double tdiff_find;
  double tdiff_per_numstep;
  uint64_t work_sz;
  double tmeta_start;
  double hashrate;
  unsigned ismeta;
  unsigned metafactor;
};

//////////////////////////////////////////////////////////

// last thing called by makebid before quitting
// copies proof and indices to caller
void finaliseBid(timing_data &tdata, std::vector<uint32_t> best_indices, uint64_t chainHash,
    const std::shared_ptr<Packet_ServerBeginRound> roundInfo,
    std::vector<uint32_t> &solution, uint32_t *pProof);

// trivial case of maxindices == 1. Scrappy implementation for correctness
void direct_idx1_search(uint64_t work_size, bigint_t &mbest,
                   std::vector<uint32_t> &best_indices, std::minstd_rand &gen,
                   std::uniform_int_distribution<uint32_t> &dis,
                   const std::shared_ptr<Packet_ServerBeginRound> roundInfo,
                   uint64_t chainHash);

// maxindices == 2 case, simple linear scan as we have already found
// the optimal offset, just permuting carry bit noise
void idx2_scan(uint64_t work_size, uint32_t best_offset, bigint_t &mbest,
               std::vector<uint32_t> &best_indices, std::minstd_rand &gen,
               std::uniform_int_distribution<uint32_t> &dis2,
               const std::shared_ptr<Packet_ServerBeginRound> roundInfo,
               uint64_t chainHash);

//////////////////////////////////////////////////////////

class EndpointMiner : public EndpointClient {
private:
  unsigned m_knownRounds;
  std::map<std::string, unsigned> m_knownCoins;

  // OpenCL Context
  std::vector<cl::Platform> platforms;
  std::vector<cl::Device> devices;
  cl::Device device;
  cl::Context context;
  cl::Program program;
  cl::CommandQueue queue;

  cl::Buffer cBuffer;

  cl::Kernel pass2Kernel;
  cl::Buffer pass2Hashes, pass2Indices;

  unsigned pass2Size = 1 << 25; // max work size
  uint32_t *pass2Hash;
  uint32_t *pass2Index; // base index
  uint32_t *pass2Pairing;

public:
  EndpointMiner(std::string clientId, std::string minerId,
                std::unique_ptr<Connection> &conn, std::shared_ptr<ILog> &log)
      : EndpointClient(clientId, minerId, conn, log) {

    program = setupOpenCL(platforms, devices, device, context, log);
    queue = cl::CommandQueue(context, device);

    cBuffer = cl::Buffer(context, CL_MEM_READ_ONLY, 4 * sizeof(uint32_t));

    pass2Kernel = cl::Kernel(program, "poolhash_pair_tophalf");
    pass2Hashes = cl::Buffer(context, CL_MEM_WRITE_ONLY,
                             pass2Size * sizeof(uint32_t) * 8);
    pass2Indices =
        cl::Buffer(context, CL_MEM_READ_ONLY, pass2Size * sizeof(uint32_t));

    pass2Kernel.setArg(0, pass2Indices);
    pass2Kernel.setArg(1, pass2Hashes);
    pass2Kernel.setArg(2, cBuffer);

    pass2Hash = new uint32_t[pass2Size * 8]();
    pass2Index = new uint32_t[pass2Size]();
    pass2Pairing = new uint32_t[pass2Size]();
  }

  ~EndpointMiner() {
    delete[] pass2Hash;
    delete[] pass2Index;
    delete[] pass2Pairing;
  }
  EndpointMiner(const EndpointMiner &) = delete;
  void operator=(const EndpointMiner &) = delete;

  virtual void MakeBid(
      const std::shared_ptr<Packet_ServerBeginRound> roundInfo,
      const std::shared_ptr<Packet_ServerRequestBid> request,
      double period,       // How long this bidding period will last
      double skewEstimate, // An estimate of the time difference between us and
                           // the server (positive -> we are ahead)
      std::vector<uint32_t> &solution, // solution -> index array
      uint32_t *pProof                 // proof
      ) {

    //////////////////////////////////////////////////////////

    Log(Log_Info, "Maxindices: %u", roundInfo->maxIndices);

    uint32_t maxidx = roundInfo->maxIndices;
    uint32_t lg2idx = 0;
    while (maxidx != 0) {
      maxidx >>= 1;
      lg2idx++;
    }
    lg2idx -= 1;


    Log(Log_Info, "Log2_index: %u", lg2idx);
    Log(Log_Info, "Hashsteps: %u", roundInfo->hashSteps);

    //////////////////////////////////////////////////////////

    // timing setup
    Log(Log_Info, "MakeBid - start, total period=%lg.", period);

    static timing_data tdata = { .tdiff_find = 0.1, .tdiff_per_numstep = 0.01, .work_sz = 0, .tmeta_start=0, .hashrate=(1<<18) };

    double tdiff_expected = tdata.tdiff_per_numstep*roundInfo->hashSteps;

    // accounts for uncertainty in network conditions
    double tSafetyMargin = 0.5;

    double tFinish =
        request->timeStampReceiveBids * 1e-9 + skewEstimate - tSafetyMargin;

    double timeframe = period - tSafetyMargin - tdiff_expected;
    fprintf(stderr, "[-] time frame: %lg\n",timeframe);

    fprintf(stderr, "[-] Expected time for tdiff_find: %lg\n", tdiff_expected);

    double t1, t2;

    //////////////////////////////////////////////////////////

    // TODO: incremental hasher temporarily disabled to work with local server
    // static size_t previousChainSize = 0;
    // uint64_t chainHash = fnvIterative::getInstance()(
    //    (const char *)&roundInfo->chainData[previousChainSize],
    //    roundInfo->chainData.size() - previousChainSize);
    // previousChainSize = roundInfo->chainData.size();

    hash::fnv<64> hasher;
    uint64_t chainHash = hasher((const char *)&roundInfo->chainData[0],
                                roundInfo->chainData.size());


    //////////////////////////////////////////////////////////

    bigint_t mbest; // best proof found
    wide_ones(BIGINT_WORDS, mbest.limbs);
    std::vector<uint32_t> best_indices; // indices of best proof

    std::random_device rd;
    std::minstd_rand gen(rd());
    std::uniform_int_distribution<uint32_t> dis;

    //////////////////////////////////////////////////////////

    // calculate work done in this round

    const unsigned ptvct_sz = 1 << 17; // initial point vector for diff search

    uint32_t exp_metaN_factor = (1+std::min(lg2idx - 2, 2u));
    tdata.metafactor = exp_metaN_factor;
    double ddt = tdata.hashrate * timeframe / double(exp_metaN_factor);

    unsigned metapass_sz;
    if (ddt < 0){
      lg2idx = 1; // fake short work
    }
    tdata.ismeta = 0; //hack

    metapass_sz = unsigned(ddt * 0.8); // safety factor
    metapass_sz = std::max(metapass_sz, (1u<<16));
    metapass_sz = std::min(metapass_sz, pass2Size);
    fprintf(stderr, "[!] running with size: %u\n", metapass_sz);

    tdata.work_sz = metapass_sz;

    //////////////////////////////////////////////////////////

    // trivial case of maxindices == 1
    if (lg2idx == 0) {
      direct_idx1_search((1 << 17), mbest, best_indices, gen, dis, roundInfo,
                         chainHash);
      finaliseBid(tdata ,best_indices, chainHash, roundInfo, solution, pProof);
      return;
    }

    //////////////////////////////////////////////////////////

    // top 2 words of the proof with the corresponding index
    // do not need entire 8 words for the initial offset search
    struct point_top {
      uint64_t msdw;
      uint32_t indx;

      bool operator<(point_top const &other) const { return msdw < other.msdw; }
    };

    // initial point vector for differential search
    std::vector<point_top> pts(ptvct_sz);

    t1 = now() * 1e-9;

#ifdef TBB
    for (unsigned i = 0; i < ptvct_sz; ++i) {
      pts[i].indx = dis(gen);
    }

    tbb::parallel_for(0u, ptvct_sz, [&](unsigned i) {
      bigint_t temphash = PoolHashMiner(roundInfo.get(),
                                        pts[i].indx, chainHash);
      pts[i].msdw =
          uint64_t(temphash.limbs[6]) | (uint64_t(temphash.limbs[7]) << 32);
    });
#else
    // generate point vector from a set of random indices
    for (unsigned i = 0; i < ptvct_sz; ++i) {
      pts[i].indx = dis(gen);
      bigint_t temphash = PoolHashMiner(roundInfo.get(),
                                        pts[i].indx, chainHash);
      pts[i].msdw =
          uint64_t(temphash.limbs[6]) | (uint64_t(temphash.limbs[7]) << 32);
    }
#endif

    //////////////////////////////////////////////////////////
#ifdef TBB
    tbb::parallel_sort(pts.begin(), pts.end());
#else
    // sort to later consider closest pairs
    std::sort(pts.begin(), pts.end());
#endif

    //////////////////////////////////////////////////////////

    uint64_t best_diff = 0xFFFFFFFFFFFFFFFFULL;
    uint32_t best_offset = 0;

    // find best difference between generated points (most leading zeroes)
    // the offset between the two indices corresponding to this pair of points
    // is the best offset to use for generation of 4+ index proofs
    // (exploits a differential weakness in the hash)
    for (auto it = pts.begin(); it != pts.end() - 1; ++it) {
      uint64_t curr_dw = it->msdw;
      uint64_t next_dw = (it + 1)->msdw;
      uint32_t curr_indx = it->indx;
      uint32_t next_indx = (it + 1)->indx;

      if (curr_indx == next_indx) {
        // Log(Log_Verbose, "[*] Skipped identical index sample in diff
        // finder");
        continue;
      }

      // arithmetic difference for differential attack
      uint64_t diff =
          (curr_dw > next_dw) ? curr_dw - next_dw : next_dw - curr_dw;

      if (diff < best_diff) {
        best_diff = diff;
        best_offset = (curr_indx > next_indx) ? curr_indx - next_indx
                                              : next_indx - curr_indx;
      }
    }

    Log(Log_Info, "Best diff: %016" PRIx64 "", best_diff);
    Log(Log_Info, "Best offset: %8x", best_offset);

    pts.clear();

    t2 = now() * 1e-9;
    tdata.tdiff_find = t2-t1;
    Log(Log_Info, "[=] total diff_find : %lg", tdata.tdiff_find);

    //////////////////////////////////////////////////////////

    // random generator for valid base indices
    std::uniform_int_distribution<uint32_t> dis2(0, 0xFFFFFFFE - best_offset);

    //////////////////////////////////////////////////////////

    // done with offset search, do simple linear scan if maxindices == 2
    if (lg2idx == 1) {
      idx2_scan((1 << 17), best_offset, mbest, best_indices, gen, dis2,
                roundInfo, chainHash);
      finaliseBid(tdata, best_indices, chainHash, roundInfo, solution, pProof);
      return;
    }
    //////////////////////////////////////////////////////////

    tdata.ismeta = 1; // hack
    tdata.tmeta_start = now() * 1e-9;
    //////////////////////////////////////////////////////////

    t1 = now() * 1e-9;

    struct metapoint {
      bigint_t point;
      std::vector<uint32_t> indices;

      bool operator<(metapoint const &other) const {
        return wide_compare(BIGINT_WORDS, point.limbs, other.point.limbs) < 0;
      }
    };

    std::vector<metapoint> metaN_fb; // front buffer
    std::vector<metapoint> metaN_bb; // back buffer

    // generate metapoint vector
    for (unsigned i = 0; i < metapass_sz; ++i) {
      pass2Index[i] = dis2(gen);
      // corresponding id + best_offset generated on the gpu
      pass2Pairing[i] = i;
    }

    // generate metapoints
    std::vector<cl::Event> copyEvents(2);
    queue.enqueueWriteBuffer(cBuffer, CL_FALSE, 0, 4 * sizeof(uint32_t),
                             roundInfo.get()->c, nullptr, &copyEvents[0]);
    queue.enqueueWriteBuffer(pass2Indices, CL_FALSE, 0,
                             metapass_sz * sizeof(uint32_t), pass2Index,
                             nullptr, &copyEvents[1]);

    cl::NDRange offset(0); // Always start iterations at x=0, y=0
    cl::NDRange globalSize(
        metapass_sz); // Global size must match the original loops
    cl::NDRange localSize = cl::NullRange; // We don't care about local size

    pass2Kernel.setArg(3, cl_uint(roundInfo.get()->roundId));
    pass2Kernel.setArg(4, cl_uint(roundInfo.get()->roundSalt));
    pass2Kernel.setArg(5, cl_uint(chainHash));
    pass2Kernel.setArg(6, cl_uint(roundInfo.get()->hashSteps));
    pass2Kernel.setArg(7, cl_uint(best_offset));

    std::vector<cl::Event> kernelExecution(1);

    queue.enqueueNDRangeKernel(pass2Kernel, offset, globalSize, localSize,
                               &copyEvents, &kernelExecution[0]);
    queue.enqueueReadBuffer(pass2Hashes, CL_TRUE, 0,
                            metapass_sz * sizeof(uint32_t) * 8, pass2Hash,
                            &kernelExecution);

    //////////////////////////////////////////////////////////

    t2 = now() * 1e-9;
    Log(Log_Info, "[=] metapt gen : %lg", t2 - t1);
    t1 = now() * 1e-9;

    // sort positions in arrays of metapoints and corresponding indices by
    // metapoint value to put most similar values in adjacent bins
    auto comparer = [&](const uint32_t &a, const uint32_t &b) {
      return wide_compare(8, pass2Hash + (a * 8), pass2Hash + (b * 8)) == -1;
    };

#ifdef TBB
    tbb::parallel_sort(pass2Pairing, pass2Pairing + metapass_sz, comparer);
#else
    std::sort(pass2Pairing, pass2Pairing + metapass_sz, comparer);
#endif

    t2 = now() * 1e-9;
    Log(Log_Info, "[=] metapt sort : %lg", t2 - t1);

    //////////////////////////////////////////////////////////

    t1 = now() * 1e-9;

    // first meta pass, constructs 4-index proof tuples from pairs of metapoints
    // pushes these proofs into a vector for potential further passes
    metapoint temp_meta;
    for (unsigned i = 0; i < metapass_sz - 1; ++i) {
      std::vector<uint32_t> a(2), b(2);

      // reconstruct the original index pairs from metapoint generation
      a[0] = pass2Index[pass2Pairing[i]];
      a[1] = a[0] + best_offset;

      b[0] = pass2Index[pass2Pairing[i + 1]];
      b[1] = b[0] + best_offset;

      // check for duplicates in indices and merge into sorted array
      if (merge_and_check_uniq(temp_meta.indices, a, b) == 0) {
        // Log(Log_Verbose, "[*] Skipped identical idx in metapass");
        continue;
      }

      wide_xor(8, temp_meta.point.limbs, pass2Hash + (pass2Pairing[i] * 8),
               pass2Hash + (pass2Pairing[i + 1] * 8));

      // add to vector of metameta points for the potential next pass
      metaN_fb.push_back(temp_meta);

      // if new best proof, record value and indices
      if (wide_compare(BIGINT_WORDS, temp_meta.point.limbs, mbest.limbs) < 0) {
        mbest = temp_meta.point;
        best_indices.clear();
        best_indices.insert(best_indices.begin(), a.begin(), a.end());
        best_indices.insert(best_indices.begin(), b.begin(), b.end());
      }
    }

    Log(Log_Info, "Best indices .size(): %zu", best_indices.size());
    double TTTemp;
    TTTemp = wide_as_double(BIGINT_WORDS, mbest.limbs);
    Log(Log_Fatal, "(META ) score=%lg, leading zeros=%lg.", TTTemp,
        (256 - log(TTTemp) * 1.44269504088896340736));

    // 4 indices per metameta point for next stage

    t2 = now() * 1e-9;
    Log(Log_Info, "[=] metapt scan : %lg", t2 - t1);

    //////////////////////////////////////////////////////////

    // calculate how many more passes required (as it is always worthwhile to do
    // the extra work)
    // note: currently capped at 3 total metapasses
    uint32_t metaN_passes = std::min(lg2idx - 2, 2u);
    Log(Log_Info, "{!} %u metameta passes", metaN_passes);

    // if duplicate index culling managed to empty entire work
    // vector, we are done
    if (metaN_fb.size() < 2) {
      finaliseBid(tdata, best_indices, chainHash, roundInfo, solution, pProof);
      return;
    }

    //////////////////////////////////////////////////////////

    for (unsigned npass = 0; npass < metaN_passes; ++npass) {
      t1 = now() * 1e-9;
#ifdef TBB
      tbb::parallel_sort(metaN_fb.begin(), metaN_fb.end());
#else
      // sort to bring closest metapoints into adjacent positions
      std::sort(metaN_fb.begin(), metaN_fb.end());
#endif

      t2 = now() * 1e-9;
      Log(Log_Info, "[=] meta[%u] sort : %lg", npass, t2 - t1);
      t1 = now() * 1e-9;

      metaN_bb.clear();

      // walk over pairs of metapoints and see if their xor is the new best
      // proof
      // additionally, write these xors into back buffer for next metapass
      for (auto it = metaN_fb.begin(); it != metaN_fb.end() - 1; ++it) {

        // check for duplicates and merge two sorted index arrays
        // maintaining the sorted invariant
        if (merge_and_check_uniq(temp_meta.indices, it->indices,
                                 (it + 1)->indices) == 0) {
          // Log(Log_Verbose, "[*] Skipped identical idx in metapass");
          continue;
        }

        wide_xor(8, temp_meta.point.limbs, it->point.limbs,
                 (it + 1)->point.limbs);
        metaN_bb.push_back(temp_meta);

        // if new best proof, record value and indices
        if (wide_compare(BIGINT_WORDS, temp_meta.point.limbs, mbest.limbs) <
            0) {
          mbest = temp_meta.point;
          best_indices.clear();
          best_indices.insert(best_indices.begin(), it->indices.begin(),
                              it->indices.end());
          best_indices.insert(best_indices.begin(), (it + 1)->indices.begin(),
                              (it + 1)->indices.end());
        }
      }
      // swap buffers for next pass
      std::swap(metaN_bb, metaN_fb);

      // if duplicate index culling managed to empty entire work
      // vector, we are done
      if (metaN_fb.size() < 2) {
        finaliseBid(tdata, best_indices, chainHash, roundInfo, solution, pProof);
        return;
      }

      // logging
      Log(Log_Info, "Best indices .size(): %zu", best_indices.size());
      TTTemp = wide_as_double(BIGINT_WORDS, mbest.limbs);
      Log(Log_Fatal, "(META%d) score=%lg, leading zeros=%lg.", npass, TTTemp,
          (256 - log(TTTemp) * 1.44269504088896340736));
      t2 = now() * 1e-9;
      Log(Log_Info, "[=] meta[%u] scan : %lg", npass, t2 - t1);
    }

    //////////////////////////////////////////////////////////

    Log(Log_Info, "Finished metasearch");

    //////////////////////////////////////////////////////////

    finaliseBid(tdata, best_indices, chainHash, roundInfo, solution, pProof);
    return;
  }
};

// last thing called by makebid before quitting
// copies proof and indices to caller
void finaliseBid(timing_data &tdata,std::vector<uint32_t> best_indices, uint64_t chainHash,
                 const std::shared_ptr<Packet_ServerBeginRound> roundInfo,
                 std::vector<uint32_t> &solution, uint32_t *pProof) {

  // recalculate timing information for feedback
  tdata.tdiff_per_numstep = tdata.tdiff_find/roundInfo->hashSteps;
  fprintf(stderr, "[-] diff time per step: %lg\n", tdata.tdiff_per_numstep);

  if (tdata.ismeta){
    double tnow = now() * 1e-9;
    double ttaken = tnow - tdata.tmeta_start;
    fprintf(stderr, "[***] real time taken : %lg\n", ttaken);
    tdata.hashrate = (tdata.work_sz * double(tdata.metafactor)) / std::max(ttaken, 0.1);
    fprintf(stderr, "[***] new hashrate : %lg\n", tdata.hashrate);
  }

  // reconstruct best proof
  bigint_t proof;
  for (auto idx : best_indices) {
    bigint_t hash = PoolHashMiner(roundInfo.get(), idx, chainHash);
    wide_xor(8, proof.limbs, proof.limbs, hash.limbs);
  }

  // sort indices as they were generated out of order
  std::sort(best_indices.begin(), best_indices.end());

  // return results
  solution = best_indices;
  wide_copy(BIGINT_WORDS, pProof, proof.limbs);
}

// trivial case of maxindices == 1. Scrappy implementation for correctness
void
direct_idx1_search(uint64_t work_size, bigint_t &mbest,
                   std::vector<uint32_t> &best_indices, std::minstd_rand &gen,
                   std::uniform_int_distribution<uint32_t> &dis,
                   const std::shared_ptr<Packet_ServerBeginRound> roundInfo,
                   uint64_t chainHash) {
  for (unsigned i = 0; i < work_size; ++i) {
    uint32_t id = dis(gen);
    bigint_t temphash = PoolHashMiner(roundInfo.get(), id, chainHash);

    if (wide_compare(BIGINT_WORDS, temphash.limbs, mbest.limbs) < 0) {
      mbest = temphash;
      best_indices.clear();
      best_indices.push_back(id);
    }
  }
}

// maxindices == 2 case, simple linear scan as we have already found
// the optimal offset, just permuting carry bit noise
void idx2_scan(uint64_t work_size, uint32_t best_offset, bigint_t &mbest,
               std::vector<uint32_t> &best_indices, std::minstd_rand &gen,
               std::uniform_int_distribution<uint32_t> &dis2,
               const std::shared_ptr<Packet_ServerBeginRound> roundInfo,
               uint64_t chainHash) {
  for (unsigned i = 0; i < work_size; ++i) {
    uint32_t id = dis2(gen);
    uint32_t id2 = id + best_offset;

    bigint_t xoredw;
    bigint_t temphash = PoolHashMiner(roundInfo.get(), id, chainHash);
    bigint_t temphash2 = PoolHashMiner(roundInfo.get(), id2, chainHash);

    wide_xor(BIGINT_WORDS, xoredw.limbs, temphash.limbs, temphash2.limbs);

    if (wide_compare(BIGINT_WORDS, xoredw.limbs, mbest.limbs) < 0) {
      mbest = xoredw;
      best_indices.clear();
      best_indices.push_back(id);
      best_indices.push_back(id2);
    }
  }
}

}; // bitecoin

#endif
