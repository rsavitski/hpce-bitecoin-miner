#ifndef bitecoin_miner_endpoint_hpp
#define bitecoin_miner_endpoint_hpp

#define TBB

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
#include "tbb/parallel_sort.h"
#endif

#include "bitecoin_protocol.hpp"
#include "bitecoin_endpoint.hpp"
#include "bitecoin_endpoint_client.hpp"
#include "bitecoin_hashing_miner.hpp"
#include "cl_boilerplate.hpp"

#define __STDC_FORMAT_MACROS
#include <inttypes.h>
#include <cstdio>

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

class EndpointMiner : public EndpointClient {
private:
  EndpointMiner(EndpointMiner &) = delete;
  void operator=(const EndpointMiner &) = delete;

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

  unsigned pass2Size = 1 << 22;
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

  virtual void MakeBid(
      const std::shared_ptr<Packet_ServerBeginRound> roundInfo, // Information
                                                                // about this
                                                                // particular
                                                                // round
      const std::shared_ptr<Packet_ServerRequestBid> request,   // The specific
                                                                // request we
                                                                // received
      double period,       // How long this bidding period will last
      double skewEstimate, // An estimate of the time difference between us and
                           // the server (positive -> we are ahead)
      std::vector<uint32_t> &solution, // Our vector of indices describing the
                                       // solution
      uint32_t *pProof // Will contain the "proof", which is just the value
      ) {
    double tSafetyMargin =
        0.5;              // accounts for uncertainty in network conditions
    tSafetyMargin += 0.3; // from binning //TODO

    double tFinish =
        request->timeStampReceiveBids * 1e-9 + skewEstimate - tSafetyMargin;

    double t1, t2;

    Log(Log_Verbose, "MakeBid - start, total period=%lg.", period);

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

    struct point_top {
      uint64_t msdw;
      uint32_t indx;

      bool operator<(point_top const &other) const { return msdw < other.msdw; }
    };

    const unsigned ptvct_sz = 1 << 17;

    // point vector
    std::vector<point_top> pts;

    std::random_device rd;
    std::minstd_rand gen(rd());
    std::uniform_int_distribution<uint32_t> dis;

    t1 = now() * 1e-9;
    // generate point vector
    for (unsigned i = 0; i < ptvct_sz; ++i) {
      uint32_t id = dis(gen);
      point_top pt;

      bigint_t temphash = PoolHashMiner(roundInfo.get(), id, chainHash);

      pt.indx = id;
      pt.msdw =
          uint64_t(temphash.limbs[6]) | (uint64_t(temphash.limbs[7]) << 32);

      pts.push_back(pt);
    }

    //////////////////////////////////////////////////////////

    std::sort(pts.begin(), pts.end());

    //////////////////////////////////////////////////////////

    uint64_t best_diff = 0xFFFFFFFFFFFFFFFFULL;
    uint32_t best_offset = 0;

    // note: .end()-1 for proper edge case
    for (auto it = pts.begin(); it != pts.end() - 1; ++it) {
      uint64_t curr_dw = it->msdw;
      uint64_t next_dw = (it + 1)->msdw;
      uint32_t curr_indx = it->indx;
      uint32_t next_indx = (it + 1)->indx;

      if (curr_indx == next_indx) {
        Log(Log_Verbose, "[*] Skipped identical index sample in diff finder");
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
    Log(Log_Info, "Best diff: %016 PRIx64", best_diff);
    Log(Log_Info, "Best offset: %8x", best_offset);

    pts.clear();

    t2 = now() * 1e-9;
    Log(Log_Info, "<x> hashsteps: %u", roundInfo->hashSteps);
    Log(Log_Info, "[=] total diff_find : %lg", t2 - t1); // TODO

    // done with offset search

    //////////////////////////////////////////////////////////

    t1 = now() * 1e-9;

    struct metapoint {
      bigint_t point;
      std::vector<uint32_t> indices;

      bool operator<(metapoint const &other) const {
        return wide_compare(BIGINT_WORDS, point.limbs, other.point.limbs) < 0;
      }
    };

    // metapoint vector
    const unsigned metaptvct_sz =
        1 << 18; // TODO: autotune, maybe golden diff finder too

    std::vector<metapoint> metaN_fb; // front buffer
    std::vector<metapoint> metaN_bb; // back buffer

    std::uniform_int_distribution<uint32_t> dis2(0, 0xFFFFFFFE - best_offset);

    // generate metapoint vector
    for (unsigned i = 0; i < pass2Size; ++i) {
      pass2Index[i] = dis2(gen);
      // uint32_t id2 = id + best_offset;
      pass2Pairing[i] = i;
    }
    // OpenCL Time
    std::vector<cl::Event> copyEvents(2);
    queue.enqueueWriteBuffer(cBuffer, CL_FALSE, 0, 4 * sizeof(uint32_t),
                             roundInfo.get()->c, nullptr, &copyEvents[0]);
    queue.enqueueWriteBuffer(pass2Indices, CL_FALSE, 0,
                             pass2Size * sizeof(uint32_t), pass2Index, nullptr,
                             &copyEvents[1]);

    cl::NDRange offset(0); // Always start iterations at x=0, y=0
    cl::NDRange globalSize(
        pass2Size); // Global size must match the original loops
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
                            pass2Size * sizeof(uint32_t) * 8, pass2Hash,
                            &kernelExecution);
    //////////////////////////////////////////////////////////

    t2 = now() * 1e-9;
    Log(Log_Info, "[=] metapt gen : %lg", t2 - t1); // TODO

    t1 = now() * 1e-9;

    auto comparer = [&](const uint32_t &a, const uint32_t &b) {
      return wide_compare(8, pass2Hash + (a * 8), pass2Hash + (b * 8)) == -1;
    };

#ifdef TBB
    tbb::parallel_sort(pass2Pairing, pass2Pairing + pass2Size, comparer);
#else
    std::sort(pass2Pairing, pass2Pairing + pass2Size, comparer);
#endif
    t2 = now() * 1e-9;
    Log(Log_Info, "[=] metapt sort : %lg", t2 - t1); // TODO

    //////////////////////////////////////////////////////////

    bigint_t mbest;
    wide_ones(BIGINT_WORDS, mbest.limbs);
    std::vector<uint32_t> best_indices;

    t1 = now() * 1e-9;

    metapoint temp_meta;
    for (unsigned i = 0; i < pass2Size - 1; ++i) {
      std::vector<uint32_t> a(2), b(2);
      a[0] = pass2Index[pass2Pairing[i]];
      a[1] = a[0] + best_offset;

      b[0] = pass2Index[pass2Pairing[i + 1]];
      b[1] = b[0] + best_offset;
      if (merge_and_check_uniq(temp_meta.indices, a, b) == 0) {
        Log(Log_Verbose, "[*] Skipped identical idx in metapass");
        continue;
      }

      wide_xor(8, temp_meta.point.limbs, pass2Hash + (pass2Pairing[i] * 8),
               pass2Hash + (pass2Pairing[i + 1] * 8));
      metaN_fb.push_back(temp_meta);

      // bigint_t xoredw;
      // wide_xor(8, xoredw.limbs, it->point.limbs, (it + 1)->point.limbs);
      // TODO: retain for maxindices = 4 case
      // if (wide_compare(BIGINT_WORDS, xoredw.limbs, mbest.limbs) < 0) {
      //  mbest = xoredw;
      //  best_indices.clear();
      //  best_indices.insert(best_indices.begin(), it->indices.begin(),
      //                      it->indices.end());
      //  best_indices.insert(best_indices.begin(), (it + 1)->indices.begin(),
      //                      (it + 1)->indices.end());
      //}
    }
    // for (auto item : tempvct)
    //  Log(Log_Info, "thing: %u", item);
    // Log(Log_Verbose, "Stamp");  // TODO

    // 4 indices per metapt at this time

    t2 = now() * 1e-9;
    Log(Log_Info, "[=] metapt scan : %lg", t2 - t1); // TODO

    //////////////////////////////////////////////////////////

    Log(Log_Info, "Maxindices: %u", roundInfo->maxIndices);
    Log(Log_Info, "Hashsteps: %u", roundInfo->hashSteps);
    uint32_t maxidx = roundInfo->maxIndices;
    uint32_t lg2idx = 0;
    while (maxidx != 0) {
      maxidx >>= 1;
      lg2idx++;
    }
    lg2idx -= 1;
    Log(Log_Info, "Log2_index: %u", lg2idx);

    if (lg2idx <= 2) {
      throw std::runtime_error("TODO: less than 2 metapasses required"); // TODO
    }

    uint32_t metaN_passes = lg2idx - 2;
    Log(Log_Info, "{!} %u metameta passes", metaN_passes);

    //////////////////////////////////////////////////////////

    for (unsigned npass = 0; npass < metaN_passes; ++npass) {
      t1 = now() * 1e-9;

      std::sort(metaN_fb.begin(), metaN_fb.end());

      t2 = now() * 1e-9;
      Log(Log_Info, "[=] meta[%u] sort : %lg", npass, t2 - t1); // TODO
      t1 = now() * 1e-9;

      metaN_bb.clear(); // TAG

      for (auto it = metaN_fb.begin(); it != metaN_fb.end() - 1; ++it) {
        if (merge_and_check_uniq(temp_meta.indices, it->indices,
                                 (it + 1)->indices) == 0) {
          Log(Log_Verbose, "[*] Skipped identical idx in metapass");
          continue;
        }
        wide_xor(8, temp_meta.point.limbs, it->point.limbs,
                 (it + 1)->point.limbs);
        metaN_bb.push_back(temp_meta);

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
      Log(Log_Info, "Best indices .size(): %zu", best_indices.size());
      std::swap(metaN_bb, metaN_fb); // TAG
      t2 = now() * 1e-9;
      Log(Log_Info, "[=] meta[%u] scan : %lg", npass, t2 - t1); // TODO
    }

    //////////////////////////////////////////////////////////

    // t2 = now() * 1e-9;
    // Log(Log_Info, "[=] N extra passes : %lg", t2-t1);  // TODO

    Log(Log_Info, "Finished metasearch");

    // Log(Log_Info, "Salt: %" PRIx64 "", roundInfo->roundSalt);
    // Log(Log_Info, "c[0]: %x", roundInfo->c[0]);
    // Log(Log_Info, "c[1]: %x", roundInfo->c[1]);
    // Log(Log_Info, "c[2]: %x", roundInfo->c[2]);
    // Log(Log_Info, "c[3]: %x", roundInfo->c[3]);

    //////////////////////////////////////////////////////////

    //  double t = now() * 1e-9;  // Work out where we are against the deadline
    //  double timeBudget = tFinish - t;
    //  Log(Log_Debug, "Finish trial %d, time remaining =%lg seconds.", nTrials,
    //      timeBudget);

    std::vector<uint32_t> bestSolution(best_indices.begin(),
                                       best_indices.end());

    bigint_t proof;
    for (auto idx : bestSolution) {
      bigint_t hash = PoolHashMiner(roundInfo.get(), idx, chainHash);
      wide_xor(8, proof.limbs, proof.limbs, hash.limbs);
    }

    double score = wide_as_double(BIGINT_WORDS, proof.limbs);
    double leadingzeros = 256 - log(score) * 1.44269504088896340736;
    Log(Log_Fatal, "score=%lg, leading zeros=%lg.", score, leadingzeros);

    std::sort(bestSolution.begin(), bestSolution.end());

    solution = bestSolution;
    wide_copy(BIGINT_WORDS, pProof, proof.limbs);

    Log(Log_Verbose, "MakeBid - finish.");
  }
};

}; // bitecoin

#endif
