#ifndef bitecoin_miner_endpoint_hpp
#define bitecoin_miner_endpoint_hpp

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
#include <functional>

#include "tbb/parallel_for.h"
#define __CL_ENABLE_EXCEPTIONS
#include "CL/cl.hpp"

#include "bitecoin_protocol.hpp"
#include "bitecoin_endpoint.hpp"
#include "bitecoin_endpoint_client.hpp"
#include "bitecoin_hashing_miner.hpp"
#include "cl_boilerplate.hpp"

#define __STDC_FORMAT_MACROS
#include <inttypes.h>
#include <cstdio>

namespace bitecoin {

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
  cl::Buffer pass2Word1, pass2Word2, pass2Indices;

  // Buffers: Approx 1 GB of memory usage now
  unsigned pass2Size = 1 << 24;
  uint64_t *pass2MSW;   // 2 most significant word (MSW followed by 2nd MSW)
  uint64_t *pass2TW;    // 3rd, 4th MS words
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
    pass2Word1 =
        cl::Buffer(context, CL_MEM_WRITE_ONLY, pass2Size * sizeof(uint64_t));
    pass2Word2 =
        cl::Buffer(context, CL_MEM_WRITE_ONLY, pass2Size * sizeof(uint64_t));
    pass2Indices =
        cl::Buffer(context, CL_MEM_READ_ONLY, pass2Size * sizeof(uint32_t));

    pass2MSW = new uint64_t[pass2Size]();
    pass2TW = new uint64_t[pass2Size]();
    pass2Index = new uint32_t[pass2Size]();
    pass2Pairing = new uint32_t[pass2Size]();
  }

  ~EndpointMiner() {
    delete[] pass2MSW;
    delete[] pass2TW;
    delete[] pass2Index;
    delete[] pass2Pairing;
  }

  EndpointMiner(const EndpointMiner &obj) = delete;

  /* Here is a default implementation of make bid.
          I would suggest that you override this method as a starting point.
  */
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
        0.5; // accounts for uncertainty in network conditions
    /* This is when the server has said all bids must be produced by, plus the
            adjustment for clock skew, and the safety margin
    */

    tSafetyMargin += 0.3; // from binning //TODO

    double tFinish =
        request->timeStampReceiveBids * 1e-9 + skewEstimate - tSafetyMargin;

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

    struct point_top {
      uint64_t msdw;
      uint32_t indx;

      bool operator<(point_top const &other) const { return msdw < other.msdw; }
    };

    const unsigned ptvct_sz = 1 << 24;

    // point vector
    std::vector<point_top> pts(ptvct_sz);

    std::random_device rd;
    std::minstd_rand gen(rd());
    std::uniform_int_distribution<uint32_t> dis;

    // generate point vector
    double t1 = now() * 1e-9;
    for (unsigned i = 0; i < ptvct_sz; ++i) {
      uint32_t id = dis(gen);
      point_top pt;
      pt.indx = id;
      pts[i] = pt;
    }
    // generate point vector
    tbb::parallel_for(0u, ptvct_sz, [&](unsigned i) {
      bigint_t temphash =
          PoolHashMiner(roundInfo.get(), pts[i].indx, chainHash);
      pts[i].msdw =
          uint64_t(temphash.limbs[6]) | (uint64_t(temphash.limbs[7]) << 32);

      // fprintf(stderr, "id: %x\n", id);
      // fprintf(stderr, "before: %8x %8x\n", temphash.limbs[7],
      // temphash.limbs[6]);
      // uint64_t sdh = uint64_t(temphash.limbs[6]) |
      // (uint64_t(temphash.limbs[7]) << 32);
      // fprintf(stderr, "after: %8x %8x\n\n", sdh>>32 , sdh & 0xffffFFFF);
    });
    double t2 = now() * 1e-9;
    Log(Log_Info, "Time taken %lf", t2 - t1);

    // for (auto pt : pts) {
    //  fprintf(stderr, "idx : %8x\n", pt.indx);
    //  fprintf(stderr, "msdw: %" PRIx64 "\n", pt.msdw);
    //}

    std::sort(pts.begin(), pts.end());

    // fprintf(stderr, "--------\n\n");
    // for (auto pt : pts) {
    //  fprintf(stderr, "idx : %8x\n", pt.indx);
    //  fprintf(stderr, "msdw: %" PRIx64 "\n", pt.msdw);
    //}
    // fprintf(stderr, "--------\n\n");

    uint64_t best_diff = 0xFFFFFFFFFFFFFFFFULL;
    uint32_t best_offset = 0;

    // note: .end()-1 for proper edge case
    for (auto it = pts.begin(); it != pts.end() - 1; ++it) {
      uint64_t curr_dw = it->msdw;
      uint64_t next_dw = (it + 1)->msdw;
      uint32_t curr_indx = it->indx;
      uint32_t next_indx = (it + 1)->indx;

      if (curr_indx == next_indx) {
        Log(Log_Verbose, "Skipped identical index sample in diff finder");
        continue;
      }

      uint64_t diff =
          (curr_dw > next_dw) ? curr_dw - next_dw : next_dw - curr_dw;

      // fprintf(stderr, "msdw_curr: %16" PRIx64 "\n", curr_dw);
      // fprintf(stderr, "msdw_next: %16" PRIx64 "\n", next_dw);
      // fprintf(stderr, "diff     : %16" PRIx64 "\n", diff);

      if (diff < best_diff) {
        best_diff = diff;
        best_offset = (curr_indx > next_indx) ? curr_indx - next_indx
                                              : next_indx - curr_indx;
        // fprintf(stderr, "FOUND better diff\n");
      }
      // fprintf(stderr, "\n");
    }
    // Log(Log_Fatal, "Best diff: %016" PRIx64 "", best_diff);
    Log(Log_Fatal, "Best offset: %8x", best_offset);

    struct metapoint_top {
      uint64_t msdw; // 2 most significant word (MSW followed by 2nd MSW)
      uint64_t tdw;  // 3rd, 4th MS words
      uint32_t indx; // base index

      bool operator<(metapoint_top const &other) const {
        if (msdw == other.msdw) {
          return tdw < other.tdw;
        } else
          return msdw < other.msdw;
      }
    };

    // metapoint vector
    std::vector<metapoint_top> metapts(pass2Size);

    std::uniform_int_distribution<uint32_t> dis2(0, 0xFFFFFFFE - best_offset);

    t1 = now() * 1e-9;
    // generate point vector
    for (unsigned i = 0; i < pass2Size; ++i) {
      pass2Index[i] = dis2(gen);
      pass2Pairing[i] = i;
    }

    // OpenCL Time
    queue.enqueueWriteBuffer(cBuffer, CL_TRUE, 0, 4 * sizeof(uint32_t),
                             roundInfo.get()->c);
    queue.enqueueWriteBuffer(pass2Indices, CL_TRUE, 0,
                             pass2Size * sizeof(uint32_t), pass2Index);

    // tbb::parallel_for(0u, pass2Size, [&](unsigned i) {
    //   bigint_t temphash =
    //       PoolHashMiner(roundInfo.get(), pass2Index[i], chainHash);
    //   bigint_t temphash2 = PoolHashMiner(
    //       roundInfo.get(), pass2Index[i] + best_offset, chainHash);

    //   uint32_t msw = temphash.limbs[7] ^ temphash2.limbs[7];
    //   uint32_t msw2 = temphash.limbs[6] ^ temphash2.limbs[6];
    //   uint32_t msw3 = temphash.limbs[5] ^ temphash2.limbs[5];
    //   uint32_t msw4 = temphash.limbs[4] ^ temphash2.limbs[4];

    //   pass2MSW[i] = uint64_t(msw2) | (uint64_t(msw) << 32);
    //   pass2TW[i] = uint64_t(msw4) | (uint64_t(msw3) << 32);
    // });
    t2 = now() * 1e-9;
    Log(Log_Info, "Time taken %lf", t2 - t1);
    // fprintf(stderr, "--------\n\n");
    // for (auto pt : metapts) {
    //  fprintf(stderr, "idx : %8x\n", pt.indx);
    //  fprintf(stderr, "msdw: %016" PRIx64 "\n", pt.msdw);
    //  fprintf(stderr, "tdw: %8x\n", pt.tdw);
    //}
    // fprintf(stderr, "--------\n\n");

    std::sort(pass2Pairing, pass2Pairing + pass2Size,
              [&](const uint32_t &a, const uint32_t &b) {
      if (pass2MSW[a] == pass2MSW[b]) {
        return pass2TW[a] < pass2TW[b];
      } else {
        return pass2MSW[a] < pass2MSW[b];
      }
    });

    // fprintf(stderr, "--------\n\n");
    // for (auto pt : metapts) {
    //  fprintf(stderr, "idx : %8x\n", pt.indx);
    //  fprintf(stderr, "msdw: %016" PRIx64 "\n", pt.msdw);
    //  fprintf(stderr, "tdw: %8x\n", pt.tdw);
    //}
    // fprintf(stderr, "--------\n\n");

    metapoint_top mbest_diff;
    mbest_diff.msdw = 0xFFFFFFFFFFFFFFFFULL;
    mbest_diff.tdw = 0xFFFFFFFF;

    uint32_t metaidx[2];

    for (unsigned i = 0; i < pass2Size - 1; ++i) {
      uint32_t curr_indx = pass2Index[pass2Pairing[i]];
      uint32_t next_indx = pass2Index[pass2Pairing[i + 1]];

      if (curr_indx == next_indx || curr_indx + best_offset == next_indx ||
          next_indx + best_offset == curr_indx) {
        Log(Log_Verbose,
            "Skipped identical index sample in second pass search");
        continue;
      }

      metapoint_top xored;
      xored.msdw = pass2MSW[pass2Pairing[i]] ^ pass2MSW[pass2Pairing[i + 1]];
      xored.tdw = pass2TW[pass2Pairing[i]] ^ pass2TW[pass2Pairing[i + 1]];

      if (xored < mbest_diff) {
        // fprintf(stderr, "FOUND better metadiff\n");

        mbest_diff = xored;
        metaidx[0] = curr_indx;
        metaidx[1] = next_indx;
      }
    }

    // fprintf(stderr, "--------\n\n");
    // fprintf(stderr, "idx : %8x\n", metaidx[0]);
    // fprintf(stderr, "idx : %8x\n", metaidx[1]);
    // fprintf(stderr, "msdw: %016" PRIx64 "\n", mbest_diff.msdw);
    // fprintf(stderr, "tdw: %8x\n", mbest_diff.tdw);
    // fprintf(stderr, "--------\n\n");
    Log(Log_Fatal, "Finished metasearch");

    ///////////////////////////////////////////

    unsigned nTrials = 0;
    uint32_t r = 0;
    while (1) {
      ++nTrials;

      Log(Log_Debug, "Trial %d.", nTrials);

      // TODO

      double t = now() * 1e-9; // Work out where we are against the deadline
      double timeBudget = tFinish - t;
      Log(Log_Debug, "Finish trial %d, time remaining =%lg seconds.", nTrials,
          timeBudget);

      if (timeBudget <= 0)
        break; // We have run out of time, send what we have
    }
    Log(Log_Info, "nTrials: %u", nTrials);
    // Log(Log_Info, "Effective hashrate: %lf",
    //    ((double)nTrials * BIN_SIZE * BIN_SIZE * BIN_SIZE) / (t2 - t1));

    std::vector<uint32_t> bestSolution;
    bestSolution.push_back(metaidx[0]);
    bestSolution.push_back(metaidx[0] + best_offset);
    bestSolution.push_back(metaidx[1]);
    bestSolution.push_back(metaidx[1] + best_offset);

    bigint_t proof;
    for (auto idx : bestSolution) {
      bigint_t hash = PoolHashMiner(roundInfo.get(), idx, chainHash);
      wide_xor(8, proof.limbs, proof.limbs, hash.limbs);
    }

    double score = wide_as_double(BIGINT_WORDS, proof.limbs);
    double leadingzeros = 256 - log(score) * 1.44269504088896340736;
    Log(Log_Info, "score=%lg, leading zeros=%lg.", score, leadingzeros);

    std::sort(bestSolution.begin(), bestSolution.end());

    solution = bestSolution;
    wide_copy(BIGINT_WORDS, pProof, proof.limbs);

    Log(Log_Verbose, "MakeBid - finish.");
  }
};

}; // bitecoin

#endif
