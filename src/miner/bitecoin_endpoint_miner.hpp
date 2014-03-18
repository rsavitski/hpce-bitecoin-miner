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

#include "bitecoin_protocol.hpp"
#include "bitecoin_endpoint.hpp"
#include "bitecoin_endpoint_client.hpp"
#include "bitecoin_hashing_miner.hpp"

#define __STDC_FORMAT_MACROS
#include <inttypes.h>
#include <cstdio>

namespace bitecoin
{

class EndpointMiner : public EndpointClient
{
 private:
  EndpointMiner(EndpointMiner &) = delete;
  void operator=(const EndpointMiner &) = delete;

  unsigned m_knownRounds;
  std::map<std::string, unsigned> m_knownCoins;

 public:
  EndpointMiner(std::string clientId, std::string minerId,
                std::unique_ptr<Connection> &conn, std::shared_ptr<ILog> &log)
      : EndpointClient(clientId, minerId, conn, log)
  {
  }

  /* Here is a default implementation of make bid.
          I would suggest that you override this method as a starting point.
  */
  virtual void MakeBid(
      const std::shared_ptr<Packet_ServerBeginRound> roundInfo,  // Information
                                                                 // about this
                                                                 // particular
                                                                 // round
      const std::shared_ptr<Packet_ServerRequestBid> request,    // The specific
                                                                 // request we
                                                                 // received
      double period,        // How long this bidding period will last
      double skewEstimate,  // An estimate of the time difference between us and
                            // the server (positive -> we are ahead)
      std::vector<uint32_t> &solution,  // Our vector of indices describing the
                                        // solution
      uint32_t *pProof  // Will contain the "proof", which is just the value
      )
  {
    double tSafetyMargin =
        0.5;  // accounts for uncertainty in network conditions
    /* This is when the server has said all bids must be produced by, plus the
            adjustment for clock skew, and the safety margin
    */

    tSafetyMargin += 0.3;  // from binning //TODO

    double tFinish =
        request->timeStampReceiveBids * 1e-9 + skewEstimate - tSafetyMargin;

    Log(Log_Verbose, "MakeBid - start, total period=%lg.", period);

    /*
            We will use this to track the best solution we have created so far.
    */
    bigint_t bestProof;
    wide_ones(BIGINT_WORDS, bestProof.limbs);

    // TODO: incremental hasher temporarily disabled to work with local server
    // static size_t previousChainSize = 0;
    // uint64_t chainHash = fnvIterative::getInstance()(
    //    (const char *)&roundInfo->chainData[previousChainSize],
    //    roundInfo->chainData.size() - previousChainSize);
    // previousChainSize = roundInfo->chainData.size();

    hash::fnv<64> hasher;
    uint64_t chainHash = hasher((const char *)&roundInfo->chainData[0],
                                roundInfo->chainData.size());

    // TODO: kill
    std::vector<uint32_t> bestSolution(3u);
    // const unsigned BIN_SIZE = 465;
    // uint32_t indx[3][BIN_SIZE];
    // bigint_t hashes[3][BIN_SIZE];

    struct point_top
    {
      uint64_t msdw;
      uint32_t indx;

      bool operator<(point_top const &other) const { return msdw < other.msdw; }
    };

    const unsigned ptvct_sz = 1 << 16;

    // point vector
    std::vector<point_top> pts;

    std::random_device rd;
    std::minstd_rand gen(rd());
    std::uniform_int_distribution<uint32_t> dis;

    // generate point vector
    for (unsigned i = 0; i < ptvct_sz; ++i) {
      uint32_t id = dis(gen);
      point_top pt;

      bigint_t temphash = PoolHashMiner(roundInfo.get(), id, chainHash);

      pt.indx = id;
      pt.msdw =
          uint64_t(temphash.limbs[6]) | (uint64_t(temphash.limbs[7]) << 32);

      pts.push_back(pt);

      // fprintf(stderr, "id: %x\n", id);
      // fprintf(stderr, "before: %8x %8x\n", temphash.limbs[7],
      // temphash.limbs[6]);
      // uint64_t sdh = uint64_t(temphash.limbs[6]) |
      // (uint64_t(temphash.limbs[7]) << 32);
      // fprintf(stderr, "after: %8x %8x\n\n", sdh>>32 , sdh & 0xffffFFFF);
    }

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
        //fprintf(stderr, "FOUND better diff\n");
      }
      // fprintf(stderr, "\n");
    }
    Log(Log_Fatal, "Best diff: %016" PRIx64 "", best_diff);
    Log(Log_Fatal, "Best offset: %8x", best_offset);

    //while (1)
    //  ;

    double t1 = now() * 1e-9;

    unsigned nTrials = 0;

    while (1) {
      ++nTrials;

      Log(Log_Debug, "Trial %d.", nTrials);

      // for (unsigned i = 0; i < BIN_SIZE; ++i) {
      //  indx[0][i] = (retained + i) & 0x3fffffff;
      //  hashes[0][i] = PoolHashMiner(roundInfo.get(), indx[0][i], chainHash);
      //}
      // for (unsigned i = 0; i < BIN_SIZE; ++i) {
      //  indx[1][i] = (1 << 30) | ((retained + i) & 0x3fffffff);
      //  hashes[1][i] = PoolHashMiner(roundInfo.get(), indx[1][i], chainHash);
      //}
      // for (unsigned i = 0; i < BIN_SIZE; ++i) {
      //  indx[2][i] = (1 << 31) | ((retained + i) & 0x3fffffff);
      //  hashes[2][i] = PoolHashMiner(roundInfo.get(), indx[2][i], chainHash);
      //}
      // retained += BIN_SIZE;

      // TODO: consider trying 1-tuples as well

      // for (unsigned i = 0; i < BIN_SIZE; ++i) {
      //  for (unsigned j = 0; j < BIN_SIZE; ++j) {
      //    bigint_t ij_acc;
      //    wide_xor(8, ij_acc.limbs, hashes[0][i].limbs, hashes[1][j].limbs);
      //    // for (unsigned p = 0; p < 8; ++p) {
      //    //   ij_acc.limbs[p] = hashes[0][i].limbs[p] ^
      // hashes[1][j].limbs[p];
      //    // }
      //    for (unsigned k = 0; k < BIN_SIZE; ++k) {
      //      bigint_t proof;
      //      //          wide_xor(8, proof.limbs, ij_acc.limbs,
      //      // hashes[2][k].limbs);
      //      for (unsigned p = 0; p < 8; ++p) {
      //        proof.limbs[p] = ij_acc.limbs[p] ^ hashes[2][k].limbs[p];
      //      }
      //      if (wide_compare(BIGINT_WORDS, proof.limbs, bestProof.limbs) < 0)
      // {
      //        bestProof = proof;
      //        bestSolution[0] = indx[0][i];
      //        bestSolution[1] = indx[1][j];
      //        bestSolution[2] = indx[2][k];
      //      }
      //    }
      //  }
      //}

      double t = now() * 1e-9;  // Work out where we are against the deadline
      double timeBudget = tFinish - t;
      Log(Log_Debug, "Finish trial %d, time remaining =%lg seconds.", nTrials,
          timeBudget);

      if (timeBudget <= 0)
        break;  // We have run out of time, send what we have
    }
    Log(Log_Info, "nTrials: %u", nTrials);
    double t2 = now() * 1e-9;
    // Log(Log_Info, "Effective hashrate: %lf",
    //    ((double)nTrials * BIN_SIZE * BIN_SIZE * BIN_SIZE) / (t2 - t1));

    solution = bestSolution;
    wide_copy(BIGINT_WORDS, pProof, bestProof.limbs);

    Log(Log_Verbose, "MakeBid - finish.");
  }
};

};  // bitecoin

#endif
