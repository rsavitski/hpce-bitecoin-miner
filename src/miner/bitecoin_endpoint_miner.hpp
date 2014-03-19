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

    // TODO: incremental hasher temporarily disabled to work with local server
    // static size_t previousChainSize = 0;
    // uint64_t chainHash = fnvIterative::getInstance()(
    //    (const char *)&roundInfo->chainData[previousChainSize],
    //    roundInfo->chainData.size() - previousChainSize);
    // previousChainSize = roundInfo->chainData.size();

    Log(Log_Verbose, "Maxindices: %u", roundInfo->maxIndices);

    hash::fnv<64> hasher;
    uint64_t chainHash = hasher((const char *)&roundInfo->chainData[0],
                                roundInfo->chainData.size());

    struct point_top
    {
      uint64_t msdw;
      uint32_t indx;

      bool operator<(point_top const &other) const { return msdw < other.msdw; }
    };

    //Log(Log_Verbose, "Stamp after fnv");

    const unsigned ptvct_sz = 1 << 16;

    Log(Log_Verbose, "Stamp"); //TODO
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
    Log(Log_Verbose, "Stamp"); //TODO

    std::sort(pts.begin(), pts.end());

    Log(Log_Verbose, "Stamp"); //TODO
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

      // arithmetic difference for differential attack
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
    Log(Log_Fatal, "Best diff: %016" PRIx64 "", best_diff);
    Log(Log_Fatal, "Best offset: %8x", best_offset);

    pts.clear();
    Log(Log_Verbose, "Stamp before metapt gen"); //TODO

    // done with offset search
    //////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////

    struct metapoint_top
    {
      bigint_t point;
      uint32_t indx;  // base index

      bool operator<(metapoint_top const &other) const
      {
          return wide_compare(BIGINT_WORDS, point.limbs, other.point.limbs) < 0;
      }
    };

    // metapoint vector
    const unsigned metaptvct_sz =
        1 << 18;  // TODO: autotune, maybe golden diff finder too
    std::vector<metapoint_top> metapts;

    std::uniform_int_distribution<uint32_t> dis2(0, 0xFFFFFFFE - best_offset);

    // generate metapoint vector
    for (unsigned i = 0; i < metaptvct_sz; ++i) {
      uint32_t id = dis2(gen);
      metapoint_top pt;

      bigint_t temphash = PoolHashMiner(roundInfo.get(), id, chainHash);
      bigint_t temphash2 =
          PoolHashMiner(roundInfo.get(), id + best_offset, chainHash);

      wide_xor(8, pt.point.limbs, temphash.limbs, temphash2.limbs);
      pt.indx = id;

      metapts.push_back(pt);
    }

    Log(Log_Verbose, "Stamp after metapt gen"); //TODO
    // fprintf(stderr, "--------\n\n");
    // for (auto pt : metapts) {
    //  fprintf(stderr, "idx : %8x\n", pt.indx);
    //  fprintf(stderr, "msdw: %016" PRIx64 "\n", pt.msdw);
    //  fprintf(stderr, "tdw: %8x\n", pt.tdw);
    //}
    // fprintf(stderr, "--------\n\n");

    std::sort(metapts.begin(), metapts.end());
    Log(Log_Verbose, "Stamp: metapts sorted"); //TODO

    // fprintf(stderr, "--------\n\n");
    // for (auto pt : metapts) {
    //  fprintf(stderr, "idx : %8x\n", pt.indx);
    //  fprintf(stderr, "msdw: %016" PRIx64 "\n", pt.msdw);
    //  fprintf(stderr, "tdw: %8x\n", pt.tdw);
    //}
    // fprintf(stderr, "--------\n\n");

    bigint_t mbest;
    wide_ones(BIGINT_WORDS, mbest.limbs);

    uint32_t metaidx[2];

    for (auto it = metapts.begin(); it != metapts.end() - 1; ++it) {
      uint32_t curr_indx = it->indx;
      uint32_t next_indx = (it + 1)->indx;

      if (curr_indx == next_indx || curr_indx + best_offset == next_indx ||
          next_indx + best_offset == curr_indx) {
        Log(Log_Verbose,
            "Skipped identical index sample in second pass search");
        continue;
      }

      bigint_t xoredw;
      wide_xor(8, xoredw.limbs, it->point.limbs, (it+1)->point.limbs);

      if (wide_compare(BIGINT_WORDS, xoredw.limbs, mbest.limbs) < 0){
        mbest = xoredw;
        metaidx[0] = curr_indx;
        metaidx[1] = next_indx;
      }
    }
    Log(Log_Verbose, "Stamp"); //TODO

    // fprintf(stderr, "--------\n\n");
    // fprintf(stderr, "idx : %8x\n", metaidx[0]);
    // fprintf(stderr, "idx : %8x\n", metaidx[1]);
    // fprintf(stderr, "msdw: %016" PRIx64 "\n", mbest_diff.msdw);
    // fprintf(stderr, "tdw: %8x\n", mbest_diff.tdw);
    // fprintf(stderr, "--------\n\n");
    Log(Log_Fatal, "Finished metasearch");

    Log(Log_Info, "Maxindices: %u", roundInfo->maxIndices);
    uint32_t maxidx = roundInfo->maxIndices;
    uint32_t lg2idx = 0;
    while (maxidx != 0){
      maxidx >>=1;
      lg2idx++;
    }
    lg2idx -=1;
    Log(Log_Info, "Log: %u", lg2idx);
      
    ///////////////////////////////////////////

    double t1 = now() * 1e-9;
    unsigned nTrials = 0;
    uint32_t r = 0;
    while (1) {
      ++nTrials;

      Log(Log_Debug, "Trial %d.", nTrials);

      // TODO

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

};  // bitecoin

#endif
