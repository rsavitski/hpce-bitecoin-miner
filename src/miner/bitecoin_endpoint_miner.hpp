#ifndef bitecoin_miner_endpoint_hpp
#define bitecoin_miner_endpoint_hpp

#include <cstdio>
#include <cstdarg>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <cassert>

#include <vector>
#include <memory>
#include <map>

#include "bitecoin_protocol.hpp"
#include "bitecoin_endpoint.hpp"
#include "bitecoin_endpoint_client.hpp"
#include "bitecoin_hashing_miner.hpp"

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
    double tFinish =
        request->timeStampReceiveBids * 1e-9 + skewEstimate - tSafetyMargin;

    Log(Log_Verbose, "MakeBid - start, total period=%lg.", period);

    /*
            We will use this to track the best solution we have created so far.
    */
    std::vector<uint32_t> bestSolution(roundInfo->maxIndices);
    bigint_t bestProof;
    wide_ones(BIGINT_WORDS, bestProof.limbs);

    static size_t previousChainSize = 0;
    uint64_t chainHash = fnvIterative::getInstance()(
        (const char *)&roundInfo->chainData[previousChainSize],
        roundInfo->chainData.size() - previousChainSize);
    previousChainSize = roundInfo->chainData.size();

    // do only 2 indices
    uint32_t capped_num = std::min(roundInfo->maxIndices, 2u);
    std::vector<uint32_t> indices(capped_num);
    indices[0] = 1 + (rand() % 10);
    indices[1] = indices[0] + 1 + (rand() % 10);

    unsigned nTrials = 0;
    while (1) {
      ++nTrials;

      Log(Log_Debug, "Trial %d.", nTrials);

      // TODO: why +1 to both terrible? salt being zero? c being ffff?
      indices[0] += 1 + (rand() % 10);
      // indices[1] = indices[0] + 1 + (rand() % 10);
      indices[1] = 0xffffffff - indices[0];

      bigint_t acc;
      for (unsigned i = 0; i < indices.size(); i++) {
        // Calculate the hash for this specific point
        bigint_t point = PoolHashMiner(roundInfo.get(), indices[i], chainHash);

        // Combine the hashes of the points together using xor
        wide_xor(8, acc.limbs, acc.limbs, point.limbs);
      }
      bigint_t proof = acc;

      // bigint_t proof = HashMiner(roundInfo.get(), indices.size(),
      // &indices[0]);
      double score = wide_as_double(BIGINT_WORDS, proof.limbs);
      //Log(Log_Debug, "    Score=%lg", score);

      if (wide_compare(BIGINT_WORDS, proof.limbs, bestProof.limbs) < 0) {
        double worst =
            pow(2.0, BIGINT_LENGTH * 8);  // This is the worst possible score
        Log(Log_Verbose,
            "    Found new best, nTrials=%d, score=%lg, ratio=%lg.", nTrials,
            score, worst / score);
        bestSolution = indices;
        bestProof = proof;
      }

      double t = now() * 1e-9;  // Work out where we are against the deadline
      double timeBudget = tFinish - t;
      Log(Log_Debug, "Finish trial %d, time remaining =%lg seconds.", nTrials,
          timeBudget);

      if (timeBudget <= 0)
        break;  // We have run out of time, send what we have
    }

    solution = bestSolution;
    wide_copy(BIGINT_WORDS, pProof, bestProof.limbs);

    Log(Log_Verbose, "MakeBid - finish.");
  }
};

};  // bitecoin

#endif
