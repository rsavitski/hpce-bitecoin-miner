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

    tSafetyMargin += 0.3;  // from binning //TODO

    double tFinish =
        request->timeStampReceiveBids * 1e-9 + skewEstimate - tSafetyMargin;

    Log(Log_Verbose, "MakeBid - start, total period=%lg.", period);

    /*
            We will use this to track the best solution we have created so far.
    */
    bigint_t bestProof;
    wide_ones(BIGINT_WORDS, bestProof.limbs);

    //    static size_t previousChainSize = 0;
    //    uint64_t chainHash = fnvIterative::getInstance()(
    //        (const char *)&roundInfo->chainData[previousChainSize],
    //        roundInfo->chainData.size() - previousChainSize);
    //    previousChainSize = roundInfo->chainData.size();
    hash::fnv<64> hasher;
    uint64_t chainHash = hasher((const char *)&roundInfo->chainData[0],
                                roundInfo->chainData.size());

    // std::vector<uint32_t> bestSolution(roundInfo->maxIndices);
    std::vector<uint32_t> bestSolution(2u);
    const unsigned BIN_SIZE = 0xFFFFffff - 0x94632009;//100;//1<<32;
    //uint32_t indx[3][BIN_SIZE];
    //bigint_t hashes[3][BIN_SIZE];

    double temp_time = 0;
    double t1 = now() * 1e-9;

    unsigned nTrials = 0;

    // unsigned retained = 0;
    while (1) {
      ++nTrials;

      Log(Log_Debug, "Trial %d.", nTrials);

      //for (unsigned i = 0; i < BIN_SIZE; ++i) {
      //  indx[0][i] = i;
      //  hashes[0][i] = PoolHashMiner(roundInfo.get(), indx[0][i], chainHash);
      //  indx[1][i] = i + 0x94632009;//0x5178c98e;//0x94632009;//0xf233c9d;
      //  hashes[1][i] = PoolHashMiner(roundInfo.get(), indx[1][i], chainHash);
      //}

      //for (unsigned i = 0; i < BIN_SIZE; ++i) {
      //  indx[0][i] = i;
      //  hashes[0][i] = PoolHashMiner(roundInfo.get(), indx[0][i], chainHash);
      //  indx[1][i] = indx[0][i] + 0x94632009;//0x5178c98e;//0x94632009;//0xf233c9d;
      //  hashes[1][i] = PoolHashMiner(roundInfo.get(), indx[1][i], chainHash);
      //}

      uint32_t bitcounts[32] = {0};


      for (unsigned i = 0; i < 1<<24/*BIN_SIZE*/; ++i) {
        
        uint32_t indx1 = i;
        bigint_t hash1  = PoolHashMiner(roundInfo.get(), indx1, chainHash);
        uint32_t indx2 = i+0x94632009;;
        bigint_t hash2  = PoolHashMiner(roundInfo.get(), indx2, chainHash);
        
        
          bigint_t bestProof;
          wide_xor(8, bestProof.limbs, hash1.limbs, hash2.limbs);


          uint32_t v = bestProof.limbs[7];
          unsigned int c; // c accumulates the total bits set in v
          for (c = 0; v; c++)
          {
              v &= v - 1; // clear the least significant bit set
          }
          bitcounts[c]++;
          //fprintf(stderr, "val: %8x", bestProof.limbs[7]);
          //fprintf(stderr, "BITS SET: %u\n", c);
          //continue;

          if (i%10000 == 0)
            fprintf(stderr, "%d\n", i);

          //bigint_t diff;
          //wide_sub(8, diff.limbs, hashes[0][i].limbs, hashes[1][i].limbs);
          //  fprintf( stderr, "%6i ::: ", i);
          //for (int z=7; z>=0; --z)
          //  fprintf(stderr, "%8x ", diff.limbs[z]);
          //  fprintf(stderr, "\n");
      }
      for (int i=0; i<32; ++i){
        fprintf(stderr, "bin[%2d]: 0x%X\n", i,bitcounts[i]);
      }

      while(1);
      // for (unsigned i = 0; i < BIN_SIZE; ++i) {
      //  for (unsigned j = 0; j < BIN_SIZE; ++j) {
      //    bigint_t ij_acc;
      //    wide_xor(8, ij_acc.limbs, hashes[0][i].limbs, hashes[1][j].limbs);
      //    bigint_t proof = ij_acc;
      //    if (wide_compare(BIGINT_WORDS, proof.limbs, bestProof.limbs) < 0) {
      //      bestProof = proof;
      //      bestSolution[0] = indx[0][i];
      //      bestSolution[1] = indx[1][j];
      //      // bestSolution[2] = indx[2][k];
      //    }
      //  }
      //}

      double t = now() * 1e-9;  // Work out where we are against the deadline
      double timeBudget = tFinish - t;
      Log(Log_Debug, "Finish trial %d, time remaining =%lg seconds.", nTrials,
          timeBudget);
      Log(Log_Debug, "time taken for %u: %lf", BIN_SIZE, t - temp_time);
      temp_time = t;

      if (timeBudget <= 0)
        break;  // We have run out of time, send what we have
    }

    for (int z = 0; z < 8; ++z) {
      Log(Log_Fatal, "proof[%d]: %8x", z, bestProof.limbs[z]);
    }
    Log(Log_Fatal, "best proof: %lg",
        wide_as_double(BIGINT_WORDS, bestProof.limbs));
    Log(Log_Fatal, "ind[0]: %8x", bestSolution[0]);
    Log(Log_Fatal, "ind[1]: %8x", bestSolution[1]);
    Log(Log_Fatal, "sub   : %8x", bestSolution[1] - bestSolution[0]);

    Log(Log_Info, "nTrials: %u", nTrials);
    double t2 = now() * 1e-9;
    Log(Log_Info, "Effective hashrate: %lf",
        ((double)nTrials * BIN_SIZE * BIN_SIZE * BIN_SIZE) / (t2 - t1));

    solution = bestSolution;
    wide_copy(BIGINT_WORDS, pProof, bestProof.limbs);

    Log(Log_Verbose, "MakeBid - finish.");
  }
};

};  // bitecoin

#endif
