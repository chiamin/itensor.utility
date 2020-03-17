#ifndef __ITENSORUTILITY_H_CMC__
#define __ITENSORUTILITY_H_CMC__
#include <vector>
#include "itensor/all.h"
using namespace std;
using namespace itensor;

// For ITensor3

vector<ITensor> get_REs (const MPS& mps)
{
    int N = length (mps);
    vector<ITensor> Rs (N+1);
    Rs.at(N) = ITensor(1);
    for(int i = N-1; i >= 1; i--)
    {
        Rs.at(i) = Rs.at(i+1) * mps(i+1) * dag(prime(mps(i+1), "Link"));
    }
    return Rs;
}


// ==========
// Contract the transfer matrix
inline void contract_transfer (itensor::ITensor& E, const itensor::ITensor& A)
{
    E *= A;
    E *= dag(prime(A,"Link"));
}

inline void contract_transfer (ITensor& E, const ITensor& A, const ITensor& op)
{
    E *= A;
    E *= op;
    E.noPrime ("Site");
    E *= dag(prime(A,"Link"));
}

inline void contract_transfer (ITensor& E, const MPS& mps, int i)
{
    E *= mps(i);
    E *= dag(prime(mps(i),"Link"));
}

// Contract the transfer matrix
template <typename SitesT>
inline void contract_transfer (ITensor& E, const MPS& mps, int i, const SitesT& sites, string op)
{
    E *= mps(i);
    E *= sites.op(op,i);
    E.noPrime ("Site");
    E *= dag(prime(mps(i),"Link"));
}
#endif
