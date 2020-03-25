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

// For fermionic tensors
ITensor parity_sign_tensor (const Index& ii)
{
    Index iip = prime(dag(ii));
    ITensor s (ii, iip);

    int block_ipre = 0;
    for(int i = 1; i <= nblock(ii); i++)
    // For each QN block
    {
        Real a = (qn(ii,i).val("Nf") % 2 == 1 ? -1. : 1.);
        int bsize = blocksize (ii, i);
        for(int j = 1; j <= bsize; j++)
        // For each element in the block
        {
            int k = j + block_ipre;
            s.set (ii=k, iip=k, a);
        }
        block_ipre += bsize;
    }
    return dag(s);
}


// ==========
// Contract the transfer matrix

inline void contract_transfer (ITensor& E, const ITensor& A)
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


void contract_transfer_matrix (ITensor& re, const SiteSet& sites, const MPS& mps, int i, const vector<string>& ops, Direction close, bool fermionic)
// The operators are applied in the reversed order in <ops>
{
    assert (close != BothDir);

    ITensor A = mps(i);
    for(int j = ops.size()-1; j >= 0; j--)
    {
        const auto& op = ops.at(j);
        A *= sites.op(op,i);
        A.noPrime();
    }
    if (fermionic)
    {
        if (ops.size() % 2 != 0 && i != 1)
        {
            // Fermionic sign
            auto iL = leftLinkIndex (mps, i);
            auto F = parity_sign_tensor (iL);
            A *= F;
            A.noPrime();
        }
    }

    ITensor Ap = dag(mps(i));
    if (close == NoDir)
        Ap.prime ("Link");
    else if (close == Fromleft)
        Ap.prime (rightLinkIndex (mps, i));
    else if (close == Fromright)
        Ap.prime (leftLinkIndex (mps, i));
    re *= A;
    re *= Ap;
}

#endif
