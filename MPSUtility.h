#ifndef __MPSUTILITY_H_CMC__
#define __MPSUTILITY_H_CMC__
#include "itensor/mps/mpo.h"
using namespace std;
using namespace itensor;

void print_MPO_tensors (const MPO& mpo, int i)
{
    cout << "site " << i << endl;
    auto il = leftLinkIndex (mpo, i);
    auto ir = rightLinkIndex (mpo, i);
    for(int i1 = 1; i1 <= il.dim(); i1++)
        for(int i2 = 1; i2 <= ir.dim(); i2++)
        {
            auto Wij = mpo(i);
            auto is = findIndex (Wij, "Site,0");
            if (il) Wij *= setElt(dag(il)=i1);
            if (ir) Wij *= setElt(dag(ir)=i2);
            if (norm(Wij) != 0.)
            {
                if (!il)
                    cout << "iright = " << i2 << endl;
                else if (!ir)
                    cout << "ileft = " << i1 << endl;
                else
                    cout << "ileft, iright = " << i1 << " " << i2 << endl;
                Wij.permute ({is, prime(is)});
                PrintData(Wij);
            }
        }
}

template <typename MPSTYPE>
inline Index leftIndex (const MPSTYPE& mps, int i)
{
    if (i == 1)
        return uniqueIndex (mps(1), mps(2), "Link");
    else
        return leftLinkIndex (mps, i);
}

template <typename MPSTYPE>
inline Index rightIndex (const MPSTYPE& mps, int i)
{
    int N = length (mps);
    if (i == N)
        return uniqueIndex (mps(N), mps(N-1), "Link");
    else
        return rightLinkIndex (mps, i);
}

Real inner_imps (const MPS& mps1, const MPS& mps2)
{
    int N = length(mps1);
    auto L = delta (dag(leftIndex(mps1,1)), prime(leftIndex(mps2,1)));
    for(int i = 1; i < N; i++)
    {
        L *= mps1(i);
        L *= prime(dag(mps2(i)),"Link");
    }
    L *= mps1(N);
    auto il = leftIndex (mps1, N);
    L *= prime(dag(mps2(N)),{il});
    return elt(L);
}

void renewLinkInds (MPS& mps, int ibeg, int iend)
{
    auto ir0 = rightLinkIndex (mps, ibeg);
    auto ir = sim(ir0);
    mps.ref(ibeg).replaceInds ({ir0}, {ir});
    for(int i = ibeg+1; i < iend; i++)
    {
        // Left index
        auto il0 = leftLinkIndex (mps, i);
        mps.ref(i).replaceInds ({il0}, {dag(ir)});
        // Right index
        ir0 = rightLinkIndex (mps, i);
        ir = sim(ir0);
        mps.ref(i).replaceInds ({ir0}, {ir});
    }
    auto il0 = leftLinkIndex (mps, iend);
    mps.ref(iend).replaceInds ({il0}, {dag(ir)});
}

#endif
