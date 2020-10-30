#ifndef __ITENSORUTILITY_H_CMC__
#define __ITENSORUTILITY_H_CMC__
#include <vector>
#include "itensor/all.h"
using namespace std;
using namespace itensor;

// For ITensor3

// Contract <ten1> and <ten2> by the tags in <tags1> and <tags2>
void auto_contractEqual_by_tag (ITensor& ten1, const ITensor& ten2, const vector<string>& tags1, const vector<string>& tags2)
{
    for(int i = 0; i < tags1.size(); i++)
    {
        auto tag1 = tags1.at(i);
        auto tag2 = tags2.at(i);
        auto ii1  = findIndex (ten1, tag1);
        auto ii2  = findIndex (ten2, tag2);
        if (ii1 != ii2)
            ten1 *= dag(delta(ii1,ii2));
    }
    ten1 *= ten2;
}

inline ITensor auto_contract_by_tag (ITensor ten1, const ITensor& ten2, const vector<string>& tags1, const vector<string>& tags2)
{
    auto_contractEqual_by_tag (ten1, ten2, tags1, tags2);
    return ten1;
}

// Add <ten1> and <ten2> by the tags in <tags1> and <tags2>
void auto_addEqual_by_tag (ITensor& ten1, const ITensor& ten2, const vector<string>& tags1, const vector<string>& tags2)
{
    for(int i = 0; i < tags1.size(); i++)
    {
        auto tag1 = tags1.at(i);
        auto tag2 = tags2.at(i);
        auto ii1  = findIndex (ten1, tag1);
        auto ii2  = findIndex (ten2, tag2);
        ten1.replaceInds ({ii1},{ii2});
    }
    ten1 += ten2;
}

// Add <ten1> and <ten2> by the tags in <tags1> and <tags2>
inline ITensor auto_add_by_tag (ITensor ten1, const ITensor& ten2, const vector<string>& tags1, const vector<string>& tags2)
{
    auto_addEqual_by_tag (ten1, ten2, tags1, tags2);
    return ten1;
}

inline ITensor Identity (const Index& ii)
{
    ITensor id (dag(ii), prime(ii));
    for(int i = 1; i <= ii.dim(); i++)
        id.set (i,i,1.);
    return id;
}

inline ITensor Identity (const Index& ii, const Index& iip)
{
    ITensor id (ii, iip);
    for(int i = 1; i <= ii.dim(); i++)
        id.set (i,i,1.);
    return id;
}

inline int dim (const ITensor& T)
{
    int d = 1;
    for(const auto& ii : T.inds())
        d *= ii.dim();
    return d;
}

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

vector<ITensor> get_LEs (const MPS& mps)
{
    int N = length (mps);
    vector<ITensor> Ls (N+1);
    Ls.at(1) = ITensor(1);
    for(int i = 2; i <= N; i++)
    {
        Ls.at(i) = Ls.at(i-1) * mps(i-1) * dag(prime(mps(i-1), "Link"));
    }
    return Ls;
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
        Real a = 1.;
        if (qn(ii,i).hasName("Nf"))
            a = (qn(ii,i).val("Nf") % 2 == 1 ? -1. : 1.);
        else if (qn(ii,i).hasName("Pf"))
            a = (qn(ii,i).val("Pf") % 2 == 1 ? -1. : 1.);

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

inline ITensor merge_onsite_operators (const SiteSet& sites, int i, const vector<string>& ops)
{
    ITensor op_all (1.);
    for(int j = ops.size()-1; j >= 0; j--)
    {
        const auto& op = ops.at(j);
        if (op == "") continue;
        op_all *= prime (sites.op(op,i));
        op_all.mapPrime(1,0);
        op_all.mapPrime(2,1);
    }
    return op_all;
}

inline void apply_onsite_ops (ITensor& A, const SiteSet& sites, int i, const vector<string>& ops)
{
    assert (hasIndex (A, sites(i)));

    ITensor op_all = merge_onsite_operators (sites, i, ops);
    A *= op_all;
    A.mapPrime(1,0);
}

inline void apply_fermionic_sign (ITensor& A, const Index& iL)
{
    // Fermionic sign
    auto F = parity_sign_tensor (iL);
    A *= F;
    A.noPrime();
}

void contract_transfer_matrix (ITensor& re, const SiteSet& sites, const MPS& mps, int i, const vector<string>& ops, Direction close, bool fermionic)
// The operators are applied in the reversed order in <ops>
{
    assert (close != BothDir);

    ITensor A = mps(i);
    apply_onsite_ops (A, sites, i, ops);
    if (fermionic)
    {
        if (ops.size() % 2 != 0 && i != 1)
        {
            // Fermionic sign
            auto iL = leftLinkIndex (mps, i);
            apply_fermionic_sign (A, iL);
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
