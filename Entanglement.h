#ifndef __ENTANGLEMENT_H_CMC__
#define __ENTANGLEMENT_H_CMC__
#include "ContainerUtility.h"
#include "GeneralUtility.h"

Real EntangEntropy (const Spectrum& spec)
{
    Real S = 0.;
    for(int i = 1; i <= spec.size(); i++)
    {
        auto p = spec.eig(i);
        S += - p * log(p);
    }
    return S;
}

Real EntangEntropy (const ITensor& rho)
{
    mycheck (order(rho) == 2, "density matrix must has two indices");
    mycheck (rho.inds()(1) == prime(dag(rho.inds()(2))), "indices not match");
    auto i1 = rho.inds()(1);
    ITensor U(i1), S, V;
    auto spec = svd (rho, U, S, V, {"Cutoff",1e-12});
    return EntangEntropy (spec);
}

ITensor get_rho_brute_force (const MPS& psi, vector<int>& sites)
{
    std::sort (sites.begin(), sites.end());

    // Left
    ITensor L(1.);
    int imid = sites.at(sites.size()/2);
    for(int i = 1; i <= length(psi); i++)
    {
        ITensor A = psi(i);
        ITensor Adag;
        if (iut::in_vector (sites, i))
            Adag = prime(dag(A));
        else
            Adag = prime(dag(A),"Link");

        L *= A;
        L *= Adag;
        if (i == imid)
            break;
    }
    // Right
    ITensor R(1.);
    for(int i = length(psi); i >= 1; i--)
    {
        if (i == imid)
            break;
        ITensor A = psi(i);
        ITensor Adag;
        if (iut::in_vector (sites, i))
            Adag = prime(dag(A));
        else
            Adag = prime(dag(A),"Link");

        R *= A;
        R *= Adag;
    }
    ITensor rho = L * R;
    mycheck (order(rho) == 2*sites.size(), "size not match");
    return rho;
}

Real get_entang_entropy_brute_force (const MPS& psi, vector<int>& sites)
{
    auto rho = get_rho_brute_force (psi, sites);

    // Combine indices
    vector<Index> phys_inds;
    for(int k : sites)
    {
        phys_inds.push_back (siteIndex (psi,k));
    }
    auto [comb, cindex] = combiner (phys_inds);
    rho *= comb;
    rho *= prime(dag(comb));

    auto EE = EntangEntropy (rho);
    return EE;
}
#endif
