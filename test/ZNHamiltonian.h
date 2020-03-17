#ifndef __ZNHAMILTONIAN_H_CMC__
#define __ZNHAMILTONIAN_H_CMC__
#include "itensor/all_basic.h"
#include "itensor/mps/sweeps.h"
#include "itensor/mps/autompo.h"
#include "itensor/mps/mps.h"
#include "itensor/mps/dmrg.h"
#include "StringUtility.h"
#include "GaugeLadder.h"
using namespace itensor;
using namespace std;

template <typename SitesT>
MPS empty_state (const SitesT& sites)
{
    InitState init (sites);
    for(int i = 1; i <= length(sites); i++)
        init.set (i,"0");
    return MPS (init);
}

template <typename SitesT>
MPS random_state (const SitesT& sites, int ZN, int Npar, unsigned long int seed=0)
{
    if (seed == 0)
    {
        std::random_device rd;
        seed = rd();
    }
    cout << __FUNCTION__ << ": seed = " << seed << endl;

    vector<int> notfull (length(sites));
    for(size_t i = 0; i < notfull.size(); i++)
        notfull.at(i) = i+1;

    vector<int> states_tmp (length(sites)+1,0);
    std::mt19937_64 gen(seed);
    while (Npar > 0)
    {
        std::uniform_int_distribution<> dis (0, notfull.size()-1);
        int j = dis(gen);
        int site = notfull.at(j);

        auto& state = states_tmp.at(site);
        state++;
        Npar--;
        if (state == ZN-1)
            notfull.erase (notfull.begin()+j);
    }

    vector<string> states;
    for(auto const st : states_tmp)
        states.push_back (to_string(st));

    cout << "Random state:" << endl;
    InitState initState (sites);
    for(int i = 1; i <= length(sites); i++)
    {
        cout << states.at(i) << " ";
        initState.set (i,states.at(i));
    }
    cout << endl;
    return MPS (initState);
}

void add_H (AutoMPO& ampo, Real coef, vector<string> ops, vector<int> sites, bool HC=true)
{
    // Generate strings for Hermitian conjugate
    vector<string> ops_dag = ops;
    for(auto& op : ops_dag)
    {
        if (find (op,"Dag"))
            eraseSubStr (op,"Dag");
        else
            op += "Dag";
    }

    if (ops.size() == 1)
    {
        ampo += coef+0i,ops.at(0),sites.at(0);
        if (HC) ampo += coef+0i,ops_dag.at(0),sites.at(0);
    }
    else if (ops.size() == 2)
    {
        ampo += coef+0i,ops.at(0),sites.at(0),ops.at(1),sites.at(1);
        if (HC) ampo += coef+0i,ops_dag.at(1),sites.at(1),ops_dag.at(0),sites.at(0);
    }
    else if (ops.size() == 3)
    {
        ampo += coef+0i,ops.at(0),sites.at(0),ops.at(1),sites.at(1),ops.at(2),sites.at(2);
        if (HC) ampo += coef+0i,ops_dag.at(2),sites.at(2),ops_dag.at(1),sites.at(1),ops_dag.at(0),sites.at(0);
    }
    else if (ops.size() == 4)
    {
        ampo += coef+0i,ops.at(0),sites.at(0),ops.at(1),sites.at(1),ops.at(2),sites.at(2),ops.at(3),sites.at(3);
        if (HC) ampo += coef+0i,ops_dag.at(3),sites.at(3),ops_dag.at(2),sites.at(2),ops_dag.at(1),sites.at(1),ops_dag.at(0),sites.at(0);
    }
    else
    {
        cout << __FUNCTION__ << ": " << ops.size() << "-site operator not yet defined" << endl;
        throw;
    }
}

tuple <vector<string>, vector<int>> vertex_operator (const NeighborData& vertex)
{
    vector<string> ops;
    vector<int> pis;
    if (vertex.l() != -1)
    {
        ops.push_back ("Tau");
        pis.push_back (vertex.l());
    }
    if (vertex.u() != -1)
    {
        ops.push_back ("TauDag");
        pis.push_back (vertex.u());
    }
    if (vertex.d() != -1)
    {
        ops.push_back ("Tau");
        pis.push_back (vertex.d());
    }
    if (vertex.r() != -1)
    {
        ops.push_back ("TauDag");
        pis.push_back (vertex.r());
    }
    return make_tuple (ops, pis);
}

template <typename LatticeT, typename SitesT>
AutoMPO ZN_effective (const LatticeT& latt, const SitesT& sites, Real g1, Real g2, Real lambda1, Real lambda2)
{
    int N = itensor::length(sites);
    if (latt.Nf() != N)
    {
        cout << __FUNCTION__ << ": Error: size not match" << endl;
        throw;
    }

    AutoMPO ampo (sites);

    for(int i = 1; i <= N; ++i)
    {
        add_H (ampo, -g1, {"Tau"}, {i});
        add_H (ampo, -lambda1, {"Sig"}, {i});
    }
    // Plaque terms
    auto plaques = latt.plaques();
    for(auto pq : plaques)
    {
        vector<string> ops;
        vector<int> pis;
        if (pq.l() != -1)
        {
            ops.push_back ("SigDag");
            pis.push_back (pq.l());
        }
        if (pq.u() != -1)
        {
            ops.push_back ("SigDag");
            pis.push_back (pq.u());
        }
        if (pq.d() != -1)
        {
            ops.push_back ("Sig");
            pis.push_back (pq.d());
        }
        if (pq.r() != -1)
        {
            ops.push_back ("Sig");
            pis.push_back (pq.r());
        }
        add_H (ampo, -1./g2, ops, pis);
    }

    // Vertex terms
    auto vertecies = latt.vertecies();
    for(auto vt : vertecies)
    {
        auto [ops, pis] = vertex_operator (vt);
        add_H (ampo, -1./lambda2, ops, pis);
    }
  return ampo;
}
#endif
