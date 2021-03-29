#ifndef __GATEGONTAINER_H_CMC__
#define __GATEGONTAINER_H_CMC__
#include <initializer_list>
#include "itensor/all.h"
#include "IUtility.h"

struct GateBase
{
    public:
        GateBase (const ITensor& gate, const Index& ii)
        : _gate (gate)
        , _inds ({ii, prime(ii)})
        , _size (1)
        {}
        GateBase (const ITensor& gate, const Index& i1, const Index& i2)
        : _gate (gate)
        , _inds ({i1, i2, prime(i1), prime(i2)})
        , _size (2)
        {}
        template <typename GenerateFunction, typename... FuncArgs>
        GateBase (GenerateFunction GateGen, const Index& i1, FuncArgs... args)
        : _gate (GateGen(i1,args...))
        , _inds ({i1, prime(i1)})
        , _size (1)
        {}
        template <typename GenerateFunction, typename... FuncArgs>
        GateBase (GenerateFunction GateGen, const Index& i1, const Index& i2, FuncArgs... args)
        : _gate (GateGen(i1,i2,args...))
        , _inds ({i1, i2, prime(i1), prime(i2)})
        , _size (2)
        {}

        template <typename IndsType>
        ITensor gen (const IndsType& inds) const
        {
            if (_inds == inds)
                return _gate;
            else
                return replaceInds (_gate, _inds, inds);
        }

        int size () const { return _size; }

    private:
        ITensor _gate;
        vector<Index> _inds;
        int _size;
};

struct ApplyGateInfo
{
    string name;
    int pos;
    bool move_oc;
};

class GateContainer
{
    public:
        void new_gate (const string& name, const ITensor& gate, const Index& ii)
        {
            _gates.emplace (name, GateBase(gate, ii));
        }
        void new_gate (const string& name, const ITensor& gate, const Index& i1, const Index& i2)
        {
            _gates.emplace (name, GateBase(gate, i1, i2));
        }
        template <typename GenerateFunction, typename... FuncArgs>
        void new_gate (const string& name, GenerateFunction GateGen, FuncArgs... args)
        {
            _gates.emplace (name, GateBase(GateGen, args...));
        }

        void add (const string& name, int pos, bool move_oc=true)
        {
            ApplyGateInfo info {name, pos, move_oc};
            _info.push_back (info);
        }
        void add_front (const string& name, int pos, bool move_oc=true)
        {
            ApplyGateInfo inf   o {name, pos, move_oc};
            _info.insert (_info.begin(), info);
        }

        void clear () { _gates.clear(); }
        template <typename MPSType>
        tuple<int,int> apply (MPSType& psi, const Args& args=Args::global(),
                              int leftEdge=std::numeric_limits<int>::min(),
                              int rightEdge=std::numeric_limits<int>::max());

    private:
        unordered_map <string,GateBase> _gates;
        vector<ApplyGateInfo>           _info;
};

template <typename MPSType>
tuple<int,int>
GateContainer :: apply
(MPSType& psi, const Args& args, int leftEdge, int rightEdge)
{
    auto iss = IndexSet();
    if constexpr (is_same_v<MPSType,MPS>)
        iss = siteInds (psi);
    else // MPO
    {
        vector<Index> ii;
        ii.reserve (length(psi));
        for(int i = 1; i <= length(psi); i++)
        {
            auto is = findIndex (psi(i), "Site,0");
            ii.push_back (dag(is));
        }
        iss = IndexSet (ii);
    }
    for(int j = 0; j < _info.size(); j++)
    {
        auto info = _info.at(j);
        auto const& gatebase = _gates.at (info.name);
        // Site to apply
        int i = info.pos;
        if (gatebase.size() == 1)
        {
            if (i < leftEdge || i > rightEdge) continue;
            // Generate gate
            vector<Index> inds {iss(i), prime(iss(i))};
            ITensor gate = gatebase.gen (inds);
            // Apply gate
            int oc = orthoCenter (psi);
            if constexpr (is_same_v <MPSType, MPS>)
            {
                psi.ref(i) *= gate;
                psi.ref(i).noPrime("Site");
            }
            else // MPO
            {
                psi.ref(i) *= prime(gate);
                psi.ref(i) *= dag(gate);
                psi.ref(i).prime(-1,"Site");
            }
            psi.leftLim (oc-1);
            psi.rightLim (oc+1);
        }
        else if (gatebase.size() == 2)
        {
            //if (i+1 == leftEdge) leftEdge++;
            //if (i == rightEdge) rightEdge--;
            if (i < leftEdge || i > rightEdge) continue;
            // Generate gate
            vector<Index> inds {iss(i), iss(i+1), prime(iss(i)), prime(iss(i+1))};
            ITensor gate = gatebase.gen (inds);
            // Move orthogonality center
            int oc = orthoCenter (psi);
            if (info.move_oc && i != oc && i+1 != oc)
            {
                psi.position (abs(oc-i) < abs(oc-i-1) ? i : i+1);
                oc = orthoCenter (psi);
            }
            // Set direction
            Direction dir = Fromleft;
            if (j != _info.size()-1)
            {
                int inext = _info.at(j+1).pos;
                dir = (inext > i ? Fromleft : Fromright);
            }
            else if (j != 0)
            {
                int ipre = _info.at(j-1).pos;
                dir = (i > ipre ? Fromleft : Fromright);
            }
            // Apply gate
            if constexpr (is_same_v <MPSType, MPS>)
            {
                apply_gate (psi, i, gate, dir, args);
            }
            else // MPO
            {
                assert (oc == i || oc == i+1);
                auto AA = psi(i) * psi(i+1);
                AA *= prime(gate);
                AA *= dag(gate);
                AA.prime(-1,"Site");
                auto iiL = uniqueInds (psi(i), psi(i+1));
                tie (psi.ref(i), psi.ref(i+1)) = denmatDecomp (AA, iiL, dir, args);
                int new_oc = (dir == Fromleft ? i+1 : i);
                psi.leftLim (new_oc-1);
                psi.rightLim (new_oc+1);
            }
        }
        else
        {
            cout << __FUNCTION__ << ": Not support gate size: " << gatebase.size() << endl;
            throw;
        }
    }
    return {leftEdge, rightEdge};
}
#endif