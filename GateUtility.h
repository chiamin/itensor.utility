#ifndef __GATEUTILITY_H_CMC__
#define __GATEUTILITY_H_CMC__
#include "IUtility.h"

template <typename SiteType>
class GateContainer
{
    public:
        GateContainer (int N=2, const Args& args)
        : _sites (N, args)
        {}

        GateContainer (const SiteType& sites, int N=2, const Args& args)
        : _sites (N, args)
        , _out_sites (sites)
        {}

        void set_sites (const SiteType& sites) { _out_sites = sites; }

        void add (string name, string op)
        {
            _gates[name] = _sites.op (op, 1);
        }
        void add (string name, string op1, string op2)
        {
            _gates[name] = _sites.op (op1, 1) * _sites.op (op2, 2);
        }

        void set_exp_t (string name, Real t)
        {
            _exp_gates[name] = expHermitian (_gates[name], -1_i*t);
        }

        int n () const
        {
            return (order(_gate) / 2);
        }
        ITensor get (string name, int i) const
        {
            auto T = _gates.at(name)
            return replace_inds (T, i);
        }
        ITensor get_exp (string name, int i) const
        {
            auto T = _exp_gates.at(name);
            return replace_inds (T, i);
        }
        ITensor replace_inds (const ITensor& T, int i) const
        {
            if (order(T) == 1)
            {
                auto i1 = _sites(1);
                auto i2 = _out_sites(i);
                return replaceInds (T, {i1, prime(i1)}, {i2, prime(i2)});
            }
            else if (order(T) == 2)
            {
                auto i1 = _sites(1);
                auto i2 = _sites(2);
                auto j1 = _out_sites(i);
                auto j2 = _out_sites(i+1);
                return replaceInds (T, {i1, prime(i1), i2, prime(i2)}, {j1, prime(j1), j2, prime(j2)});
            }
        }

    private:
        unsorted_map <string,ITensor> _gates, _exp_gates;
        SiteType                      _sites, _out_sites;
};
#endif
