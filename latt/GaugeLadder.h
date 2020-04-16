#ifndef __GAUGELADDER_H__
#define __GAUGELADDER_H__
#include <vector>
#include <string>
#include "SquareLattice.h"
using namespace std;

class ChargeSite
{
    public:
        ChargeSite () {}
        ChargeSite (int i_, int x_, int y_, int lc, int rc, int uc, int dc, int lf, int rf, int uf, int df)
        : _i(i_), _x(x_), _y(y_)
        , _nb_flux ({lf,df,rf,uf})
        , _nb_i_charges ({lc,dc,rc,uc})
        {}

        int i () const { return _i; }
        int x () const { return _x; }
        int y () const { return _y; }
        // flux
        int left_flux () const { return _nb_flux.at(0); }
        int down_flux () const { return _nb_flux.at(1); }
        int right_flux () const { return _nb_flux.at(2); }
        int up_flux () const { return _nb_flux.at(3); }
        const vector<int>& nb_flux () const { return _nb_flux; }
        // charge
        int left_charge () const { return _nb_i_charges.at(0); }
        int down_charge () const { return _nb_i_charges.at(1); }
        int right_charge () const { return _nb_i_charges.at(2); }
        int up_charge () const { return _nb_i_charges.at(3); }
        const vector<int>& nb_i_charge () const { return _nb_i_charges; }

    private:
        int _x, _y, _i;
        vector<int> _nb_flux, _nb_i_charges;
};

class NeighborData
{
    public:
        NeighborData (const vector<int>& lrud)
        : p (lrud)
        {}
        int l () const { return p.at(0); }
        int r () const { return p.at(1); }
        int u () const { return p.at(2); }
        int d () const { return p.at(3); }
    private:
        vector<int> p;
};

class GaugeLadder
{
    public:
        GaugeLadder (int L, bool lrough_, bool rrough_);

        int  L      () const { return _L; }
        int  Nf     () const { return _N_flux; }
        int  Nc     () const { return _N_charge; }
        bool lrough () const { return _lrough; }
        bool rrough () const { return _rrough; }

        const vector<ChargeSite>&  charges   () const { return _charges; }
              vector<NeighborData> plaques   () const;
              vector<NeighborData> vertecies () const;

        const ChargeSite& charge_dn (int x) const { return _charges.at(x*2-1); }
        const ChargeSite& charge_up (int x) const { return _charges.at(x*2); }

    private:
        int _L, _N_flux, _N_charge;
        vector<ChargeSite>  _charges;
        vector<vector<int>> _plaques;
        bool _lrough, _rrough;

        using SL = SquareLattice;
};

GaugeLadder :: GaugeLadder (int L, bool lrough_, bool rrough_)
: _L (L)
, _N_charge (L*2)
, _lrough (lrough_)
, _rrough (rrough_)
{
    bool xpbc = false, ypbc = false;
    _charges.resize (_N_charge+1);

    int i_flux = 1;
    for(int x = 1; x <= L; x++)
    {
        for(int y = 1; y <= 2; y++)
        {
            int i_charge = SL::mix_ind (x, y, L, 2);

            vector<int> flux, charges;
            for(auto dir : vector<string>{"l","d","u","r"}) // the order defines the labelling
            {
                // Check if the neighbor NOT in the lattice. If so, <nb_i_charge> = -1.
                // Store the direction, that the neighbor can find the current site, in <opp_dir>
                int x2 = x, y2 = y;
                int opp_dir;
                bool nb_in_latt = true;
                if (dir == "l") // left
                {
                    x2 = SL::apply_bc (x-1, L, xpbc);
                    if (x2 < 0) nb_in_latt = false;
                    opp_dir = 2;                            // right; defined in <ChargeSite>
                }
                else if (dir == "u") // up
                {
                    y2 = SL::apply_bc (y+1, 2, ypbc);
                    if (y2 < 0) nb_in_latt = false;
                    opp_dir = 1;                            // down; defined in <ChargeSite>
                }
                else if (dir == "d") // down
                {
                    y2 = SL::apply_bc (y-1, 2, ypbc);
                    if (y2 < 0) nb_in_latt = false;
                    opp_dir = 3;                            // up; defined in <ChargeSite>
                }
                else if (dir == "r") // right
                {
                    x2 = SL::apply_bc (x+1, L, xpbc);
                    if (x2 < 0) nb_in_latt = false;
                    opp_dir = 0;                            // left; defined in <ChargeSite>
                }
                // If neighbor NOT in the lattice,
                // 1. keep <nb_i_charge> = -1
                // 2. set <nb_flux> = -1, except for the rough boundaries
                int nb_i_charge, nb_flux;
                if (!nb_in_latt)
                {
                    nb_i_charge = -1;
                    if ((x == 1 && _lrough && dir == "l")
                    ||  (x == L && _rrough && dir == "r"))
                    {
                        nb_flux = i_flux++;
                    }
                    else
                    {
                        nb_flux = -1;
                    }
                }
                // If neighbor in the lattice,
                // 1. store the neighboring charge in <nb_i_charge>
                // 2. store the neighboring flux in <nb_flux>
                else
                {
                    // Neighboring charge
                    nb_i_charge = SL::mix_ind (x2, y2, L, 2);
                    // Find the neighboring flux from the neighboring charge
                    const vector<int>& nb_fluxs = _charges.at(nb_i_charge).nb_flux();
                    int flux_from_nb = (nb_fluxs.size() == 0 ? -1 : nb_fluxs.at(opp_dir));
                    if (flux_from_nb != -1)
                    // flux has been labelled from the neighboring charge, use the old label
                    {
                        nb_flux = flux_from_nb;
                    }
                    else // add a new label to the flux
                    {
                        nb_flux = i_flux++;
                    }
                }
                charges.push_back (nb_i_charge);
                flux.push_back (nb_flux);
            }
            _charges.at(i_charge) = ChargeSite (i_charge, x, y,
                                                charges.at(0), charges.at(3), charges.at(2), charges.at(1),
                                                flux.at(0), flux.at(3), flux.at(2), flux.at(1));
        }
    }
    _N_flux = i_flux-1;
}

// Return the gauge fluxes in the plaques
// The labeling of charge and gauge flux
//
//        (4)--2--(5)--4--(8)--6--(11)--8--(14)
//             |       |       |        |
// rough      (2)     (6)     (9)     (12)       rough
//             |       |       |        |
//        (1)--1--(3)--3--(7)--5--(10)--7--(13)
//
//            2--(3)--4--(6)--6--(9)--8
//            |       |       |       |
// smooth    (1)     (4)     (7)    (10)  smooth
//            |       |       |       |
//            1--(2)--3--(5)--5--(8)--7
vector<NeighborData> GaugeLadder :: plaques () const
{
    bool xpbc = false;
    vector<NeighborData> pqs;
    int y = 1;
    vector<int> lrud (4);
    int L = this->_L;

    // left most
    int iL = SL::mix_ind (1, y, L, 2);
    const auto& charge_L = this->charges().at(iL);
    if (charge_L.left_flux() != -1) // left rough boundary
    {
        lrud.at(0) = -1;
        lrud.at(1) = charge_L.up_flux();
        int up_i = SL::mix_ind (1, y+1, L, 2);
        lrud.at(2) = this->charges().at(up_i).left_flux();
        lrud.at(3) = charge_L.left_flux();
        pqs.push_back (NeighborData(lrud));
    }
    for(int x = 2; x <= L; x++)
    {
        int i = SL::mix_ind (x, y, L, 2);
        // left
        int xn = SL::apply_bc (x-1, L, xpbc);
        if (xn > 0)
        {
            int left_i = SL::mix_ind (xn, y, L, 2);
            lrud.at(0) = this->charges().at(left_i).up_flux();
        }
        else
        {
            lrud.at(0) = -1;
        }
        // right
        lrud.at(1) = this->charges().at(i).up_flux();
        // up
        int up_i = SL::mix_ind (x, y+1, L, 2);
        lrud.at(2) = this->charges().at(up_i).left_flux();
        // down
        lrud.at(3) = this->charges().at(i).left_flux();

        pqs.push_back (NeighborData(lrud));
    }
    // right most
    int iR = SL::mix_ind (L, y, L, 2);
    const auto& charge_R = this->charges().at(iR);
    if (charge_R.right_flux() != -1) // right rough boundary
    {
        lrud.at(0) = charge_R.up_flux();
        int up_i = SL::mix_ind (L, y+1, L, 2);
        lrud.at(1) = -1;
        lrud.at(2) = this->charges().at(up_i).right_flux();
        lrud.at(3) = charge_R.right_flux();
        pqs.push_back (NeighborData(lrud));
    }
    return pqs;
}

// Return the gauge fluxes in the vertices in the left, right, up, down order
// If a direction does not have gauge flux, the value is -1
// Otherwise the value is the label index
vector<NeighborData> GaugeLadder :: vertecies () const
{
    vector<NeighborData> vts;
    vector<int> lrud (4);
    for(int i = 1; i <= this->Nc(); i++)
    {
        const auto& charge = this->charges().at(i);
        lrud.at(0) = charge.left_flux();
        lrud.at(1) = charge.right_flux();
        lrud.at(2) = charge.up_flux();
        lrud.at(3) = charge.down_flux();
        vts.push_back (NeighborData(lrud));
    }
    return vts;
}
#endif
