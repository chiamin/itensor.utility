#include "itensor/all.h"
#include "ITensorUtility.h"
#include "ZNHamiltonian.h"
using namespace std;
using namespace itensor;

int main()
{
    Z3 sites (4);
    MPS psi = random_state (sites, 3, 4);
    ITensor E (1);
    for(int i = 1; i <= 4; i++)
        contract_transfer (E, psi, i);
    cout << elt(E) << "  " << inner(psi,psi) << endl;
}
