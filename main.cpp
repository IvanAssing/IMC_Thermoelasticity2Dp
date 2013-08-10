#include <iostream>

#include "imc_dfm.h"
#include "diffusion2dp.h"



int main()
{
    Constant2D zero(0.0q);

    Boundary2D bn(Dirichlet, &zero);

    return 0;
}

