#include <iostream>

#include "imc_dfm.h"
#include "diffusion2dp.h"
#include "gaussseidel.h"



int main()
{

    tFloat *x = new tFloat[2];

    GaussSeidel sys(x, 2);

    sys(1, 2.0, 4.0);

    return 0;
}

