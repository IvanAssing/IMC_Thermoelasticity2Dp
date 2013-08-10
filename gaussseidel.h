#ifndef GAUSSSEIDEL_H
#define GAUSSSEIDEL_H

#include "imc_dfm.h"

class GaussSeidel
{
    public:
        tFloat *x, *ax, *b;
        tInteger bwidth, neqmax, neq, imax;
        tInteger **AIndex, *neIndex;
        tFloat **A;
        tFloat itol;

        GaussSeidel(tFloat *results, tInteger equationMax, tInteger iterationMax = 100, tFloat iterationTolerance = 1.0e-25q, tInteger bandWidth = 5);
        void operator()(tInteger equation, tInteger index, tFloat value);
        void operator()(tInteger equation, tFloat a, tFloat b);
        void solver();
};

#endif // GAUSSSEIDEL_H
