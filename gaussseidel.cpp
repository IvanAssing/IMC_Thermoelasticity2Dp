#include "gaussseidel.h"

GaussSeidel::GaussSeidel(tFloat *results, tInteger equationMax, tInteger iterationMax, tFloat iterationTolerance, tInteger bandWidth)
    :x(results), neqmax(equationMax), imax(iterationMax), itol(iterationTolerance), bwidth(bandWidth)
{
    A = new tFloat*[neqmax];
    for(tInteger i = 0; i<neqmax; i++)
        A[i] = new tFloat[bwidth-1];

    AIndex = new tInteger*[neqmax];
    for(tInteger i = 0; i<neqmax; i++)
        AIndex[i] = new tInteger[bwidth-1];


    ax = new tFloat[neqmax];
    b = new tFloat[neqmax];

    for(tInteger i = 0; i<neqmax; i++)
        ax[i] = 1.0q, b[i] = 0.0q;

    neIndex = new tInteger[neqmax];
    for(tInteger i = 0; i<neqmax; i++)
        neIndex[i] = 0;

    neq = 0;
}

void GaussSeidel::operator()(tInteger i, tInteger j, tFloat value)
{
    if(i<0 || i>neqmax || j<0 || j>neqmax)
        throw std::string("error: invalid index");
    else{
        A[i][neIndex[i]] = value;
        AIndex[i][neIndex[i]] = j;
        neIndex[i]++;
    }
}

void GaussSeidel::operator()(tInteger i, tFloat _a, tFloat _b)
{
    if(i<0 || i>neqmax)
        throw std::string("error: invalid index");
    else{
        ax[i] = _a;
        b[i] = _b;
    }
}

void GaussSeidel::solver()
{
    tInteger it = 0;
    tFloat sum, residual;

    do{
        it++;

        // Solver
        for(tInteger i = 0; i<neqmax; i++){
            sum = b[i];
            for(tInteger j=0; j<neIndex[i]; j++)
                sum += A[i][j]*x[AIndex[i][j]];
            x[i] = sum/ax[i];
        }

        // Residuo
        residual = 0.0q;
        for(tInteger i = 0; i<neqmax; i++){
            sum = -b[i]+ax[i]*x[i];
            for(tInteger j=0; j<neIndex[i]; j++)
                sum += A[i][j]*x[AIndex[i][j]];
            residual = sum*sum;
        }

        residual = sqrtq(residual);

    }while(residual > itol || it < imax);
}
