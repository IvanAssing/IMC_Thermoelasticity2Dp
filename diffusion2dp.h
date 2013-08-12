#ifndef DIFFUSION2DP_H
#define DIFFUSION2DP_H

#include "imc_dfm.h"
#include "functor2d.h"
#include "node2d.h"
#include "boundary2d.h"
#include "data.h"

class Diffusion2Dp
{
    public:
        Diffusion2DData *data; // Dados do problema
        Node2D *nodes; // Lista de nós
        tInteger nx, ny; // Número de divisões em x e y
        tFloat lx, ly, hx, hy; // Tamanho de malha
        Boundary2D *ccS, *ccN, *ccE, *ccW; // Condições de contorno

        tFloat *T, Tm; // Vetor solução

        Diffusion2Dp(tFloat lengthX, tFloat lengthY, tInteger nx, tInteger ny, Diffusion2DData *data,
                     Boundary2D *ccSouth, Boundary2D *ccNorth, Boundary2D *ccEast, Boundary2D *ccWest);

        void solver(tInteger iterationMax = 1000, tFloat iterationTolerance = 1.0e-28q, bool plotlog = false); // Discretização + solver do sistema linear

        tInteger direction(tInteger position, DirectionType dir);
        tInteger position(tInteger i, tInteger j);

        tFloat operator()(tInteger i, tInteger j);

        void plotX(tInteger p, Functor2D &analyticalSolution);
        void plotY(tInteger p, Functor2D &analyticalSolution);
};

#endif // DIFFUSION2DP_H
