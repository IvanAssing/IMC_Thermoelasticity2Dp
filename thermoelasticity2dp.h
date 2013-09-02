#ifndef THERMOELASTICITY2DP_H
#define THERMOELASTICITY2DP_H

#include "imc_dfm.h"
#include "functor2d.h"
#include "node2d.h"
#include "boundary2d.h"
#include "data.h"
#include "diffusion2dp.h"

class Thermoelasticity2Dp
{

    public:
        Thermoelasticity2DData *data; // Dados do problema
        Node2D *nodes; // Lista de nós
        tInteger nx, ny; // Número de divisões em x e y
        tFloat lx, ly, hx, hy; // Tamanho de malha
        Boundary2D *ccT, *ccU, *ccV; // Condições de contorno

        tInteger nit;

        tFloat *Lu, *Lv; // Resíduos

        tFloat *T, *U, *V; // Vetores solução

        Diffusion2Dp *thermo;

        Thermoelasticity2Dp(tFloat lengthX, tFloat lengthY, tInteger nx, tInteger ny, Thermoelasticity2DData *data,
                            Boundary2D *ccT, Boundary2D *ccU, Boundary2D *ccV);

        tInteger direction(tInteger position, DirectionType dir);
        tInteger position(tInteger i, tInteger j);

        void solver(tInteger iterationMax = 1000, tFloat iterationTolerance = 1.0e-28q, tInteger internalIterationMax = 10, bool plotlog = false); // Discretização + solver do sistema linear

        void printX(tInteger p);
        void printY(tInteger p);
};

#endif // THERMOELASTICITY2DP_H
