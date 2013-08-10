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
        Diffusion2DData *data;
        Node2D *nodes;
        tInteger nx, ny;
        tFloat lx, ly, hx, hy;
        Boundary2D bS, bN, bE, bW;


        Diffusion2Dp(tFloat lx, tFloat ly, tInteger nx, tInteger ny, Diffusion2DData *data,
                     Boundary2D bS, Boundary2D bN, Boundary2D bE, Boundary2D bW);
};

#endif // DIFFUSION2DP_H
