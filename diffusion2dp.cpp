#include "diffusion2dp.h"

Diffusion2Dp::Diffusion2Dp(tFloat _lx, tFloat _ly, tInteger _nx, tInteger _ny, Diffusion2DData *_data,
                                   Boundary2D _bS, Boundary2D _bN, Boundary2D _bE, Boundary2D _bW)
    :lx(_lx), ly(_ly), nx(_nx), ny(_ny), data(_data), bS(_bS), bN(_bN), bE(_bE), bW(_bW)
{
    hx = lx/(nx - 1.0q);
    hy = ly/(ny - 1.0q);

    nodes = new Node2D[nx*ny]();
    tInteger p = 0;

    for(tInteger j=0; j<ny; j++)
        for(tInteger i=0; i<nx; i++)
            nodes[p] = Node2D(i*hx, j*hy, i, j, p++);

}
