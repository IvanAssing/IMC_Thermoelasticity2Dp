#ifndef NODE2D_H
#define NODE2D_H

#include "imc_dfm.h"

class Node2D
{
    public:
        tFloat x, y;
        tInteger i,j,p;

        Node2D(){}
        Node2D(tFloat x, tFloat y, tInteger i, tInteger j, tInteger p);
};

#endif // NODE2D_H
