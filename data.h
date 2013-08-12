#ifndef DATA_H
#define DATA_H

#include "imc_dfm.h"
#include "functor2d.h"

class Diffusion2DData
{
    public:
        tFloat k;
        Functor2D *heatSource;


        Diffusion2DData(){
            k = 0.0;
            heatSource = new Constant2D(0.0);
        }
};


enum DirectionType{
    South,
    North,
    East,
    West
};


#endif // DATA_H
