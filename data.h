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
            k = 0.0q;
            heatSource = new Constant2D(0.0q);
        }
};

class Thermoelasticity2DData
{
    public:
        tFloat k, mi, alpha, E, T0;
        Functor2D *heatSource;


        Thermoelasticity2DData(){
            k = mi = alpha = E = T0 = 0.0q;
            heatSource = new Constant2D(0.0q);
        }
};


enum DirectionType{
    South,
    North,
    East,
    West
};


#endif // DATA_H
