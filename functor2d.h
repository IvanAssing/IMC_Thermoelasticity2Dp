#ifndef FUNCTOR2D_H
#define FUNCTOR2D_H


#include "imc_dfm.h"

// Objeto Funcão (Super Classe Abstrata)
class Functor2D
{
    public:
        Functor2D(){}
        virtual tFloat operator()(tFloat x, tFloat y) = 0;
};

class Constant2D : public Functor2D
{
    public:
        tFloat c;

        Constant2D(tFloat constant){ c = constant;}
        tFloat operator()(tFloat x, tFloat y)
        {
            return c;
        }
};


#endif // FUNCTOR2D_H
