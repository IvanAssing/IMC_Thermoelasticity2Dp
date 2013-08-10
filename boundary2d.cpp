#include "boundary2d.h"

Boundary2D::Boundary2D()
{
}

Boundary2D::Boundary2D(BoundaryCondition _type, Functor2D *_bcValue)
{
    type = _type;
    bcValue = _bcValue;
}

