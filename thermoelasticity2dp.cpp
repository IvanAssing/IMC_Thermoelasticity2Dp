#include "thermoelasticity2dp.h"

#include "gaussseidel.h"

Thermoelasticity2Dp::Thermoelasticity2Dp(tFloat _lx, tFloat _ly, tInteger _nx, tInteger _ny,
                                         Thermoelasticity2DData *_data, Boundary2D *_ccT, Boundary2D *_ccU, Boundary2D *_ccV)
    :lx(_lx), ly(_ly), nx(_nx), ny(_ny), data(_data), ccT(_ccT), ccU(_ccU), ccV(_ccV)
{
    Diffusion2DData data_thermo;
    data_thermo.heatSource = data->heatSource;
    data_thermo.k = data->k;

    thermo = new Diffusion2Dp(lx, ly, nx, ny, &data_thermo, &ccT[0], &ccT[1], &ccT[2], &ccT[3]);

    hx = thermo->hx;
    hy = thermo->hy;

    nodes = thermo->nodes;

    T = thermo->T;

    U = new tFloat[nx*ny];
    V = new tFloat[nx*ny];

    for(tInteger i=0; i<nx*ny; i++)
        U[i] = V[i] = 0.q;

}


tInteger Thermoelasticity2Dp::position(tInteger i, tInteger j)
{
    return i+j*nx;
}


tInteger Thermoelasticity2Dp::direction(tInteger position, DirectionType dir)
{
    switch (dir) {
    case North:
        if(position+nx > nx*ny)
            throw std::string("error: invalid index");
        else
            return position+nx;
        break;
    case South:
        if(position-nx < 0)
            throw std::string("error: invalid index");
        else
            return position-nx;
        break;
    case West:
        if(position-1 < 0)
            throw std::string("error: invalid index");
        else
            return position-1;
        break;
    case East:
        if(position+1 > nx*ny)
            throw std::string("error: invalid index");
        else
            return position+1;
        break;
    default:
        throw std::string("error: invalid argument");
        break;
    }
}

void Thermoelasticity2Dp::solver(tInteger iterationMax, tFloat iterationTolerance, tInteger internalIterationMax, bool plotlog)
{

    thermo->solver(iterationMax);

    GaussSeidel sysU(U, nx*ny, internalIterationMax, iterationTolerance);
    GaussSeidel sysV(V, nx*ny, internalIterationMax, iterationTolerance);

    tFloat Cmi = (1.0q + data->mi)/(1.0q - data->mi);
    tFloat a[4];

    for(tInteger j=1; j<ny-1; j++)
        for(tInteger i=1; i<nx-1; i++)
        {
            tInteger p = position(i,j);
            a[0] = a[1] = 1.0q/(hy*hy);
            a[2] = a[3] = (1.0q + Cmi)/(hx*hx);

            sysU(p, p, -(a[0]+a[1]+a[2]+a[3]));

            sysU(p, direction(p, South), a[0]);
            sysU(p, direction(p, North), a[1]);
            sysU(p, direction(p, East), a[2]);
            sysU(p, direction(p, West), a[3]);

//            sysU(p,
//                 Cmi*((V[direction(direction(p, North), East)] - V[direction(direction(p, North), West)] -
//                 V[direction(direction(p, South), East)] + V[direction(direction(p, South), West)])/(4.0q*hx*hy)
//                    - (data->alpha/hx)*(T[direction(p, East)] - T[direction(p, West)])));
        }

    for(tInteger j=1; j<ny-1; j++)
        for(tInteger i=1; i<nx-1; i++)
        {
            tInteger p = position(i,j);
            a[0] = a[1] = (1.0q + Cmi)/(hy*hy);
            a[2] = a[3] = 1.0q/(hx*hx);

            sysV(p, p, -(a[0]+a[1]+a[2]+a[3]));

            sysV(p, direction(p, South), a[0]);
            sysV(p, direction(p, North), a[1]);
            sysV(p, direction(p, East), a[2]);
            sysV(p, direction(p, West), a[3]);

//            sysV(p,
//                 Cmi*((U[direction(direction(p, North), East)] - U[direction(direction(p, North), West)] -
//                 U[direction(direction(p, South), East)] + U[direction(direction(p, South), West)])/(4.0q*hx*hy)
//                    - (data->alpha/hy)*(T[direction(p, North)] - T[direction(p, South)])));
        }

        // Condições de Contorno
        // Sul
        for(tInteger i=0; i<nx; i++)
            if(ccU[0].type == Dirichlet)
                sysU(i, ccU[0].bcValue->operator ()(nodes[i].x, nodes[i].y));
            else{
    //            sys(i, i, -3.0q);
    //            sys(i, direction(i, North), 4.0q);
    //            sys(i, direction(direction(i, North), North), -1.0q);
    //            sys(i, -2.0q*hy*ccS->bcValue->operator ()(nodes[i].x, nodes[i].y)/data->k);
            }

        // Norte
        for(tInteger i=(ny-1)*nx; i<nx*ny; i++)
            if(ccU[1].type == Dirichlet)
                sysU(i, ccU[1].bcValue->operator ()(nodes[i].x, nodes[i].y));
            else{
    //            sys(i, i, 3.0q);
    //            sys(i, direction(i, South), -4.0q);
    //            sys(i, direction(direction(i, South), South), 1.0q);
    //            sys(i, -2.0q*hy*ccN->bcValue->operator ()(nodes[i].x, nodes[i].y)/data->k);
            }
        // East
        for(tInteger i=1; i<ny; i++)
            if(ccU[2].type == Dirichlet)
                sysU(i*nx-1, ccU[2].bcValue->operator ()(nodes[i*nx-1].x, nodes[i*nx-1].y));
            else{
    //            int p = i*nx-1;
    //            sys(p, p, 3.0q);
    //            sys(p, direction(p, West), -4.0q);
    //            sys(p, direction(direction(p, West), West), 1.0q);
    //            sys(p, -2.0q*hx*ccE->bcValue->operator ()(nodes[p].x, nodes[p].y)/data->k);
            }

        // West
        for(tInteger i=0; i<ny; i++)
            if(ccU[3].type == Dirichlet)
                sysU(i*nx, ccU[3].bcValue->operator ()(nodes[i*nx].x, nodes[i*nx].y));
            else{
    //            int p = i*nx;
    //            sys(p, p, -3.0q);
    //            sys(p, direction(p, East), 4.0q);
    //            sys(p, direction(direction(p, East), East), -1.0q);
    //            sys(p, -2.0q*hx*ccW->bcValue->operator ()(nodes[p].x, nodes[p].y)/data->k);
            }

        // Condições de Contorno
        // Sul
        for(tInteger i=0; i<nx; i++)
            if(ccV[0].type == Dirichlet)
                sysV(i, ccV[0].bcValue->operator ()(nodes[i].x, nodes[i].y));
            else{
    //            sys(i, i, -3.0q);
    //            sys(i, direction(i, North), 4.0q);
    //            sys(i, direction(direction(i, North), North), -1.0q);
    //            sys(i, -2.0q*hy*ccS->bcValue->operator ()(nodes[i].x, nodes[i].y)/data->k);
            }

        // Norte
        for(tInteger i=(ny-1)*nx; i<nx*ny; i++)
            if(ccV[1].type == Dirichlet)
                sysV(i, ccV[1].bcValue->operator ()(nodes[i].x, nodes[i].y));
            else{
    //            sys(i, i, 3.0q);
    //            sys(i, direction(i, South), -4.0q);
    //            sys(i, direction(direction(i, South), South), 1.0q);
    //            sys(i, -2.0q*hy*ccN->bcValue->operator ()(nodes[i].x, nodes[i].y)/data->k);
            }
        // East
        for(tInteger i=1; i<ny; i++)
            if(ccV[2].type == Dirichlet)
                sysV(i*nx-1, ccV[2].bcValue->operator ()(nodes[i*nx-1].x, nodes[i*nx-1].y));
            else{
    //            int p = i*nx-1;
    //            sys(p, p, 3.0q);
    //            sys(p, direction(p, West), -4.0q);
    //            sys(p, direction(direction(p, West), West), 1.0q);
    //            sys(p, -2.0q*hx*ccE->bcValue->operator ()(nodes[p].x, nodes[p].y)/data->k);
            }

        // West
        for(tInteger i=0; i<ny; i++)
            if(ccV[3].type == Dirichlet)
                sysV(i*nx, ccV[3].bcValue->operator ()(nodes[i*nx].x, nodes[i*nx].y));
            else{
    //            int p = i*nx;
    //            sys(p, p, -3.0q);
    //            sys(p, direction(p, East), 4.0q);
    //            sys(p, direction(direction(p, East), East), -1.0q);
    //            sys(p, -2.0q*hx*ccW->bcValue->operator ()(nodes[p].x, nodes[p].y)/data->k);
            }


    tFloat sum/*, residual*/;
    nit = 0;
    Lu = new tFloat[iterationMax];
    Lv = new tFloat[iterationMax];

    do{
        // Solver U
        for(tInteger j=1; j<ny-1; j++)
            for(tInteger i=1; i<nx-1; i++)
            {
                tInteger p = position(i,j);
                sysU(p,
                     -Cmi*((V[direction(direction(p, North), East)] - V[direction(direction(p, North), West)] -
                     V[direction(direction(p, South), East)] + V[direction(direction(p, South), West)])/(4.0q*hx*hy)
                        - (data->alpha/hx)*(T[direction(p, East)] - T[direction(p, West)])));

            }
        sysU.solver();


        // Solver V
        for(tInteger j=1; j<ny-1; j++)
            for(tInteger i=1; i<nx-1; i++)
            {
                tInteger p = position(i,j);
                sysV(p,
                     -Cmi*((U[direction(direction(p, North), East)] - U[direction(direction(p, North), West)] -
                     U[direction(direction(p, South), East)] + U[direction(direction(p, South), West)])/(4.0q*hx*hy)
                        - (data->alpha/hy)*(T[direction(p, North)] - T[direction(p, South)])));
            }
        sysV.solver();

        // Residuo
        Lu[nit] = sysU.residual(U);
        Lv[nit] = sysV.residual(V);
//        Lu[nit] = 1.0;
//        Lv[nit] = 1.0;

        std::cout<<"\n"<<nit<<"\t"<<print(Lu[nit])<<"\t"<<print(Lv[nit]);

    }while(Lu[nit] > iterationTolerance && Lv[nit] > iterationTolerance && nit++ < iterationMax);
}


void Thermoelasticity2Dp::printX(tInteger p)
{

    for(tInteger i=0; i<nx; i++)
        std::cout<<std::endl<<print(nodes[p*nx+i].x)<<"\t"<<print(U[p*nx+i])<<"\t"<<print(V[p*nx+i])<<"\t"<<print(T[p*nx+i]);


}

void Thermoelasticity2Dp::printY(tInteger p)
{

    for(tInteger i=0; i<ny; i++)
        std::cout<<std::endl<<print(nodes[p+nx*i].y)<<"\t"<<print(U[p+nx*i])<<"\t"<<print(V[p+nx*i])<<"\t"<<print(T[p+nx*i]);


}
