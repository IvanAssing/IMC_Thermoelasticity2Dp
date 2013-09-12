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

    imax = iterationMax;
    itol = iterationTolerance;

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

        std::cout<<"\n Lu["<<nit<<"] = "<<print(Lu[nit])<<"    Lv["<<nit<<"] = "<<print(Lv[nit]);

    }while(Lu[nit] > iterationTolerance && Lv[nit] > iterationTolerance && nit++ < iterationMax);

    // RESULTADOS SECUNDÁRIOS
    ex = new tFloat[nx*ny];
    ey = new tFloat[nx*ny];
    exy = new tFloat[nx*ny];
    sx = new tFloat[nx*ny];
    sy = new tFloat[nx*ny];
    sxy = new tFloat[nx*ny];

    tFloat *dudx = new tFloat[nx*ny];
    tFloat *dudy = new tFloat[nx*ny];
    tFloat *dvdx = new tFloat[nx*ny];
    tFloat *dvdy = new tFloat[nx*ny];

    // Derivadas
    for(tInteger j=1; j<ny-1; j++)
        for(tInteger i=1; i<nx-1; i++)
        {
            tInteger p = position(i,j);
            dudx[p] = (U[direction(p,East)] - U[direction(p,West)])/(2.0q*hx);
            dvdx[p] = (V[direction(p,East)] - V[direction(p,West)])/(2.0q*hx);
            dudy[p] = (U[direction(p,North)] - U[direction(p,South)])/(2.0q*hy);
            dvdy[p] = (V[direction(p,North)] - V[direction(p,South)])/(2.0q*hy);
        }

    // Sul
    for(tInteger i=1; i<nx; i++)
    {
        dudx[i] = (U[direction(i,East)] - U[direction(i,West)])/(2.0q*hx);
        dvdx[i] = (V[direction(i,East)] - V[direction(i,West)])/(2.0q*hx);
        dvdy[i] = (4.0q*V[direction(i,North)] - 3.0q*V[i] - V[direction(direction(i,North),North)])/(2.0q*hy);
        dudy[i] = (4.0q*U[direction(i,North)] - 3.0q*U[i] - U[direction(direction(i,North),North)])/(2.0q*hy);
    }

    // Norte
    for(tInteger i=(ny-1)*nx; i<nx*ny; i++)
    {
        dudx[i] = (U[direction(i,East)] - U[direction(i,West)])/(2.0q*hx);
        dvdx[i] = (V[direction(i,East)] - V[direction(i,West)])/(2.0q*hx);
        dvdy[i] = (-4.0q*V[direction(i,South)] + 3.0q*V[i] + V[direction(direction(i,South),South)])/(2.0q*hy);
        dudy[i] = (-4.0q*U[direction(i,South)] + 3.0q*U[i] + U[direction(direction(i,South),South)])/(2.0q*hy);
    }

    // Leste
    for(tInteger ii=2; ii<ny; ii++)
    {
        tInteger i = ii*nx-1;
        dudx[i] = (4.0q*U[direction(i,East)] - 3.0q*U[i] - U[direction(direction(i,East),East)])/(2.0q*hx);
        dvdx[i] = (4.0q*V[direction(i,East)] - 3.0q*V[i] - V[direction(direction(i,East),East)])/(2.0q*hx);
        dudy[i] = (U[direction(i,North)] - U[direction(i,South)])/(2.0q*hy);
        dvdy[i] = (V[direction(i,North)] - V[direction(i,South)])/(2.0q*hy);
    }

    // Oeste
    for(tInteger ii=1; ii<ny; ii++)
    {
        tInteger i = ii*nx;
        dudx[i] = (-4.0q*U[direction(i,West)] + 3.0q*U[i] + U[direction(direction(i,West),West)])/(2.0q*hx);
        dvdx[i] = (-4.0q*V[direction(i,West)] + 3.0q*V[i] + V[direction(direction(i,West),West)])/(2.0q*hx);
        dudy[i] = (U[direction(i,North)] - U[direction(i,South)])/(2.0q*hy);
        dvdy[i] = (V[direction(i,North)] - V[direction(i,South)])/(2.0q*hy);
    }

    for(tInteger p=0; p<nx*ny; p++)
    {
        ex[p] = dudx[p];
        ey[p] = dvdy[p];
        exy[p] = (dudy[p] + dvdx[p])/2.0q;

        sx[p] = data->E/(1.0q+data->mi)*(
                    data->mi/(1.0q-data->mi)*(dudx[p]+dvdy[p]) + dudx[p]
                    - Cmi*data->alpha*(T[p]-data->T0));
        sy[p] = data->E/(1.0q+data->mi)*(
                    data->mi/(1.0q-data->mi)*(dudx[p]+dvdy[p]) + dvdy[p]
                    - Cmi*data->alpha*(T[p]-data->T0));

        sxy[p] = 0.5q*data->E/(1.0q+data->mi)*(dudy[p] + dvdx[p]);
    }

    delete dudx;
    delete dudy;
    delete dvdx;
    delete dvdy;
}


void Thermoelasticity2Dp::printX(tInteger p)
{

    for(tInteger i=0; i<nx; i++)
        std::cout<<std::endl<<p*nx+i<<"\t"<<print2(nodes[p*nx+i].x)<<"\t"<<print2(nodes[p*nx+i].y)<<"\t"<<print(U[p*nx+i])<<"\t"<<print(V[p*nx+i])<<"\t"<<print(T[p*nx+i]);


}

void Thermoelasticity2Dp::printY(tInteger p)
{

    for(tInteger i=0; i<ny; i++)
        std::cout<<std::endl<<p+nx*i<<"\t"<<print2(nodes[p+nx*i].x)<<"\t"<<print2(nodes[p+nx*i].y)<<"\t"<<print(U[p+nx*i])<<"\t"<<print(V[p+nx*i])<<"\t"<<print(T[p+nx*i]);
}


void Thermoelasticity2Dp::plotY(tInteger p)
{
    const std::string cmd_filename = "plotconfig.gnu";
    const std::string pic_filename = "ploty.png";
    const std::string dat1_filename = "data1.txt";
    const std::string dat2_filename = "data2.txt";

    // Solução numérica
    std::ofstream file1(dat1_filename.c_str());
    for(tInteger i=0; i<ny; i++)
        file1<<QtoD(nodes[p+nx*i].y)<<"\t"<<QtoD(U[p+nx*i])<<std::endl;
    file1.close();

    std::ofstream file2(dat2_filename.c_str());
    for(tInteger i=0; i<ny; i++)
        file2<<QtoD(nodes[p+nx*i].y)<<"\t"<<QtoD(V[p+nx*i])<<std::endl;
    file2.close();


    std::ofstream file3(cmd_filename.c_str());
    file3 <<
             "set terminal pngcairo enhanced font \"arial,12\" size 1600, 1000 \n"
             "set output '" << pic_filename <<"'\n"
             "set key inside right top vertical Right noreverse enhanced autotitles box linetype -1 linewidth 1.000\n"
             "set grid\n"
             "set title \"TERMOELASTICIDADE BIDIMENSIONAL EM REGIME PERMANENTE\\nResolução com MDF / Aproximação com CDS-2\"\n"
             "set xlabel 'y'\n"
             "set ylabel \"u( "<<QtoD(p*hx)<<", y), v( "<<QtoD(p*hx)<<", y)\"\n"
             "plot '" <<dat2_filename<<"' t\"u( "<<QtoD(p*hx)<<", y)\" with linespoints lt 2 lc 2 pt 4 lw 1, "
             "'" <<dat1_filename<<"' t\"v( "<<QtoD(p*hx)<<", y)\" with linespoints lt 2 lc 1 pt 3 lw 1";
    file3.close();

    const std::string cmd1 = "gnuplot " + cmd_filename; // Gráfico com GNUPLOT
    const std::string cmd2 = "eog " + pic_filename; // Visualizador de imagem

    std::system(cmd1.c_str());
    std::system(cmd2.c_str());
}

void Thermoelasticity2Dp::plotX(tInteger p)
{
    const std::string cmd_filename = "plotconfig.gnu";
    const std::string pic_filename = "plotx.png";
    const std::string dat1_filename = "data1.txt";
    const std::string dat2_filename = "data2.txt";

    // Solução numérica
    std::ofstream file1(dat1_filename.c_str());
    for(tInteger i=0; i<nx; i++)
        file1<<QtoD(nodes[p*nx+i].x)<<"\t"<<QtoD(U[p*nx+i])<<std::endl;
    file1.close();

    std::ofstream file2(dat2_filename.c_str());
    for(tInteger i=0; i<nx; i++)
        file2<<QtoD(nodes[p*nx+i].x)<<"\t"<<QtoD(V[p*nx+i])<<std::endl;
    file2.close();


    std::ofstream file3(cmd_filename.c_str());
    file3 <<
             "set terminal pngcairo enhanced font \"arial,12\" size 1600, 1000 \n"
             "set output '" << pic_filename <<"'\n"
             "set key inside right top vertical Right noreverse enhanced autotitles box linetype -1 linewidth 1.000\n"
             "set grid\n"
             "set title \"TERMOELASTICIDADE BIDIMENSIONAL EM REGIME PERMANENTE\\nResolução com MDF / Aproximação com CDS-2\"\n"
             "set xlabel 'x'\n"
             "set ylabel \"u( x, "<<QtoD(p*hy)<<"), v( x, "<<QtoD(p*hy)<<")\"\n"
             "plot '" <<dat2_filename<<"' t\"u( x, "<<QtoD(p*hy)<<")\" with linespoints lt 2 lc 2 pt 4 lw 1, "
             "'" <<dat1_filename<<"' t\"v( x, "<<QtoD(p*hy)<<")\" with linespoints lt 2 lc 1 pt 3 lw 1";
    file3.close();

    const std::string cmd1 = "gnuplot " + cmd_filename; // Gráfico com GNUPLOT
    const std::string cmd2 = "eog " + pic_filename; // Visualizador de imagem

    std::system(cmd1.c_str());
    std::system(cmd2.c_str());
}


void Thermoelasticity2Dp::plotIterationLog()
{
    const std::string cmd_filename = "plotconfig_it.gnu";
    const std::string pic_filename = "iteractive_it.png";
    const std::string dat1_filename = "data_itu.txt";
        const std::string dat2_filename = "data_itv.txt";

    std::ofstream file1(dat1_filename.c_str());
    for(tInteger i=0; i<nit; i++)
        file1<<i<<"\t"<<QtoD(Lu[i]/Lu[0])<<std::endl;
    file1.close();

    std::ofstream file2(dat2_filename.c_str());
    for(tInteger i=0; i<nit; i++)
        file2<<i<<"\t"<<QtoD(Lv[i]/Lv[0])<<std::endl;
    file2.close();

    std::ofstream file3(cmd_filename.c_str());
    file3 <<
             "set terminal pngcairo enhanced font \"arial,12\" size 1600, 1000 \n"
             "set output '" << pic_filename <<"'\n"
             "#set key inside right top vertical Right noreverse enhanced autotitles box linetype -1 linewidth 1.000\n"
             "set grid\n"
             "set logscale y\n"
             "set format y \"10^{%L}\" \n"
             "set lmargin 10 \n"
             "set title \"DESEMPENHO DE ITERAÇÃO \\n itmax = "<<imax<<"  itol = "<<QtoD(itol)<<"  nit = "<<nit-1<<"\"\n"
             "set ylabel \"L^{n}/L^{0}\" \n"
             "set xlabel 'Número de iterações'\n"

             "plot '" <<dat1_filename<<"' t\"Lu\" with lines lt 2 lc 1 lw 2, "
             "'" <<dat2_filename<<"' t\"Lv\" with lines lt 2 lc 2 lw 2";

    file3.close();

    const std::string cmd1 = "gnuplot " + cmd_filename; // Gráfico com GNUPLOT
    const std::string cmd2 = "eog " + pic_filename; // Visualizador de imagem

    std::system(cmd1.c_str());
    std::system(cmd2.c_str());
}


void Thermoelasticity2Dp::printXStrain(tInteger p)
{

    for(tInteger i=0; i<nx; i++)
        std::cout<<std::endl<<p*nx+i<<"\t"<<print2(nodes[p*nx+i].x)<<"\t"<<print2(nodes[p*nx+i].y)<<"\t"<<print(ex[p*nx+i])<<"\t"<<print(ey[p*nx+i])<<"\t"<<print(exy[p*nx+i]);



}

void Thermoelasticity2Dp::printYStrain(tInteger p)
{

    for(tInteger i=0; i<ny; i++)
        std::cout<<std::endl<<p+nx*i<<"\t"<<print2(nodes[p+nx*i].x)<<"\t"<<print2(nodes[p+nx*i].y)<<"\t"<<print(ex[p+nx*i])<<"\t"<<print(ey[p+nx*i])<<"\t"<<print(exy[p+nx*i]);
}



void Thermoelasticity2Dp::printXStress(tInteger p)
{

    for(tInteger i=0; i<nx; i++)
        std::cout<<std::endl<<p*nx+i<<"\t"<<print2(nodes[p*nx+i].x)<<"\t"<<print2(nodes[p*nx+i].y)<<"\t"<<print(sx[p*nx+i])<<"\t"<<print(sy[p*nx+i])<<"\t"<<print(sxy[p*nx+i]);



}

void Thermoelasticity2Dp::printYStress(tInteger p)
{

    for(tInteger i=0; i<ny; i++)
        std::cout<<std::endl<<p+nx*i<<"\t"<<print2(nodes[p+nx*i].x)<<"\t"<<print2(nodes[p+nx*i].y)<<"\t"<<print(sx[p+nx*i])<<"\t"<<print(sy[p+nx*i])<<"\t"<<print(sxy[p+nx*i]);
}



