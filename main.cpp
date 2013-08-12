#include <iostream>

#include <QApplication>

#include "imc_dfm.h"
#include "diffusion2dp.h"
#include "gaussseidel.h"

#include "graphics.h"

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);


    //MANAGE_EXCEPTIONS
    //{
    //        tFloat *x = new tFloat[11];


    //        GaussSeidel sys(x, 11, 10000, 1.e-33Q, 3);

    //        for(tInteger i=1;i<10;i++)
    //        {
    //            sys(i, i, 2.0q);
    //            sys(i, i-1, -1.0q);
    //            sys(i, i+1, -1.0q);
    //            sys(i, 0.0q);
    //        }

    //        sys(10,10, 1.0q);
    //        sys(10, 1.0q);

    //       sys.solver();

    //        for(int i=0; i<11;i++)
    //            std::cout<<"\n"<<print(x[i]);


    Diffusion2DData data;

    data.heatSource = new Constant2D(0.0q);
    data.k = 1.0q;

    Constant2D cc0(0.0);
    Constant2D cc1(1.0);
    Constant2D cc2(-1.0);

    Sine ccS(1.0q, M_PIq, 0.0q);

    Boundary2D south(Dirichlet, &cc0);
    Boundary2D north(Dirichlet, &ccS);
    Boundary2D east(Dirichlet, &cc0);
    Boundary2D west(Dirichlet, &cc0);


    //    Boundary2D south(Neumann, &cc1);
    //    Boundary2D north(Neumann, &cc1);
    //    Boundary2D east(Dirichlet, &cc0);
    //    Boundary2D west(Dirichlet, &cc0);

    tFloat lx = 1.0q;
    tFloat ly = 1.0q;
    tInteger nx = 21;
     tInteger ny = 21;

    Diffusion2Dp mesh(lx, ly, nx, ny, &data, &south, &north, &east, &west);

    mesh.solver(100000, 1.0e-28q, true);

    SFAS as(1.0q, 1.0q);

    //mesh.plotX(5, as);
    //mesh.plotY(5, as);



    tFloat *erro = new tFloat[mesh.nx*mesh.ny];

    for(int i=0; i<mesh.nx*mesh.ny; i++)
        erro[i] = fabsq(as(mesh.nodes[i].x, mesh.nodes[i].y)-mesh.T[i]);

    Graphics w;
    w.mesh = &mesh;
    w.X = mesh.T;
    //w.show();

    tFloat TmAS = 2.0q*lx * (coshq(M_PIq*ly/lx) - 1.0q) / (M_PIq*M_PIq*ly*sinhq(M_PIq*ly/lx));

    std::cout<<"\n\nTemperatura Média"<<std::setfill(' ');
    std::cout<<std::setw(15)<<std::right<<"\nNumérica: "<<print(mesh.Tm);
    std::cout<<std::setw(15)<<std::right<<"\nAnalítica: "<<print(TmAS);
    std::cout<<std::setw(15)<<std::right<<"\nErro: "<<print(TmAS - mesh.Tm)<<std::endl;


    //}
    //END_EXCEPTIONS

    //                Graphics w;
    //                //w.mesh = &mesh;
    //                w.show();

    return a.exec();
}


//int main()
//{
//    MANAGE_EXCEPTIONS
//    {
//        //        tFloat *x = new tFloat[11];


//        //        GaussSeidel sys(x, 11, 10000, 1.e-33Q, 3);

//        //        for(tInteger i=1;i<10;i++)
//        //        {
//        //            sys(i, i, 2.0q);
//        //            sys(i, i-1, -1.0q);
//        //            sys(i, i+1, -1.0q);
//        //            sys(i, 0.0q);
//        //        }

//        //        sys(10,10, 1.0q);
//        //        sys(10, 1.0q);

//        //       sys.solver();

//        //        for(int i=0; i<11;i++)
//        //            std::cout<<"\n"<<print(x[i]);


//        Diffusion2DData data;

//        data.heatSource = new Constant2D(0.0q);
//        data.k = 1.0q;

//        Constant2D cc(0.0);
//        Sine cc1(1.0q, M_PIq, 0.0q);

//        Boundary2D south(Dirichlet, &cc);
//        Boundary2D north(Dirichlet, &cc1);
//        Boundary2D east(Dirichlet, &cc);
//        Boundary2D west(Dirichlet, &cc);


//        Diffusion2Dp mesh(1.0q, 1.0, 11, 11, &data, &south, &north, &east, &west);

//        mesh.solver();

//        //QApplication a1();
//        Graphics a(this);
//        a.show();

//    }
//    END_EXCEPTIONS

//            return 0;
//}

