#include <iostream>

#include <QApplication>

#include "imc_dfm.h"
#include "diffusion2dp.h"
#include "gaussseidel.h"

#include "graphics.h"

#include "thermoelasticity2dp.h"

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);

    // MANAGE_EXCEPTIONS

    //Diffusion2DData data;
    Thermoelasticity2DData data;

    data.heatSource = new Constant2D(0.0q); // Sem geração de calor
    data.k = 1.0q;
    data.alpha = 0.1q;
    data.mi = 0.3q;
    data.E = 1.0q;



    Constant2D cc0(0.0);
    Constant2D cc1(1.0);
    Constant2D cc2(-1.0);

    Sine ccS(1.0q, M_PIq, 0.0q);

    // Condições de contorno do problema 7.1
    Boundary2D south(Dirichlet, &cc0);
    Boundary2D north(Dirichlet, &ccS);
    Boundary2D east(Dirichlet, &cc0);
    Boundary2D west(Dirichlet, &cc0);

    Boundary2D ccU[4], ccV[4], ccT[4];


    ccU[0] = ccU[1] = ccU[2] = ccU[3] = Boundary2D(Dirichlet, &cc0);
    ccV[0] = ccV[1] = ccV[2] = ccV[3] = Boundary2D(Dirichlet, &cc0);
    ccT[0] = ccT[2] = ccT[3] = Boundary2D(Dirichlet, &cc0);
    ccT[1] = north;

    tFloat lx = 1.0q;
    tFloat ly = 1.0q;
    tInteger nx = 41;
    tInteger ny = 41;

    tFloat TmAS = 2.0q*lx * (coshq(M_PIq*ly/lx) - 1.0q) / (M_PIq*M_PIq*ly*sinhq(M_PIq*ly/lx));
    //Diffusion2Dp mesh(lx, ly, nx, ny, &data, &south, &north, &east, &west);


    Thermoelasticity2Dp mesh(lx, ly, nx, ny, &data, ccT, ccU, ccV);

    //mesh.solver(100000, 1.0e-28q, true); // itmax, itol, plotlog?

    mesh.solver(100000, 1.0e-28q, 20);

    SFAS as(1.0q, 1.0q);

//    mesh.plotX(ny/2, as);
//    mesh.plotY(nx/2, as);
//    mesh.plot(as);

    tFloat *erro = new tFloat[mesh.nx*mesh.ny];

    // Erro númerico
    for(int i=0; i<mesh.nx*mesh.ny; i++)
        erro[i] = fabsq(as(mesh.nodes[i].x, mesh.nodes[i].y)-mesh.T[i]);

    //    // Gerar mapa de cores
    //    Graphics w1/*, w2*/;
    //    w1.mesh = &mesh;
    //    w1.X = mesh.T;
    //    w1.setWindowTitle(QString("Temperatura"));
    ////    w2.mesh = &mesh;
    ////    w2.X = erro;
    //    w1.show();
    //    //w2.show();

        //Graphics u, v, ex, ey, exy;

    //    u.mesh = &mesh;
    //    u.X = mesh.U;
    //    u.setWindowTitle(QString("Deslocamentos na direção X"));

    //    v.mesh = &mesh;
    //    v.X = mesh.V;
    //    v.setWindowTitle(QString("Deslocamentos na direção Y"));

    //    u.show();
    //    v.show();


//    ex.mesh = &mesh;
//    ex.X = mesh.ex;
//    ex.setWindowTitle(QString("Deformação da direção X"));
//    ex.show();


//    Graphics w_t(&mesh, mesh.T, QString("Temperatura"), 0);

//    w_t.show();

//    Graphics w_u(&mesh, mesh.U, QString("Deslocamentos na direção X"), 0);

//    w_u.show();

//    Graphics w_v(&mesh, mesh.V, QString("Deslocamentos na direção Y"), 0);

//    w_v.show();


//    Graphics w_ex(&mesh, mesh.ex, QString("Deformação na direção X"), 0);

//    w_ex.show();

//    Graphics w_ey(&mesh, mesh.ey, QString("Deformação na direção Y"), 0);

//    w_ey.show();

//    Graphics w_exy(&mesh, mesh.exy, QString("Distorção XY"), 0);

//    w_exy.show();

    Graphics w_sx(&mesh, mesh.sx, QString("Tensão normal na direção X"), 0);

    w_sx.show();

    Graphics w_sy(&mesh, mesh.sy, QString("Tensão normal na direção Y"), 0);

    w_sy.show();

    Graphics w_sxy(&mesh, mesh.sxy, QString("Tensão de cisalhamento XY"), 0);

    w_sxy.show();

    std::cout<<"\n\nX = 0.5";
    mesh.printX(ny/2);

    std::cout<<"\n\nY = 0.5";
    mesh.printY(nx/2);

    mesh.plotX(ny/2);
    mesh.plotY(nx/2);

    mesh.plotIterationLog();



    //END_EXCEPTIONS

    return a.exec();
}
