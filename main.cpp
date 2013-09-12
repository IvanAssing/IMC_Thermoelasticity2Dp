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

    Thermoelasticity2DData data;

    // Dados
    data.heatSource = new Constant2D(0.0q); // Sem geração de calor
    data.k = 1.0q;
    data.alpha = 0.1q;
    data.mi = 0.3q;
    data.E = 1.0q;

    // Condições de contorno do problema
    Constant2D cc0(0.0);
    Sine ccS(1.0q, M_PIq, 0.0q);
    Boundary2D north(Dirichlet, &ccS);
    Boundary2D ccU[4], ccV[4], ccT[4];

    ccU[0] = ccU[1] = ccU[2] = ccU[3] = Boundary2D(Dirichlet, &cc0);
    ccV[0] = ccV[1] = ccV[2] = ccV[3] = Boundary2D(Dirichlet, &cc0);
    ccT[0] = ccT[2] = ccT[3] = Boundary2D(Dirichlet, &cc0);
    ccT[1] = north;

    // Parametros da malha
    tFloat lx = 1.0q;
    tFloat ly = 1.0q;
    tInteger nx = 11;
    tInteger ny = 11;

    // Malha
    Thermoelasticity2Dp mesh(lx, ly, nx, ny, &data, ccT, ccU, ccV);

    // Solver
    mesh.solver(1000000, 1.0e-28q, 10);


    // Resultados gráficos (Opengl)
    Graphics w_t(&mesh, mesh.T, QString("Temperatura"), 0);
    w_t.show();

    Graphics w_u(&mesh, mesh.U, QString("Deslocamentos na direção X"), 0);
    w_u.show();

    Graphics w_v(&mesh, mesh.V, QString("Deslocamentos na direção Y"), 0);
    w_v.show();

    Graphics w_ex(&mesh, mesh.ex, QString("Deformação na direção X"), 0);
    w_ex.show();

    Graphics w_ey(&mesh, mesh.ey, QString("Deformação na direção Y"), 0);
    w_ey.show();

    Graphics w_exy(&mesh, mesh.exy, QString("Distorção XY"), 0);
    w_exy.show();

    Graphics w_sx(&mesh, mesh.sx, QString("Tensão normal na direção X"), 0);
    w_sx.show();

    Graphics w_sy(&mesh, mesh.sy, QString("Tensão normal na direção Y"), 0);
    w_sy.show();

    Graphics w_sxy(&mesh, mesh.sxy, QString("Tensão de cisalhamento XY"), 0);
    w_sxy.show();

    // Resultados numéricos
    std::cout<<"\n\n - Deslocamentos";
    mesh.printX(ny/2);
    mesh.printY(nx/2);

    std::cout<<"\n\n - Deformações";
    mesh.printXStrain(ny/2);
    mesh.printYStrain(nx/2);

    std::cout<<"\n\n - Tensões";
    mesh.printXStress(ny/2);
    mesh.printYStress(nx/2);


    // Resultados gráficos
    mesh.plotX(ny/2);
    mesh.plotY(nx/2);
    mesh.plotIterationLog();

    //END_EXCEPTIONS

    return a.exec();
}
