#include "gaussseidel.h"

GaussSeidel::GaussSeidel(tFloat *results, tInteger equationMax, tInteger iterationMax, tFloat iterationTolerance, tInteger bandWidth)
    :x(results), neqmax(equationMax), itol(iterationTolerance), imax(iterationMax), bwidth(bandWidth)
{
    neIndex = new tInteger[neqmax];
    for(tInteger i = 0; i<neqmax; i++)
        neIndex[i] = 0;

    neq = 0;
    nit = 0;


    AIndex = new tInteger*[neqmax];
    for(tInteger i = 0; i<neqmax; i++)
        AIndex[i] = new tInteger[bwidth-1];


    ax = new tFloat[neqmax];
    b = new tFloat[neqmax];
    L = new tFloat[imax];

    for(tInteger i = 0; i<neqmax; i++)
        ax[i] = 1.0q, b[i] = 0.0q;


    for(tInteger i = 0; i<imax; i++)
        L[i] = 0.0q;

    A = new tFloat*[neqmax];
    for(tInteger i = 0; i<neqmax; i++)
        A[i] = new tFloat[bwidth-1];


}

void GaussSeidel::operator()(tInteger i, tInteger j, tFloat value)
{
    if(i<0 || i>neqmax || j<0 || j>neqmax)
        throw std::string("error: invalid index");
    else if(i==j)
        ax[i] = value;
    else{
        A[i][neIndex[i]] = value;
        AIndex[i][neIndex[i]] = j;
        neIndex[i]++;
    }
}

void GaussSeidel::operator()(tInteger i, tFloat _b)
{
    if(i<0 || i>neqmax)
        throw std::string("error: invalid index");
    else
        b[i] = _b;
}

void GaussSeidel::solver()
{
    tFloat sum, residual;

    do{
        // Solver
        for(tInteger i = 0; i<neqmax; i++){
            sum = b[i];
            for(tInteger j=0; j<neIndex[i]; j++)
                sum -= A[i][j]*x[AIndex[i][j]];
            x[i] = sum/ax[i];
        }

        // Residuo
        residual = 0.0q;
        for(tInteger i = 0; i<neqmax; i++){
            sum = b[i]-ax[i]*x[i];
            for(tInteger j=0; j<neIndex[i]; j++)
                sum -= A[i][j]*x[AIndex[i][j]];
            residual += sum*sum;
        }
        L[nit] = sqrtq(residual);

        std::cout<<"\n"<<nit<<"\t"<<print(L[nit]);

    }while(L[nit] > itol && nit++ < imax);
}


tFloat GaussSeidel::residual(tFloat *_x)
{
    // Residuo
    tFloat residual = 0.0q;
    for(tInteger i = 0; i<neqmax; i++){
        tFloat sum = b[i]-ax[i]*_x[i];
        for(tInteger j=0; j<neIndex[i]; j++)
            sum -= A[i][j]*_x[AIndex[i][j]];
        residual += sum*sum;
    }

    return sqrtq(residual);
}


void GaussSeidel::plotIterationLog()
{
    const std::string cmd_filename = "plotconfig_it.gnu";
    const std::string pic_filename = "iteractive_it.png";
    const std::string dat1_filename = "data_it.txt";

    // Solução numérica
    std::ofstream file1(dat1_filename.c_str());
    for(tInteger i=0; i<nit; i++)
        file1<<i<<"\t"<<QtoD(L[i]/L[0])<<std::endl;
    file1.close();

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

             "plot '" <<dat1_filename<<"' t\"\" with lines lt 2 lc 1 lw 2";
    file3.close();

    const std::string cmd1 = "gnuplot " + cmd_filename; // Gráfico com GNUPLOT
    const std::string cmd2 = "eog " + pic_filename; // Visualizador de imagem

    std::system(cmd1.c_str());
    std::system(cmd2.c_str());
}
