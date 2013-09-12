#ifndef IMC_DFM_H
#define IMC_DFM_H

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <string>
#include <fstream>
#include <cmath>

extern "C" {
#include "quadmath.h" // Biblioteca para Precisão Quadrúpla
}

#define QUAD_PRECISION // Ativar para habilitar Precisão Quadrúpla

#ifdef QUAD_PRECISION
typedef __float128 tFloat; // Tipo Ponto flutuante com precisão Quadrúpla (128 bits)
#else
typedef double tFloat; // Tipo Ponto flutuante com precisão Dupla (64 bits)
//typedef float tFloat; // Tipo Ponto flutuante com precisão Simples (32 bits)
#endif

typedef long int tInteger; // Tipo Inteiro (64 bits)

#define TOFLOAT(x) static_cast<tFloat>(x)

#ifdef QUAD_PRECISION
#define OUT_FLOAT_PRECISION 33 // Largura do campo de impressão
#else
#define OUT_FLOAT_PRECISION 14
#endif

#define Q_FORMAT "%.33Qe" // Mostrar número em precisão quadrúpla com 33 casas decimais

#define OUT_FLOAT_WIDTH (OUT_FLOAT_PRECISION+10)
#define OUT_TXT 10

#define MANAGE_EXCEPTIONS try {
#define END_EXCEPTIONS     }catch(std::string str){std::cout<<str;}


#define QtoD(value) static_cast<double>(value)

std::string print(tFloat value); // 33 casas
std::string print2(tFloat value); // 3 casas


#endif // IMC_DFM_H
