#include "imc_dfm.h"

std::string print(tFloat value)
{
    char str[1000];
    quadmath_snprintf(str, 1000, Q_FORMAT,value);
    return std::string(str);
}

std::string print2(tFloat value)
{
    char str[1000];
    quadmath_snprintf(str, 1000, "%.3Qe",value);
    return std::string(str);
}

