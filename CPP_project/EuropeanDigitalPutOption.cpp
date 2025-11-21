#include "EuropeanDigitalPutOption.h"
EuropeanDigitalPutOption::EuropeanDigitalPutOption(double exp, double s) 
    : EuropeanDigitalOption(exp, s) {}

double EuropeanDigitalPutOption::payoff(double z) {
    if(z <= _strike)
        return 1.0;
    else
        return 0;
}

optionType EuropeanDigitalPutOption::GetOptionType() const {
    return put;
}