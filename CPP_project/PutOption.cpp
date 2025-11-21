#include "PutOption.h"

PutOption::PutOption(double e, double s) : EuropeanVanillaOption(e, s) {}

double PutOption::payoff(double z) {
    if(_strike > z)
        return _strike - z;
    else
        return 0;
}

optionType PutOption::GetOptionType() const {
    return put;
}