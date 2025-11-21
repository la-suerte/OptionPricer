#include "CallOption.h"

CallOption::CallOption(double e, double s) : EuropeanVanillaOption(e, s) {}

double CallOption::payoff(double z) {
    if(z > _strike)
        return z - _strike;
    else
        return 0;
}

optionType CallOption::GetOptionType() const {
    return call;
}