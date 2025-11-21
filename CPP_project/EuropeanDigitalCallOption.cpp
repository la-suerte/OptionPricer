#include "EuropeanDigitalCallOption.h"

EuropeanDigitalCallOption::EuropeanDigitalCallOption(double e, double s)
    : EuropeanDigitalOption(e, s) {}

double EuropeanDigitalCallOption::payoff(double z) {
    if(z >= _strike)
        return 1.0;
    else
        return 0;
}

optionType EuropeanDigitalCallOption::GetOptionType() const {
    return call;
}