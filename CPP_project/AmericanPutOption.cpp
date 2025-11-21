#include "AmericanPutOption.h"
AmericanPutOption::AmericanPutOption(double e, double s) : AmericanOption(e, s) {}

double AmericanPutOption::payoff(double St) {
    if(St < _strike)
        return _strike - St;
    else 
        return 0.0;
}

optionType AmericanPutOption::GetOptionType() {
    return put;
}