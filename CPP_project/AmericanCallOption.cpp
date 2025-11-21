#include "AmericanCallOption.h"
AmericanCallOption::AmericanCallOption(double e, double s) : AmericanOption(e, s) {}

double AmericanCallOption::payoff(double St) {
    if(St > _strike)
        return St - _strike;
    else 
        return 0.0;
}

optionType AmericanCallOption::GetOptionType() {
    return call;
}