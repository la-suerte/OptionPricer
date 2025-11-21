#include "AsianCallOption.h"
AsianCallOption::AsianCallOption(std::vector<double> ts, double k) 
    : AsianOption(ts, ts.back(), k) {}

double AsianCallOption::payoff(double underlying) {
    if(underlying > _strike)
        return underlying - _strike;
    else
        return 0;
}