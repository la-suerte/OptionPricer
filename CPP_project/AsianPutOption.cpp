#include "AsianPutOption.h"
AsianPutOption::AsianPutOption(std::vector<double> ts, double k) 
    : AsianOption(ts, ts.back(), k) {}

double AsianPutOption::payoff(double underlying) {
    if(underlying < _strike)
        return _strike - underlying;
    else
        return 0;
}