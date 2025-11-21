#include "Option.h"

Option::Option(double e) : _expiry(e) {}

double Option::getExpiry() {
    return _expiry;
}

double Option::payoffPath(std::vector<double> d) {
    return payoff(d[d.size()-1]);
}

bool Option::isAsianOption() {
    return false;
}

bool Option::isAmericanOption() {
    return false;
}