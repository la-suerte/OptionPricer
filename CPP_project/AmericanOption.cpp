#include "AmericanOption.h"
AmericanOption::AmericanOption(double e, double s) : Option(e), _strike(s) {}

bool AmericanOption::isAmericanOption() {
    return true;
}