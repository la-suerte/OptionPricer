#include "EuropeanVanillaOption.h"
#include "EuropeanVanillaOption.h"

EuropeanVanillaOption::EuropeanVanillaOption(double e, double s) 
    : Option(e >= 0 ? e : 0), _strike(s >= 0 ? s : 0) {}

