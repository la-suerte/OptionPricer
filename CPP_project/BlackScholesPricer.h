#ifndef BLACKSCHOLESPRICER_H
#define BLACKSCHOLESPRICER_H

#include "Option.h"
#include "EuropeanVanillaOption.h"
#include "EuropeanDigitalOption.h"

class BlackScholesPricer {
    private:
        Option* option;
        double asset_price;
        double interest_rate;
        double volatility;
        double N(double a);
    public:
        BlackScholesPricer(EuropeanVanillaOption* opt, double ap, double ir, double vol);
        BlackScholesPricer(EuropeanDigitalOption* opt, double ap, double ir, double vol);
        double operator()();
        double delta();
};

#endif