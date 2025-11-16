#ifndef AMERICANPUTOPTION_H
#define AMERICANPUTOPTION_H

#include "AmericanOption.h"
#include "EuropeanVanillaOption.h"

class AmericanPutOption : public AmericanOption{
    public:
        AmericanPutOption(double e, double s);
        double payoff(double St) override;
        optionType GetOptionType();
};

#endif