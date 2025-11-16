#ifndef AMERICANCALLOPTION_H
#define AMERICANCALLOPTION_H

#include "AmericanOption.h"
#include "EuropeanVanillaOption.h"

class AmericanCallOption : public AmericanOption{
    public:
        AmericanCallOption(double e, double s);
        double payoff(double St) override;
        optionType GetOptionType();
};

#endif