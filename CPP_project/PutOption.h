#ifndef PUTOPTION_H
#define PUTOPTION_H

#include "EuropeanVanillaOption.h"

class PutOption : public EuropeanVanillaOption{
    public:
        PutOption(double e, double s);
        double payoff(double z) override;
        optionType GetOptionType() const override;
};

#endif