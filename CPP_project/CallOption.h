#ifndef CALLOPTION_H
#define CALLOPTION_H

#include "EuropeanVanillaOption.h"

class CallOption : public EuropeanVanillaOption{
    public:
        CallOption(double e, double s);
        double payoff(double z) override;
        optionType GetOptionType() const override;
};

#endif