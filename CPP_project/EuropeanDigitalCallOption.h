#ifndef EUROPEANDIGITALCALLOPTION_H
#define EUROPEANDIGITALCALLOPTION_H

#include "EuropeanDigitalOption.h"

class EuropeanDigitalCallOption : public EuropeanDigitalOption{
    public:
        EuropeanDigitalCallOption(double e, double s);
        double payoff(double z) override;
        optionType GetOptionType() const override;
};

#endif