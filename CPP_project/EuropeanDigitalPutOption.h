#ifndef EUROPEANDIGITALPUTOPTION_H
#define EUROPEANDIGITALPUTOPTION_H

#include "EuropeanDigitalOption.h"

class EuropeanDigitalPutOption : public EuropeanDigitalOption{
    public:
        EuropeanDigitalPutOption(double exp, double s);
        double payoff(double z) override;
        optionType GetOptionType() const override;
};

#endif