#ifndef EUROPEANDIGITALOPTION_H
#define EUROPEANDIGITALOPTION_H

#include "Option.h"
#include "EuropeanVanillaOption.h"

class EuropeanDigitalOption : public Option
{
    friend class EuropeanDigitalPutOption;
    friend class EuropeanDigitalCallOption;
    friend class BlackScholesPricer;

    private:
        double _strike;
    public:
        EuropeanDigitalOption(double exp, double s);
        virtual optionType GetOptionType() const = 0;
        virtual double payoff(double d) = 0;
};

#endif