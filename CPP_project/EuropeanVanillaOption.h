#ifndef EUROPEANVANILLAOPTION_H
#define EUROPEANVANILLAOPTION_H

#include "Option.h"

enum optionType{
    call,
    put,
};

class EuropeanVanillaOption : public Option{
    friend class CallOption;
    friend class PutOption;
    friend class BlackScholesPricer;

    private:
        double _strike;
    public:
        EuropeanVanillaOption(double e, double s);
        virtual optionType GetOptionType() const = 0;
        virtual double payoff(double d) = 0;
};

#endif