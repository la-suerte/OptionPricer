#ifndef AMERICANOPTION_H
#define AMERICANOPTION_H

#include "Option.h"

class AmericanOption : public Option {
    friend class AmericanCallOption;
    friend class AmericanPutOption;

    private:
        double _strike;
    public:
        AmericanOption(double e, double s);
        bool isAmericanOption() override;
        virtual double payoff(double St) = 0;
};

#endif