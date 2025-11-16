#ifndef ASIANCALLOPTION_H
#define ASIANCALLOPTION_H

#include "AsianOption.h"

class AsianCallOption : public AsianOption{
    public:
        AsianCallOption(std::vector<double> ts, double k);
        double payoff(double underlying) override;
};

#endif