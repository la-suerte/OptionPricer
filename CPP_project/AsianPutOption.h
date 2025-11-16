#ifndef ASIANPUTOPTION_H
#define ASIANPUTOPTION_H

#include "AsianOption.h"

class AsianPutOption : public AsianOption{
    public:
        AsianPutOption(std::vector<double> ts, double k);
        double payoff(double underlying) override;
};

#endif