#ifndef ASIANOPTION_H
#define ASIANOPTION_H

#include "Option.h"
#include <vector>

class AsianOption : public Option{
    friend class AsianCallOption;
    friend class AsianPutOption;
    
private:
    std::vector<double> time_steps;
    double _strike;
    
public:
    AsianOption(std::vector<double> ts, double e, double k);
    std::vector<double> getTimeSteps();
    double payoffPath(std::vector<double> St) override;
    virtual double payoff(double S) = 0;
    bool isAsianOption() override;
};

#endif