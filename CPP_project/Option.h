#ifndef OPTION_H
#define OPTION_H

#include <vector>

class Option{
    private:
        double _expiry;
    public:
        Option(double e);
        double getExpiry();
        virtual double payoff(double d) = 0;
        virtual double payoffPath(std::vector<double> d);
        virtual bool isAsianOption();
        virtual bool isAmericanOption();
};

#endif