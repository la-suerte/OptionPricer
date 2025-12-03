#ifndef BLACKSCHOLESMCPRICER_H
#define BLACKSCHOLESMCPRICER_H

#include "Option.h"
#include <vector>

class BlackScholesMCPricer{
    private:
        double estimate;
        int nbPaths;
        double sumPayoffs;
        double sumSquaredPayoffs;
        double initial_price;
        double interest_rate;
        double sigma;
        Option* option;
    public:
        BlackScholesMCPricer(Option* o, double ip, double ir, double vol);
        int getNbPaths();
        void generate(int nb_paths);
        double operator()();
        std::vector<double> confidenceInterval();
};

#endif