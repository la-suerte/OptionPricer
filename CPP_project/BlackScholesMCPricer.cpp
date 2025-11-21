#include "BlackScholesMCPricer.h"
#include <cmath>
#include <vector>
#include "AsianOption.h"
#include "MT.h"
#include <stdexcept>


BlackScholesMCPricer::BlackScholesMCPricer(Option* o, double ip, double ir, double vol) 
    : option(o), initial_price(ip), interest_rate(ir), sigma(vol), 
      nbPaths(0), sumPayoffs(0.0), sumSquaredPayoffs(0.0) {}

double BlackScholesMCPricer::getNbPaths() {
    return nbPaths;
}

void BlackScholesMCPricer::generate(int nb_paths) {
    double S0 = initial_price;
    double r = interest_rate;
    double T = option->getExpiry();
    double discount = std::exp(-r * T);
    
    for(int i = 0; i < nb_paths; i++) {
        double payoff;

        if(option->isAsianOption()) {
            AsianOption* asian = dynamic_cast<AsianOption*>(option);
            std::vector<double> timeSteps = asian->getTimeSteps();
            std::vector<double> St_path;

            double prev_t = 0.0;
            double St = S0;

            for(int t = 0; t < timeSteps.size(); t++) {
                double dt = timeSteps[t] - prev_t;
                double Z = MT::rand_norm();
                St = St * std::exp((r - sigma*sigma/2.0)*dt + sigma*std::sqrt(dt)*Z);
                St_path.push_back(St);
                prev_t = timeSteps[t];
            }
            payoff = asian->payoffPath(St_path);
        } else {
            double Z = MT::rand_norm();
            double ST = S0 * std::exp((r - sigma*sigma/2.0)*T + sigma*std::sqrt(T)*Z);
            payoff = option->payoff(ST);
        }

        double d_payoff = payoff * discount;
        sumPayoffs += d_payoff;
        sumSquaredPayoffs += d_payoff*d_payoff;
    }

    nbPaths += nb_paths;
    estimate = sumPayoffs / nbPaths;
}

double BlackScholesMCPricer::operator()() {
    if(nbPaths == 0)
        throw std::invalid_argument("ERREUR : Aucun chemin n'a encore été genere ");
    return estimate;
}

std::vector<double> BlackScholesMCPricer::confidenceInterval() {
    if(nbPaths == 0)
        throw std::invalid_argument("ERREUR : Aucun chemin n'a encore été genere ");

    double variance = sumSquaredPayoffs/nbPaths - estimate*estimate;
    double std_error = std::sqrt(variance / nbPaths);
    double margin = 1.96 * std_error;

    std::vector<double> result;
    result.push_back(estimate - margin);
    result.push_back(estimate + margin);
    return result;
}
