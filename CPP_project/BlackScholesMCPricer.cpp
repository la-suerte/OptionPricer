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

void BlackScholesMCPricer::generate(int nb_paths)
{
    double S0 = initial_price;
    double r = interest_rate;
    double T = option->getExpiry();
    double discount = std::exp(-r * T);

    double vanilla_drift = 0.0;
    double vanilla_diffusion = 0.0;
    
    if (!option->isAsianOption()) {
        vanilla_drift = (r - 0.5 * sigma * sigma) * T;
        vanilla_diffusion = sigma * std::sqrt(T);
    }

    std::vector<double> asian_drifts;
    std::vector<double> asian_diffusions;
    AsianOption* asian = nullptr;

    if (option->isAsianOption()) {
        asian = dynamic_cast<AsianOption*>(option);
        std::vector<double> timeSteps = asian->getTimeSteps();
        double prev_t = 0.0;
        
        asian_drifts.resize(timeSteps.size());
        asian_diffusions.resize(timeSteps.size());

        for(size_t t = 0; t < timeSteps.size(); t++) {
            double dt = timeSteps[t] - prev_t;
            asian_drifts[t] = (r - 0.5 * sigma * sigma) * dt;
            asian_diffusions[t] = sigma * std::sqrt(dt);
            prev_t = timeSteps[t];
        }
    }

    for(int i = 0 ; i < nb_paths ; i++)
    {
        double payoff;

        if(option->isAsianOption())
        {
            std::vector<double> St_path;
            double St = S0;

            for(size_t t = 0 ; t < asian_drifts.size() ; t++)
            {
                double Z = MT::rand_norm();
                St *= std::exp(asian_drifts[t] + asian_diffusions[t] * Z);
                St_path.push_back(St);
            }
            payoff = asian->payoffPath(St_path);
        }
        else
        {
            double Z = MT::rand_norm();
            double ST = S0 * std::exp(vanilla_drift + vanilla_diffusion * Z);
            payoff = option->payoff(ST);
        }

        double d_payoff = payoff * discount;
        sumPayoffs += d_payoff;
        sumSquaredPayoffs += d_payoff * d_payoff;
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
