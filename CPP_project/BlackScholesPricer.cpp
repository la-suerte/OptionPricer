#include "BlackScholesPricer.h"
#include <cmath>

//using namespace std;

BlackScholesPricer::BlackScholesPricer(EuropeanVanillaOption* opt, double ap, double ir, double vol)
    : option(opt), asset_price(ap), interest_rate(ir), volatility(vol) {}

BlackScholesPricer::BlackScholesPricer(EuropeanDigitalOption* opt, double ap, double ir, double vol)
    : option(opt), asset_price(ap), interest_rate(ir), volatility(vol) {}

double BlackScholesPricer::N(double a) {
    return std::erfc(-a / std::sqrt(2)) / 2;
}

double BlackScholesPricer::operator()() {
    double T = option->getExpiry();
    double S = asset_price;
    double r = interest_rate;
    double sigma = volatility;

    EuropeanVanillaOption* vanilla = dynamic_cast<EuropeanVanillaOption*>(option);
    if(vanilla) {
        double K = vanilla->_strike;
        double d1 = (log(S/K) + (r + 0.5*sigma*sigma)*T) / (sigma*std::sqrt(T));
        double d2 = d1 - sigma*std::sqrt(T);
        
        if (vanilla->GetOptionType() == call) {
            return S*N(d1) - K*std::exp(-r*T)*N(d2);
        } else {
            return K*std::exp(-r*T)*N(-d2) - S*N(-d1);
        }
    }

    EuropeanDigitalOption* digital = dynamic_cast<EuropeanDigitalOption*>(option);
    if(digital) {
        double K = digital->_strike;
        double d1 = (log(S/K) + (r + 0.5*sigma*sigma)*T) / (sigma*std::sqrt(T));
        double d2 = d1 - sigma*std::sqrt(T);
        
        if(digital->GetOptionType() == call) {
            return std::exp(-r*T)*N(d2);
        } else {
            return std::exp(-r*T)*N(-d2);
        }
    }
    return 0.0;
}

double BlackScholesPricer::delta() {
    double T = option->getExpiry();
    double S = asset_price;
    double r = interest_rate;
    double sigma = volatility;
    
    EuropeanVanillaOption* vanilla = dynamic_cast<EuropeanVanillaOption*>(option);
    if(vanilla) {
        double K = vanilla->_strike;
        double d1 = (log(S/K) + (r + 0.5*sigma*sigma)*T) / (sigma*std::sqrt(T));
        
        if(vanilla->GetOptionType() == call) {
            return N(d1);
        } else {
            return N(d1) - 1;
        }
    }
    
    EuropeanDigitalOption* digital = dynamic_cast<EuropeanDigitalOption*>(option);
    if (digital) {
        double K = digital->_strike;
        double d1 = (log(S/K) + (r + 0.5*sigma*sigma)*T) / (sigma*std::sqrt(T));
        double d2 = d1 - sigma*std::sqrt(T);
        double pdf = (1.0 / std::sqrt(2.0 * M_PI)) * std::exp(-0.5 * d2 * d2);
        double delta = (std::exp(-r*T)*pdf)/(sigma*S*std::sqrt(T));
        
        if(digital->GetOptionType() == call) {
            return delta;
        } else {
            return -delta;
        }
    }
    return 0.0;
}