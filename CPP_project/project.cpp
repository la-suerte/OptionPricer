#include "Option.h"
#include "EuropeanVanillaOption.h"
#include "CallOption.h"
#include "PutOption.h"
#include "EuropeanDigitalOption.h"
#include "EuropeanDigitalCallOption.h"
#include "EuropeanDigitalPutOption.h"
#include "AsianOption.h"
#include "AsianCallOption.h"
#include "AsianPutOption.h"
#include "AmericanOption.h"
#include "AmericanCallOption.h"
#include "AmericanPutOption.h"
#include "BlackScholesPricer.h"
#include "CRRPricer.h"
#include "MT.h"
#include "BlackScholesMCPricer.h"
#include <cmath>
#include <stdexcept>
#include <iostream>

using namespace std;


Option::Option(double e) : _expiry(e) {}

double Option::getExpiry() {
    return _expiry;
}

double Option::payoffPath(std::vector<double> d) {
    return payoff(d[d.size()-1]);
}

bool Option::isAsianOption() {
    return false;
}

bool Option::isAmericanOption() {
    return false;
}

EuropeanVanillaOption::EuropeanVanillaOption(double e, double s) 
    : Option(e >= 0 ? e : 0), _strike(s >= 0 ? s : 0) {}

CallOption::CallOption(double e, double s) : EuropeanVanillaOption(e, s) {}

double CallOption::payoff(double z) {
    if(z > _strike)
        return z - _strike;
    else
        return 0;
}

optionType CallOption::GetOptionType() const {
    return call;
}

PutOption::PutOption(double e, double s) : EuropeanVanillaOption(e, s) {}

double PutOption::payoff(double z) {
    if(_strike > z)
        return _strike - z;
    else
        return 0;
}

optionType PutOption::GetOptionType() const {
    return put;
}

EuropeanDigitalOption::EuropeanDigitalOption(double exp, double s) 
    : Option(exp), _strike(s) {}

EuropeanDigitalPutOption::EuropeanDigitalPutOption(double exp, double s) 
    : EuropeanDigitalOption(exp, s) {}

double EuropeanDigitalPutOption::payoff(double z) {
    if(z <= _strike)
        return 1.0;
    else
        return 0;
}

optionType EuropeanDigitalPutOption::GetOptionType() const {
    return put;
}

EuropeanDigitalCallOption::EuropeanDigitalCallOption(double e, double s)
    : EuropeanDigitalOption(e, s) {}

double EuropeanDigitalCallOption::payoff(double z) {
    if(z >= _strike)
        return 1.0;
    else
        return 0;
}

optionType EuropeanDigitalCallOption::GetOptionType() const {
    return call;
}

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

CRRPricer::CRRPricer(Option* option, int depth, double asset_price, double up, double down, double interest_rate)
    : option(option), N(depth), asset_price(asset_price), up(up), down(down), 
      interest_rate(interest_rate), computed(false), closed_form(false), exercised(nullptr) {
    
    if(option->isAsianOption()) {
        throw std::invalid_argument("Le CRRPricer ne peut pas pricer une option asiatique.");
    }
    
    if(!(up > interest_rate && interest_rate > down)) {
        cout << "Arbitrage possible, using default values";
        _tree = new BinaryTree<double>();
        if(option->isAmericanOption()) {
            exercised = new BinaryTree<bool>();
            exercised->setDepth(0);
        }
        _tree->setDepth(0);
        return;
    }
    
    _tree = new BinaryTree<double>();
    if(option->isAmericanOption()) {
        exercised = new BinaryTree<bool>();
        exercised->setDepth(N);
    }
    _tree->setDepth(N);
}

CRRPricer::CRRPricer(Option* option, int depth, double asset_price, double r, double volatility)
    : option(option), N(depth), asset_price(asset_price), computed(false), 
      closed_form(false), exercised(nullptr) {
    
    double h = option->getExpiry()/N;
    double exponent = (r + (std::pow(volatility,2)/2))*h + volatility*std::sqrt(h);
    this->up = std::exp(exponent) - 1;
    exponent = (r + (std::pow(volatility,2)/2))*h - volatility*std::sqrt(h);
    this->down = std::exp(exponent) - 1;
    this->interest_rate = std::exp(r*h) - 1;
    
    _tree = new BinaryTree<double>();
    if(option->isAmericanOption()) {
        exercised = new BinaryTree<bool>();
        exercised->setDepth(N);
    }
    _tree->setDepth(N);
}

CRRPricer::~CRRPricer() {
    if(_tree != nullptr) {
        delete _tree;
    }
    if(exercised != nullptr) {
        delete exercised;
    }
}

bool CRRPricer::getExercise(int a, int b) {
    return exercised->getNode(a, b);
}

double CRRPricer::max(double a, double b) {
    if(a > b)
        return a;
    return b;
}

void CRRPricer::compute() {
    double S0 = asset_price;
    double St;
    
    for(int i = 0; i <= N; i++) {
        St = S0 * pow(1+up, i) * pow(1+down, N-i);
        _tree->setNode(N, i, option->payoff(St));
    }
    
    double q = (interest_rate - down)/(up - down);

    for(int i = N-1; i >= 0; i--) {
        for(int j = 0; j <= i; j++) {
            double value = (q*_tree->getNode(i+1, j+1) + (1-q)*_tree->getNode(i+1, j))/(1+interest_rate);
            
            if(!option->isAmericanOption()) {
                _tree->setNode(i, j, value);
            } else {
                double St_current = option->payoff(S0*std::pow(1+up, j)*std::pow(1+down, i-j));
                _tree->setNode(i, j, max(value, St_current));
                
                if(value <= St_current) {
                    exercised->setNode(i, j, true);
                } else {
                    exercised->setNode(i, j, false);
                }
            }
        }
    }
    computed = true;
}

double CRRPricer::get(int a, int b) {
    return _tree->getNode(a, b);
}

int CRRPricer::fact(int N) {
    int sum = 1;
    for(int i = 1; i <= N; i++) {
        sum = sum * i;
    }
    return sum;
}

double CRRPricer::operator()(bool closed_form) {
    if(!computed) {
        compute();
    }
    
    if(closed_form) {
        double sum = 0.0;
        double q = (interest_rate - down)/(up - down);
        double discount = 1.0/pow(1+interest_rate, N);
        
        for(int i = 0; i <= N; i++) {
            // Calculate binomial coefficient more safely
            double coef = 1.0;
            for(int j = 0; j < i; j++) {
                coef *= (double)(N - j) / (double)(j + 1);
            }
            sum += coef * pow(q, i) * pow(1-q, N-i) * _tree->getNode(N, i);
        }
        
        return discount * sum;
    }
    
    return _tree->getNode(0, 0);
}

void CRRPricer::display_tree() {
    _tree->display();
}

AsianOption::AsianOption(std::vector<double> ts, double e, double k) 
    : Option(e), time_steps(ts), _strike(k) {}

std::vector<double> AsianOption::getTimeSteps() {
    return time_steps;
}

double AsianOption::payoffPath(std::vector<double> St) {
    double sum = 0.0;
    for(int i = 0; i < St.size(); i++) {
        sum += St[i];
    }
    double avg = sum/St.size();
    return payoff(avg);
}

bool AsianOption::isAsianOption() {
    return true;
}

AsianCallOption::AsianCallOption(std::vector<double> ts, double k) 
    : AsianOption(ts, ts.back(), k) {}

double AsianCallOption::payoff(double underlying) {
    if(underlying > _strike)
        return underlying - _strike;
    else
        return 0;
}

AsianPutOption::AsianPutOption(std::vector<double> ts, double k) 
    : AsianOption(ts, ts.back(), k) {}

double AsianPutOption::payoff(double underlying) {
    if(underlying < _strike)
        return _strike - underlying;
    else
        return 0;
}

std::mt19937& MT::get_generator() {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    return gen;
}

double MT::rand_unif() {
    static std::uniform_real_distribution<double> dist(0.0, 1.0);
    return dist(get_generator());
}

double MT::rand_norm() {
    static std::normal_distribution<double> dist(0.0, 1.0);
    return dist(get_generator());
}


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
            vector<double> timeSteps = asian->getTimeSteps();
            vector<double> St_path;

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

AmericanOption::AmericanOption(double e, double s) : Option(e), _strike(s) {}

bool AmericanOption::isAmericanOption() {
    return true;
}

AmericanCallOption::AmericanCallOption(double e, double s) : AmericanOption(e, s) {}

double AmericanCallOption::payoff(double St) {
    if(St > _strike)
        return St - _strike;
    else 
        return 0.0;
}

optionType AmericanCallOption::GetOptionType() {
    return call;
}


AmericanPutOption::AmericanPutOption(double e, double s) : AmericanOption(e, s) {}

double AmericanPutOption::payoff(double St) {
    if(St < _strike)
        return _strike - St;
    else 
        return 0.0;
}

optionType AmericanPutOption::GetOptionType() {
    return put;
}

#include <iostream>
#include <vector>
#include "CallOption.h"
#include "PutOption.h"
#include "EuropeanDigitalCallOption.h"
#include "EuropeanDigitalPutOption.h"
#include "AsianCallOption.h"
#include "AsianPutOption.h"
#include "BlackScholesMCPricer.h"

#include <iostream>
#include <vector>
#include "CallOption.h"
#include "PutOption.h"
#include "EuropeanDigitalCallOption.h"
#include "EuropeanDigitalPutOption.h"
#include "AmericanCallOption.h"
#include "AmericanPutOption.h"
#include "CRRPricer.h"


int main() {
    double S0(95.), K(100.), T(0.5), r(0.02), sigma(0.2);
    std::vector<Option*> opt_ptrs;
    opt_ptrs.push_back(new CallOption(T, K));
    opt_ptrs.push_back(new PutOption(T, K));
    opt_ptrs.push_back(new EuropeanDigitalCallOption(T, K));
    opt_ptrs.push_back(new EuropeanDigitalPutOption(T, K));
    opt_ptrs.push_back(new AmericanCallOption(T, K));
    opt_ptrs.push_back(new AmericanPutOption(T, K));

    CRRPricer* pricer;

    for (auto& opt_ptr : opt_ptrs) {
        pricer = new CRRPricer(opt_ptr, 150, S0, r, sigma);

        pricer->compute();

        std::cout << "price: " << (*pricer)() << std::endl << std::endl;
        delete pricer;
        delete opt_ptr;

    }
}