#include "CRRPricer.h"
#include <cmath>

CRRPricer::CRRPricer(Option* option, int depth, double asset_price, double up, double down, double interest_rate)
    : option(option), N(depth), asset_price(asset_price), up(up), down(down), 
      interest_rate(interest_rate), computed(false), closed_form(false), exercised(nullptr) {
    
    if(option->isAsianOption()) {
        throw std::invalid_argument("Le CRRPricer ne peut pas pricer une option asiatique.");
    }
    
    if(!(up > interest_rate && interest_rate > down)) {
        std::cout << "Arbitrage possible, using default values";
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