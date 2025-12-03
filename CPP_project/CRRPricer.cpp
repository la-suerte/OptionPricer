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
    double exponent = (r + (volatility*volatility/2))*h + volatility*std::sqrt(h);
    this->up = std::exp(exponent) - 1;
    exponent = (r + (volatility*volatility/2))*h - volatility*std::sqrt(h);
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
    double spot_down_mult = 1.0 + down;
    double mult_factor = (1.0 + up) / spot_down_mult; 
    double q = (interest_rate - down) / (up - down);

    double St = S0;
    for(int i = 0 ; i < N ; i++)
    {
        St*=spot_down_mult;
    }
    
    for(int i = 0 ; i <= N ; i++)
    {
        _tree->setNode(N, i, option->payoff(St));
        St *= mult_factor; 
    }
    
    St = S0;
    for(int i = 0 ; i < N ; i++)
    {
        St*=spot_down_mult;
    }

    for(int i = N-1; i >= 0; i--) {
        

        St /= spot_down_mult; 
        
        double current_spot_j = St;

        for(int j = 0; j <= i; j++) {
            double value = (q * _tree->getNode(i+1, j+1) + (1 - q) * _tree->getNode(i+1, j)) / (1 + interest_rate);
            
            if(!option->isAmericanOption()) {
                _tree->setNode(i, j, value);
            } else {
                double St_current = option->payoff(current_spot_j);
                _tree->setNode(i, j, max(value, St_current));
                
                if(value <= St_current) {
                    exercised->setNode(i, j, true);
                } else {
                    exercised->setNode(i, j, false);
                }
            }
            current_spot_j *= mult_factor; 
        }
    }
    
    computed = true;
}

double CRRPricer::get(int a, int b) {
    return _tree->getNode(a, b);
}

double CRRPricer::operator()(bool closed_form) {
    if(!computed) {
        compute();
    }
    
    if(closed_form) {
        double sum = 0.0;
        double q = (interest_rate - down) / (up - down);
        double one_minus_q = 1.0 - q;
        
        double discount = 1.0;
        double discount_step = 1.0 / (1.0 + interest_rate);
        for(int i = 0; i < N; i++) {
            discount *= discount_step;
        }


        double prob_weight = 1.0;
        for(int i = 0; i < N; i++) {
            prob_weight *= one_minus_q;
        }
        
        double q_ratio = q / one_minus_q;

        for(int i = 0; i <= N; i++) {
            sum += prob_weight * _tree->getNode(N, i);
            if (i < N) {
                prob_weight *= ((double)(N - i) / (i + 1)) * q_ratio;
            }
        }
        return discount * sum;
    }
    return _tree->getNode(0, 0);
}

void CRRPricer::display_tree() {
    _tree->display();
}