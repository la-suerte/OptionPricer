#include <iostream>
#include <vector>
#include <iomanip>   
#include <string>   
#include <cmath>    
#include <stdexcept> 
#include <random>

using namespace std;

//====================================================== PART 1 ======================================================


class Option{
    private:
        double _expiry;
    public:
        Option(double e):_expiry(e){}
        double getExpiry()
        {
            return _expiry;
        }
        virtual double payoff(double d) = 0; 

        virtual double payoffPath( std::vector<double> d){
            return payoff(d[d.size()-1]);
        }

        virtual bool isAsianOption(){return false;};

        virtual bool isAmericanOption(){return false;};
};

enum optionType{
    call,
    put,
};

class EuropeanVanillaOption : public Option{
    //friends:
        friend class CallOption;
        friend class PutOption;
        friend class BlackScholesPricer;

    private:
        double _strike;
    public:
        EuropeanVanillaOption(double e, double s) : Option(e >= 0 ? e : 0), _strike(s >= 0 ? s : 0){}
        virtual optionType GetOptionType() const = 0;

        virtual double payoff(double d) = 0;
};

class CallOption : public EuropeanVanillaOption{
    public:
        CallOption(double e, double s):EuropeanVanillaOption(e,s){}
        double payoff(double z) override
        {
            if(z>_strike)
                return z-_strike;
            else
                return 0;
        }

        optionType GetOptionType() const override{
            return call;
        }
};

class PutOption : public EuropeanVanillaOption{
    public:
        PutOption(double e, double s):EuropeanVanillaOption(e,s){}
        double payoff(double z) override
        {
            if(_strike>z)
                return _strike-z;
            else
                return 0;
        }
        optionType GetOptionType() const override{
            return put;
        }
};


//====================================================== PART 2B ======================================================

class EuropeanDigitalOption : public Option
{
    //friends:
        friend class EuropeanDigitalPutOption;
        friend class EuropeanDigitalCallOption;
        friend class BlackScholesPricer;

    private:
        double _strike;
    public:
        EuropeanDigitalOption(double exp, double s) : Option(exp),_strike(s) {}
        virtual optionType GetOptionType() const = 0;
        virtual double payoff(double d) = 0;
};

class EuropeanDigitalPutOption : public EuropeanDigitalOption{
    public:
        EuropeanDigitalPutOption(double exp, double s) : EuropeanDigitalOption(exp, s){} ;
        double payoff(double z) override
        {
            if(z<=_strike)
                return 1.0;
            else
                return 0;
        }
        optionType GetOptionType() const override
        {
            return put;
        }


};

class EuropeanDigitalCallOption : public EuropeanDigitalOption{
    public:
        EuropeanDigitalCallOption(double e, double s):EuropeanDigitalOption(e,s){}
        double payoff(double z) override
        {
            if(z>=_strike)
                return 1.0;
            else
                return 0;
        }

        optionType GetOptionType() const override
        {
            return call;
        }
};

class BlackScholesPricer {
    private:
        Option* option;
        double asset_price;
        double interest_rate;
        double volatility;
    public:
        BlackScholesPricer(EuropeanVanillaOption* opt, double ap, double ir, double vol):option(opt),asset_price(ap),interest_rate(ir),volatility(vol){}
        double N(double a)
        {
            return std::erfc(-a / std::sqrt(2)) / 2;
        }
        BlackScholesPricer(EuropeanDigitalOption* opt, double ap, double ir, double vol):option(opt),asset_price(ap),interest_rate(ir),volatility(vol){}

        double operator()()
        {
            double T = option->getExpiry();
            double S = asset_price;
            double r = interest_rate;
            double sigma = volatility;

            EuropeanVanillaOption* vanilla = dynamic_cast<EuropeanVanillaOption*>(option);
            if(vanilla)
            {
                double K = vanilla->_strike;
                double d1 = (log(S/K) + (r + 0.5*sigma*sigma)*T) / (sigma*std::sqrt(T));
                double d2 = d1 - sigma*std::sqrt(T);
                
                if (vanilla->GetOptionType() == call) {
                    return S*N(d1) - K*std::exp(-r*T)*N(d2);
                } else {  // put
                    return K*std::exp(-r*T)*N(-d2) - S*N(-d1);
                }
            }

            EuropeanDigitalOption* digital = dynamic_cast<EuropeanDigitalOption*>(option);
            if(digital)
            {
                double K = digital->_strike;
                double d1 = (log(S/K) + (r + 0.5*sigma*sigma)*T) / (sigma*std::sqrt(T));
                double d2 = d1 - sigma*std::sqrt(T);
                
                if(digital->GetOptionType() == call)
                {
                    return std::exp(-r*T)*N(d2);
                }
                else // put
                {
                    return std::exp(-r*T)*N(-d2);
                }
            }
            return 0.0;

        }

        double delta()
        {
            double T = option->getExpiry();
            double S = asset_price;
            double r = interest_rate;
            double sigma = volatility;
            EuropeanVanillaOption* vanilla = dynamic_cast<EuropeanVanillaOption*>(option);
            if(vanilla)
            {
                double K = vanilla->_strike;
                double d1 = (log(S/K) + (r + 0.5*sigma*sigma)*T) / (sigma*std::sqrt(T));
                if(vanilla->GetOptionType() == call)
                {
                    return N(d1);
                }
                else // put
                {
                    return N(d1)-1;
                }
            }
            EuropeanDigitalOption* digital = dynamic_cast<EuropeanDigitalOption*>(option);
            if (digital) {
                double K = digital->_strike;
                double d1 = (log(S/K) + (r + 0.5*sigma*sigma)*T) / (sigma*std::sqrt(T));
                double d2 = d1 - sigma*std::sqrt(T);
                double pdf = (1.0 / std::sqrt(2.0 * M_PI)) * std::exp(-0.5 * d2 * d2);
                double delta = (std::exp(-r*T)*pdf)/(sigma*S*std::sqrt(T));
                if(digital->GetOptionType() == call)
                {
                    return delta;
                }
                else // put
                {
                    return -delta;
                }
            }
            return 0.0;
        }
};

//====================================================== PART 2A ======================================================

template<typename T>

class BinaryTree{
    private:
        int _depth;
        vector<vector<T> > _tree;
    public:
        BinaryTree()
        {
            _depth = 0;
        }
        void setDepth(int d)
        {
            _depth = d;
            _tree.resize(d+1);
            for(int i = 0 ; i<=d ; i++)
            {
                _tree[i].resize(i+1);
            }
        }
        void setNode(int a, int b, T v)
        {
            _tree[a][b]=v;
        }
        T getNode(int a, int b)
        {
            return _tree[a][b];
        }
        void display()
        {
            for(int i = 0 ; i<=_depth; i++)
            {
                for(int j = 0 ; j<=i; j++)
                {
                    if(_tree[i][j]==0)
                        cout<<"0   ";
                    else if(_tree[i][j]<100)
                        cout<<_tree[i][j]<<"  ";
                    else 
                        cout<<_tree[i][j]<<" ";
                }
                cout<<"\n";
            }
            cout<<"\n\n";
            int nodeWidth = 3;  // Each number is 3 digits
            int nodeSpacing = 2; // Spaces between nodes
            
            // Calculate max width needed (bottom row)
            int maxWidth = (_depth + 1) * (nodeWidth + nodeSpacing) - nodeSpacing;
            
            for(int i = 0; i <= _depth; i++) {
                // Calculate leading spaces for this row
                int rowWidth = (i + 1) * (nodeWidth + nodeSpacing) - nodeSpacing;
                int leadingSpaces = (maxWidth - rowWidth) / 2;
                
                // Print leading dashes
                cout << string(leadingSpaces, ' ');
                
                // Print nodes
                for(int j = 0; j <= i; j++) {
                    cout << setfill('0') << setw(3) << _tree[i][j];
                    if(j < i) {
                        cout << "  ";  // 2 dashes between nodes
                    }
                }
                cout << "\n";
                
                if(i < _depth) {
                    cout << string(leadingSpaces - 1, ' ');
                    
                    for(int j = 0; j <= i; j++) {
                        cout << "/";
                        
                        // Special cases for first few rows
                        if(i == 0) {
                            cout << "   ";  //-1
                        } else if(i == 1 && j == 0) {
                            cout << "  ";  //-1
                        } else if(i == 1 && j == 1) {
                            cout << "   ";   
                        } else if(i == 2 && j == 0) {
                            cout << "   ";   
                        } else {
                            cout << "  ";    
                        }
                        
                        cout << "\\";
                        
                        if(j < i) {
                            if(i == 1) {
                                cout << " ";  
                            } else {
                                cout << " ";  
                            }
                        }
                    }
                    cout << "\n";
                }
            }
        }
    };

class CRRPricer{
    private:    
        Option* option;
        int N;
        double asset_price; //S0
        double up;
        double down;
        double interest_rate;
        BinaryTree<double>* _tree;
        bool computed;
        bool closed_form = false;
        BinaryTree<bool>* exercised = nullptr;
    public:
        CRRPricer(Option* option, int depth, double asset_price, double up,double down, double interest_rate)
        {
            if(option->isAsianOption())
            {
                throw std::invalid_argument("Le CRRPricer ne peut pas pricer une option asiatique.");
            }
            if(up>interest_rate && interest_rate>down)
            {
                this->asset_price = asset_price;
                this->up = up;
                this->down = down;
                this->interest_rate = interest_rate;
                this->N = depth;
                this->option = option;
                _tree = new BinaryTree<double>();
                if(option->isAmericanOption())
                {
                    exercised = new BinaryTree<bool>();
                    exercised->setDepth(N);
                }
                _tree->setDepth(N);
                computed = false;
            }
            else{
                cout<<"Arbitrage possible, using default values";
                this->asset_price = asset_price;
                this->up = up;
                this->down = down;
                this->interest_rate = interest_rate;
                this->N = depth;
                this->option = option;
                _tree = new BinaryTree<double>();
                if(option->isAmericanOption())
                {
                    exercised = new BinaryTree<bool>(); 
                    exercised->setDepth(N);
                }
                _tree->setDepth(N);
                computed = false;
            }
        }
        
        CRRPricer(Option* option, int depth, double asset_price, double r, double volatility)
        {
            this->option = option;
            this->N = depth;
            this->asset_price= asset_price;

            double up = 0.0;
            double h = option->getExpiry()/N;
            double exponent = (r+(std::pow(volatility,2)/2))*h + volatility* std::sqrt(h);
            this->up = std::exp(exponent) -1;
            exponent = (r+(std::pow(volatility,2)/2))*h - volatility* std::sqrt(h);
            this->down = std::exp(exponent)-1;
            this->interest_rate = std::exp(r*h)-1;
            _tree = new BinaryTree<double>();
            if(option->isAmericanOption())
            {
                exercised = new BinaryTree<bool>();
                exercised->setDepth(N);
            }
            _tree->setDepth(N);
            computed = false;

        }
        
        bool getExercise(int a, int b)
        {
            return exercised->getNode(a,b);
        }

        double max(double a, double b)
        {
            if(a>b)
                return a;
            return b;
        }
        

        void compute() 
        {
            double S0 = asset_price;
            double St;
            for(int i = 0 ; i <= N ; i++)
            {
                St = S0 * pow(1+up,i)*pow(1+down,N-i);
                _tree->setNode(N,i,option->payoff(St));
            }
            double q = ((interest_rate)-down)/(up-down);

            for(int i = N-1 ; i >= 0 ; i--)
            {
                for(int j = 0 ; j <= i ; j++)
                {
                    double value = (q*_tree->getNode(i+1,j+1)+(1-q)*_tree->getNode(i+1,j))/(1+interest_rate);
                    if(!option->isAmericanOption())
                        _tree->setNode(i,j, value);
                    else{
                        double St_current = option->payoff(S0*std::pow(1+up,j)*std::pow(1+down,i-j));
                        _tree->setNode(i,j,max(value,St_current));
                        if(value<=St_current)
                        {
                            exercised->setNode(i,j,true);
                        }
                        else
                            exercised->setNode(i,j,false);
                    }
                }
            }
            computed = true;
        }
        double get(int a, int b)
        {
            return _tree->getNode(a,b);
        }
        int fact(int N){
            int sum = 1;
            for(int i = 1 ; i<=N ; i++)
            {
                sum=sum*i;
            }
            return sum;
        }
        double operator()(bool closed_form = false)
        {
            if(!computed)
            {
                compute();
            }
            else{
                if(closed_form)
                {
                    double sum = 0.0;
                    double q = ((interest_rate)-down)/(up-down);
                    for(int i = 0 ; i<=N ; i++)
                    {
                        double coef = fact(N)/(fact(i)*fact(N-i));
                        sum+=coef*pow(q,i)*pow(1-q,N-i)*_tree->getNode(N,i);
                    }
                    double fraction = 1/pow(1+interest_rate, N);
                    return fraction*sum ; 
                }
            }
            return _tree->getNode(0,0);
        }
        ~CRRPricer() {
            if(_tree != nullptr) {
                delete _tree;
            }
            delete exercised;
        }
        void display_tree()
        {
            _tree->display();
        }

};

//====================================================== PART 3 ======================================================

class AsianOption : public Option{
    //friends:
        friend class AsianCallOption;
        friend class AsianPutOption;
    private:
        std::vector<double> time_steps;
        double _strike;
    public:
        AsianOption(std::vector<double> ts , double e, double k) : Option(e), time_steps(ts), _strike(k){}

        std::vector<double> getTimeSteps()
        {
            return time_steps;
        }
        double payoffPath(std::vector<double> St) override
        {
            double sum = 0.0;
            for(int i = 0 ; i < St.size() ; i++)
            {
                sum+=St[i];
            }
            double avg = sum/St.size();
            return payoff(avg);
        }

        virtual double payoff(double S) = 0;
        
        bool isAsianOption() override
        {
            return true;
        }

};

class AsianCallOption : public AsianOption{
    public:
        AsianCallOption(vector<double> ts, double k) : AsianOption(ts, ts.back(), k) {}

        double payoff(double underlying)
        {
            if(underlying>_strike)
                return underlying - _strike;
            else
                return 0;
        }
};

class AsianPutOption : public AsianOption{
    public:
        AsianPutOption(vector<double> ts, double k) : AsianOption(ts, ts.back(), k) {}

        double payoff(double underlying)
        {
            if(underlying<_strike)
                return _strike - underlying;
            else
                return 0;
        }
};

class MT{
    private:
        static std::mt19937& get_generator() 
        {
            static std::random_device rd;
            static std::mt19937 gen(rd());
            return gen;
        }
    
    public:
        static double rand_unif()
        {
            static std::uniform_real_distribution<double> dist(0.0, 1.0);
            return dist(get_generator());
        }
        static double rand_norm()
        {
            static std::normal_distribution<double> dist(0.0, 1.0);
            return dist(get_generator());
        }
};

class BlackScholesMCPricer{
    private:
        double estimate;
        int nbPaths = 0; // could be static if "since the beginning of the program" means regardless of the instance
        double sumPayoffs = 0.0;
        double sumSquaredPayoffs = 0.0;
        double initial_price;
        double interest_rate;
        double sigma;
        Option* option;
    public:
        BlackScholesMCPricer(Option* o, double ip, double ir, double vol) : option(o), initial_price(ip), interest_rate(ir), sigma(vol){
        }
        double getNbPaths()
        {
            return nbPaths;
        }
        void generate(int nb_paths)
        {
            double S0 = initial_price;
            double r = interest_rate;
            double T = option->getExpiry();
            double discount = std::exp(-r * T);
            for(int i = 0 ; i < nb_paths ; i++)
            {
                double payoff;

                if(option->isAsianOption())
                {
                    AsianOption* asian = dynamic_cast<AsianOption*>(option);
                    vector<double> timeSteps = asian->getTimeSteps();
                    vector<double> St_path;

                    double prev_t = 0.0;
                    double St = S0;

                    for(int t = 0 ;t < timeSteps.size() ; t++)
                    {
                        double dt = timeSteps[t] - prev_t;
                        double Z = MT::rand_norm();
                        St = St * std::exp((r-sigma*sigma/2.0)*dt + sigma*std::sqrt(dt)*Z);
                        St_path.push_back(St);
                        prev_t = timeSteps[t];
                    }
                    payoff = asian->payoffPath(St_path);

                }
                else
                {
                    double Z = MT::rand_norm();
                    double ST = S0 * std::exp((r-sigma*sigma/2.0)*T + sigma*std::sqrt(T) * Z );
                    payoff = option->payoff(ST);
                }

                double d_payoff = payoff * discount;

                sumPayoffs+=d_payoff;

                sumSquaredPayoffs += d_payoff*d_payoff;

            }

            nbPaths+=nb_paths;

            estimate = sumPayoffs / nbPaths;

        }
        double operator ()()
        {
            if(nbPaths == 0)
                throw std::invalid_argument("ERREUR : Aucun chemin n'a encore été genere ");
            return estimate;
        }

        std::vector<double> confidenceInterval()
        {
            //LLN means price is like a gaussian N(estimate , sumSquaredPayoffs/nbPaths - estimate*estimate )
            if(nbPaths == 0)
                throw std::invalid_argument("ERREUR : Aucun chemin n'a encore été genere ");

            double variance = sumSquaredPayoffs/nbPaths - estimate*estimate;
            double std_error = std::sqrt(variance / nbPaths);

            double margin = 1.96 * std_error ; // because P(-1.96<z<1.96) = 0.95

            std::vector<double> result;
            result.push_back(estimate - margin);
            result.push_back(estimate + margin);
            return result;
        }
};

//====================================================== PART 4 ======================================================

class AmericanOption : public Option {
    //friends:
        friend class AmericanCallOption;
        friend class AmericanPutOption;

    private:
        double _strike;
    public:
        AmericanOption(double e, double s) : _strike(s) , Option(e) {}
        bool isAmericanOption() override
        {
            return true; 
        };
        virtual double payoff(double St) = 0;
};

class AmericanCallOption : public AmericanOption{
    public:
        AmericanCallOption(double e, double s) : AmericanOption(e,s){}

        double payoff(double St) override
        {
            if(St>_strike)
                return St-_strike;
            else 
                return 0.0;
        }
        optionType GetOptionType()
        {
            return call;
        }
};

class AmericanPutOption : public AmericanOption{
    public:
        AmericanPutOption(double e, double s) : AmericanOption(e,s){}

        double payoff(double St) override
        {
            if(St<_strike)
                return _strike -St;
            else 
                return 0.0;
        }
        optionType GetOptionType()
        {
            return put;
        }
};

int main() {
    // Parameters
    double T = 5.0;
    double S0 = 100.0;
    double r = 0.01;
    double vol = 0.1;
    double K = 101.0;
    double R = 0.01;
    double U = 0.05;
    double D = -0.045;
    int N = 5;
    
    cout << "========================================" << endl;
    cout << "OPTION PRICING TEST" << endl;
    cout << "========================================" << endl;
    cout << "Parameters:" << endl;
    cout << "  T = " << T << endl;
    cout << "  S0 = " << S0 << endl;
    cout << "  r = " << r << endl;
    cout << "  vol = " << vol << endl;
    cout << "  K = " << K << endl;
    cout << "  R = " << R << endl;
    cout << "  U = " << U << endl;
    cout << "  D = " << D << endl;
    cout << "  N = " << N << endl;
    cout << "========================================\n" << endl;

    // ============================================================
    // CRR PRICING
    // ============================================================
    cout << "CRR PRICING (Binomial Tree)" << endl;
    cout << "----------------------------------------" << endl;
    
    // European Call
    CallOption eurCall(T, K);
    CRRPricer crrEurCall(&eurCall, N, S0, U, D, R);
    cout << "European Call:         " << crrEurCall() << endl;
    
    // European Digital Call
    EuropeanDigitalCallOption eurDigCall(T, K);
    CRRPricer crrDigCall(&eurDigCall, N, S0, U, D, R);
    cout << "European Digital Call: " << crrDigCall() << endl;
    
    // American Call
    AmericanCallOption amCall(T, K);
    CRRPricer crrAmCall(&amCall, N, S0, U, D, R);
    cout << "American Call:         " << crrAmCall() << endl;
    
    // European Put
    PutOption eurPut(T, K);
    CRRPricer crrEurPut(&eurPut, N, S0, U, D, R);
    cout << "European Put:          " << crrEurPut() << endl;
    
    // European Digital Put
    EuropeanDigitalPutOption eurDigPut(T, K);
    CRRPricer crrDigPut(&eurDigPut, N, S0, U, D, R);
    cout << "European Digital Put:  " << crrDigPut() << endl;
    
    // American Put
    AmericanPutOption amPut(T, K);
    CRRPricer crrAmPut(&amPut, N, S0, U, D, R);
    cout << "American Put:          " << crrAmPut() << endl;
    
    cout << endl;

    // ============================================================
    // BLACK-SCHOLES PRICING (Closed Form)
    // ============================================================
    cout << "BLACK-SCHOLES PRICER (Closed Form)" << endl;
    cout << "----------------------------------------" << endl;
    
    // European Call
    BlackScholesPricer bsEurCall(&eurCall, S0, r, vol);
    cout << "European Call:         " << bsEurCall() << endl;
    cout << "  Delta:               " << bsEurCall.delta() << endl;
    
    // European Digital Call
    BlackScholesPricer bsDigCall(&eurDigCall, S0, r, vol);
    cout << "European Digital Call: " << bsDigCall() << endl;
    cout << "  Delta:               " << bsDigCall.delta() << endl;
    
    // European Put
    BlackScholesPricer bsEurPut(&eurPut, S0, r, vol);
    cout << "European Put:          " << bsEurPut() << endl;
    cout << "  Delta:               " << bsEurPut.delta() << endl;
    
    // European Digital Put
    BlackScholesPricer bsDigPut(&eurDigPut, S0, r, vol);
    cout << "European Digital Put:  " << bsDigPut() << endl;
    cout << "  Delta:               " << bsDigPut.delta() << endl;
    
    cout << "\nWARNING: Black-Scholes cannot price American options (no closed form)" << endl;
    cout << endl;

    // ============================================================
    // MONTE CARLO PRICING
    // ============================================================
    cout << "MONTE CARLO PRICING" << endl;
    cout << "----------------------------------------" << endl;
    
    int num_simulations = 100000;
    cout << "Using " << num_simulations << " simulations" << endl << endl;
    
    // European Call
    BlackScholesMCPricer mcEurCall(&eurCall, S0, r, vol);
    mcEurCall.generate(num_simulations);
    vector<double> ciCall = mcEurCall.confidenceInterval();
    cout << "European Call:         " << mcEurCall() << endl;
    cout << "  95% CI:              [" << ciCall[0] << ", " << ciCall[1] << "]" << endl;
    
    // Asian Call
    vector<double> timeSteps;
    for(int i = 1; i <= 5; i++) {
        timeSteps.push_back(i * T / 5.0);
    }
    AsianCallOption asianCall(timeSteps, K);
    BlackScholesMCPricer mcAsianCall(&asianCall, S0, r, vol);
    mcAsianCall.generate(num_simulations);
    vector<double> ciAsianCall = mcAsianCall.confidenceInterval();
    cout << "Asian Call:            " << mcAsianCall() << endl;
    cout << "  95% CI:              [" << ciAsianCall[0] << ", " << ciAsianCall[1] << "]" << endl;
    
    // European Digital Call
    BlackScholesMCPricer mcDigCall(&eurDigCall, S0, r, vol);
    mcDigCall.generate(num_simulations);
    vector<double> ciDigCall = mcDigCall.confidenceInterval();
    cout << "European Digital Call: " << mcDigCall() << endl;
    cout << "  95% CI:              [" << ciDigCall[0] << ", " << ciDigCall[1] << "]" << endl;
    
    // European Put
    BlackScholesMCPricer mcEurPut(&eurPut, S0, r, vol);
    mcEurPut.generate(num_simulations);
    vector<double> ciPut = mcEurPut.confidenceInterval();
    cout << "European Put:          " << mcEurPut() << endl;
    cout << "  95% CI:              [" << ciPut[0] << ", " << ciPut[1] << "]" << endl;
    
    // European Digital Put
    BlackScholesMCPricer mcDigPut(&eurDigPut, S0, r, vol);
    mcDigPut.generate(num_simulations);
    vector<double> ciDigPut = mcDigPut.confidenceInterval();
    cout << "European Digital Put:  " << mcDigPut() << endl;
    cout << "  95% CI:              [" << ciDigPut[0] << ", " << ciDigPut[1] << "]" << endl;
    
    
    // Asian Put
    AsianPutOption asianPut(timeSteps, K);
    BlackScholesMCPricer mcAsianPut(&asianPut, S0, r, vol);
    mcAsianPut.generate(num_simulations);
    vector<double> ciAsianPut = mcAsianPut.confidenceInterval();
    cout << "Asian Put:             " << mcAsianPut() << endl;
    cout << "  95% CI:              [" << ciAsianPut[0] << ", " << ciAsianPut[1] << "]" << endl;
    
    cout << "\nWARNING: Monte Carlo cannot price American options efficiently" << endl;
    cout << "(requires LSM or similar method, not implemented)" << endl;
    cout << endl;

    // ============================================================
    // COMPARISON TABLE
    // ============================================================
    cout << "========================================" << endl;
    cout << "COMPARISON TABLE" << endl;
    cout << "========================================" << endl;
    cout << left << setw(25) << "Option Type" 
         << setw(15) << "CRR" 
         << setw(15) << "Black-Scholes" 
         << setw(15) << "Monte Carlo" << endl;
    cout << string(70, '-') << endl;
    
    cout << left << setw(25) << "European Call" 
         << setw(15) << crrEurCall() 
         << setw(15) << bsEurCall() 
         << setw(15) << mcEurCall() << endl;
    
    cout << left << setw(25) << "European Digital Call" 
         << setw(15) << crrDigCall() 
         << setw(15) << bsDigCall() 
         << setw(15) << mcDigCall() << endl;
    
    cout << left << setw(25) << "American Call" 
         << setw(15) << crrAmCall() 
         << setw(15) << "N/A" 
         << setw(15) << "N/A" << endl;
    
    cout << left << setw(25) << "European Put" 
         << setw(15) << crrEurPut() 
         << setw(15) << bsEurPut() 
         << setw(15) << mcEurPut() << endl;
    
    cout << left << setw(25) << "European Digital Put" 
         << setw(15) << crrDigPut() 
         << setw(15) << bsDigPut() 
         << setw(15) << mcDigPut() << endl;
    
    cout << left << setw(25) << "American Put" 
         << setw(15) << crrAmPut() 
         << setw(15) << "N/A" 
         << setw(15) << "N/A" << endl;
    
    cout << left << setw(25) << "Asian Call" 
         << setw(15) << "N/A" 
         << setw(15) << "N/A" 
         << setw(15) << mcAsianCall() << endl;
    
    cout << left << setw(25) << "Asian Put" 
         << setw(15) << "N/A" 
         << setw(15) << "N/A" 
         << setw(15) << mcAsianPut() << endl;
    
    cout << endl;
    
    // ============================================================
    // CONVERGENCE TEST (CRR vs Black-Scholes)
    // ============================================================
    cout << "========================================" << endl;
    cout << "CONVERGENCE TEST" << endl;
    cout << "CRR with Black-Scholes parameters (N increasing)" << endl;
    cout << "========================================" << endl;
    
    cout << left << setw(10) << "N" 
         << setw(20) << "European Call CRR" 
         << setw(20) << "BS Price" 
         << setw(15) << "Difference" << endl;
    cout << string(65, '-') << endl;
    
    double bsCallPrice = bsEurCall();
    int nValues[] = {5, 10, 50, 100, 500};
    for(int i = 0; i < 5; i++) {
        int n = nValues[i];
        CallOption testCall(T, K);
        CRRPricer crrTest(&testCall, n, S0, r, vol);
        double crrPrice = crrTest();
        cout << left << setw(10) << n 
             << setw(20) << crrPrice 
             << setw(20) << bsCallPrice 
             << setw(15) << abs(crrPrice - bsCallPrice) << endl;
    }
    
    cout << "\n========================================" << endl;
    cout << "TEST COMPLETE" << endl;
    cout << "========================================" << endl;
    
    return 0;
}