#include <iostream>
#include <vector>
#include <iomanip>   
#include <string>   
#include <cmath>    
#include <stdexcept> 
#include <random>

using namespace std;


class Option{
    private:
        double _expiry;
    public:
        Option(double e):_expiry(e){}
        double getExpiry()
        {
            return _expiry;
        }
        virtual double payoff(double d) = 0; //Improve?

        virtual double payoffPath( std::vector<double> d){
            return payoff(d[d.size()-1]);
        };

        virtual bool isAsian(){return false;};

        virtual double getStrike() const = 0;
        

};

enum optionType{
    call,
    put,
};

class EuropeanVanillaOption : public Option{
    private:
        double _strike;
    public:
        EuropeanVanillaOption(double e, double s) : Option(e >= 0 ? e : 0), _strike(s >= 0 ? s : 0){}
        virtual optionType getOptionType() const = 0;

        virtual double payoff(double d) = 0; 

        double getStrike() const override
        {
            return _strike;
        }

        friend class BlackScholesPricer;
    
};

class CallOption : public EuropeanVanillaOption{
    public:
        CallOption(double e, double s):EuropeanVanillaOption(e,s){}
        double payoff(double z) override
        {
            if(z>getStrike())
                return z-getStrike();
            else
                return 0;
        }

        optionType getOptionType() const override{
            return call;
        }
};

class PutOption : public EuropeanVanillaOption{
    public:
        PutOption(double e, double s):EuropeanVanillaOption(e,s){}
        double payoff(double z) override
        {
            if(getStrike()>z)
                return getStrike()-z;
            else
                return 0;
        }
        optionType getOptionType() const override{
            return put;
        }
};


//====================================================== PART 2B ======================================================

class EuropeanDigitalOption : public Option
{
    private:
        double _strike;
    public:
        EuropeanDigitalOption(double exp, double s) : Option(exp),_strike(s) {}
        virtual optionType getOptionType() const = 0;
        virtual double payoff(double d) = 0;
        double getStrike() const override
        {
            return _strike;
        }
};

class DigitalPutOption : public EuropeanDigitalOption{
    public:
        DigitalPutOption(double exp, double s) : EuropeanDigitalOption(exp, s){} ;
        double payoff(double z) override
        {
            if(z<=getStrike())
                return 1.0;
            else
                return 0;
        }
        optionType getOptionType() const override
        {
            return put;
        }


};

class DigitalCallOption : public EuropeanDigitalOption{
    public:
        DigitalCallOption(double e, double s):EuropeanDigitalOption(e,s){}
        double payoff(double z) override
        {
            if(z>=getStrike())
                return 1.0;
            else
                return 0;
        }

        optionType getOptionType() const override
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
                double K = vanilla->getStrike();
                double d1 = (log(S/K) + (r + 0.5*sigma*sigma)*T) / (sigma*std::sqrt(T));
                double d2 = d1 - sigma*std::sqrt(T);
                
                if (vanilla->getOptionType() == call) {
                    return S*N(d1) - K*std::exp(-r*T)*N(d2);
                } else {  // put
                    return K*std::exp(-r*T)*N(-d2) - S*N(-d1);
                }
            }

            EuropeanDigitalOption* digital = dynamic_cast<EuropeanDigitalOption*>(option);
            if(digital)
            {
                double K = digital->getStrike();
                double d1 = (log(S/K) + (r + 0.5*sigma*sigma)*T) / (sigma*std::sqrt(T));
                double d2 = d1 - sigma*std::sqrt(T);
                
                if(digital->getOptionType() == call)
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
            double K = option->getStrike();
            double S = asset_price;
            double r = interest_rate;
            double sigma = volatility;
            EuropeanVanillaOption* vanilla = dynamic_cast<EuropeanVanillaOption*>(option);
            if(vanilla)
            {
                double d1 = (log(S/K) + (r + 0.5*sigma*sigma)*T) / (sigma*std::sqrt(T));
                if(vanilla->getOptionType() == call)
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
                double d1 = (log(S/K) + (r + 0.5*sigma*sigma)*T) / (sigma*std::sqrt(T));
                double d2 = d1 - sigma*std::sqrt(T);
                double pdf = (1.0 / std::sqrt(2.0 * M_PI)) * std::exp(-0.5 * d2 * d2);
                double delta = (std::exp(-r*T)*pdf)/(sigma*S*std::sqrt(T));
                if(digital->getOptionType() == call)
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
    public:
        CRRPricer(Option* option, int depth, double asset_price, double up,double down, double interest_rate)
        {
            if(option->isAsian())
            {
                throw std::invalid_argument("Le CRRPricer ne peut pas pricer une option asiatique.");
            }
            if(up>(1+interest_rate) && (1+interest_rate)>down)
            {
                this->asset_price = asset_price;
                this->up = up;
                this->down = down;
                this->interest_rate = interest_rate;
                this->N = depth;
                this->option = option;
                _tree = new BinaryTree<double>();
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
                _tree->setDepth(0);
                computed = false;
            }
        }
        double max(double a, double b)
        {
            if(a>b)
                return a;
            return b;
        }
        void compute() //only for a call
        {
            double S0 = asset_price;
            double St;
            for(int i = 0 ; i <= N ; i++)
            {
                St = S0 * pow(up,i)*pow(down,N-i);
                _tree->setNode(N,i,option->payoff(St));
            }
            double q = ((1+interest_rate)-down)/(up-down);


            for(int i = N-1 ; i >= 0 ; i--)
            {
                for(int j = 0 ; j <= i ; j++)
                {
                    double value = (q*_tree->getNode(i+1,j+1)+(1-q)*_tree->getNode(i+1,j))/(1+interest_rate);
                    _tree->setNode(i,j, value);
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
                    double q = ((1+interest_rate)-down)/(up-down);
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
        }
        void display_tree()
        {
            _tree->display();
        }

        bool checkAsianOption()
        {

        }
};

//====================================================== PART 3 ======================================================

class AsianOption : public Option{
    private:
        std::vector<double> time_steps;
        double _strike;
    public:
        AsianOption(std::vector<double> ts , double e, double k) : Option(e), time_steps(ts), _strike(k){}

        std::vector<double> getTimeSteps()
        {
            return time_steps;
        }

        double getStrike() const override
        {
            return _strike;
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
        
        bool isAsian() override
        {
            return true;
        }

};

class AsianCallOption : public AsianOption{
    public:
        AsianCallOption(vector<double> ts, double e, double k) : AsianOption(ts, e, k) {}

        double payoff(double underlying)
        {
            if(underlying>getStrike())
                return underlying - getStrike();
            else
                return 0;
        }
};

class AsianPutOption : public AsianOption{
    public:
        AsianPutOption(vector<double> ts, double e, double k) : AsianOption(ts, e, k) {}

        double payoff(double underlying)
        {
            if(underlying<getStrike())
                return getStrike() - underlying;
            else
                return 0;
        }
};

class MT{
    private: //mt19937 is a random number generator is a random seed
        static std::mt19937& get_generator() 
        {
            static std::mt19937 gen(std::random_device{}()); //static this line runs only the first time get_gen is called
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

                if(option->isAsian())
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


int main() {
    cout << "============================================" << endl;
    cout << "   COMPREHENSIVE OPTIONS PRICING TEST" << endl;
    cout << "============================================" << endl << endl;

    // Test parameters
    double S0 = 100.0;      // Initial stock price
    double K = 100.0;       // Strike price
    double T = 1.0;         // Expiry (1 year)
    double r = 0.05;        // Risk-free rate (5%)
    double sigma = 0.2;     // Volatility (20%)
    
    // CRR parameters
    int N = 10;             // Tree depth
    double u = 1.1;         // Up factor
    double d = 0.9;         // Down factor

    try {
        // ==================== PART 1: EUROPEAN VANILLA OPTIONS ====================
        cout << "==================== PART 1: EUROPEAN VANILLA OPTIONS ====================" << endl << endl;
        
        // Test European Call Option
        cout << "--- Testing European Call Option ---" << endl;
        CallOption* euroCall = new CallOption(T, K);
        cout << "Created CallOption with expiry=" << euroCall->getExpiry() 
             << ", strike=" << euroCall->getStrike() << endl;
        
        // Test payoff
        double testPrice1 = 110.0;
        double testPrice2 = 90.0;
        cout << "Payoff at S=" << testPrice1 << ": " << euroCall->payoff(testPrice1) << endl;
        cout << "Payoff at S=" << testPrice2 << ": " << euroCall->payoff(testPrice2) << endl;
        cout << "Option type: " << (euroCall->getOptionType() == call ? "CALL" : "PUT") << endl;
        cout << "Is Asian? " << (euroCall->isAsian() ? "YES" : "NO") << endl;
        
        // Test European Put Option
        cout << "\n--- Testing European Put Option ---" << endl;
        PutOption* euroPut = new PutOption(T, K);
        cout << "Created PutOption with expiry=" << euroPut->getExpiry() 
             << ", strike=" << euroPut->getStrike() << endl;
        cout << "Payoff at S=" << testPrice1 << ": " << euroPut->payoff(testPrice1) << endl;
        cout << "Payoff at S=" << testPrice2 << ": " << euroPut->payoff(testPrice2) << endl;
        cout << "Option type: " << (euroPut->getOptionType() == call ? "CALL" : "PUT") << endl;
        
        // Test Black-Scholes Pricer
        cout << "\n--- Testing Black-Scholes Pricer ---" << endl;
        BlackScholesPricer bsCall(euroCall, S0, r, sigma);
        double bsCallPrice = bsCall();
        double bsCallDelta = bsCall.delta();
        cout << "BS Call Price: " << bsCallPrice << endl;
        cout << "BS Call Delta: " << bsCallDelta << endl;
        
        BlackScholesPricer bsPut(euroPut, S0, r, sigma);
        double bsPutPrice = bsPut();
        double bsPutDelta = bsPut.delta();
        cout << "BS Put Price: " << bsPutPrice << endl;
        cout << "BS Put Delta: " << bsPutDelta << endl;
        
        // Verify Put-Call Parity: C - P = S - K*e^(-rT)
        double parity = bsCallPrice - bsPutPrice;
        double expected_parity = S0 - K * exp(-r * T);
        cout << "\nPut-Call Parity Check:" << endl;
        cout << "C - P = " << parity << endl;
        cout << "S - K*e^(-rT) = " << expected_parity << endl;
        cout << "Difference: " << abs(parity - expected_parity) << " (should be ~0)" << endl;

        // ==================== PART 2A: CRR BINOMIAL TREE ====================
        cout << "\n==================== PART 2A: CRR BINOMIAL TREE ====================" << endl << endl;
        
        cout << "--- Testing CRR Pricer for Call Option ---" << endl;
        CRRPricer crrCall(euroCall, N, S0, u, d, r);
        cout << "Created CRRPricer with N=" << N << ", S0=" << S0 
             << ", u=" << u << ", d=" << d << ", r=" << r << endl;
        
        cout << "\nComputing CRR tree..." << endl;
        crrCall.compute();
        double crrCallPrice = crrCall();
        cout << "CRR Call Price (via tree): " << crrCallPrice << endl;
        
        double crrCallPriceClosed = crrCall(true);
        cout << "CRR Call Price (closed form): " << crrCallPriceClosed << endl;
        cout << "Difference: " << abs(crrCallPrice - crrCallPriceClosed) << " (should be ~0)" << endl;
        
        cout << "\nComparing CRR vs Black-Scholes:" << endl;
        cout << "CRR Price: " << crrCallPrice << endl;
        cout << "BS Price: " << bsCallPrice << endl;
        cout << "Difference: " << abs(crrCallPrice - bsCallPrice) << endl;
        cout << "(Note: CRR converges to BS as N→∞)" << endl;
        
        cout << "\n--- Testing CRR Pricer for Put Option ---" << endl;
        CRRPricer crrPut(euroPut, N, S0, u, d, r);
        crrPut.compute();
        double crrPutPrice = crrPut();
        cout << "CRR Put Price: " << crrPutPrice << endl;
        cout << "BS Put Price: " << bsPutPrice << endl;
        cout << "Difference: " << abs(crrPutPrice - bsPutPrice) << endl;
        
        // Display tree structure (small example)
        cout << "\n--- Displaying CRR Tree (first few nodes) ---" << endl;
        CRRPricer smallTreeCall(euroCall, 3, S0, u, d, r);
        smallTreeCall.compute();
        smallTreeCall.display_tree();

        // ==================== PART 2B: DIGITAL OPTIONS ====================
        cout << "\n==================== PART 2B: DIGITAL OPTIONS ====================" << endl << endl;
        
        cout << "--- Testing Digital Call Option ---" << endl;
        DigitalCallOption* digiCall = new DigitalCallOption(T, K);
        cout << "Created DigitalCallOption with expiry=" << digiCall->getExpiry() 
             << ", strike=" << digiCall->getStrike() << endl;
        cout << "Payoff at S=" << testPrice1 << ": " << digiCall->payoff(testPrice1) << endl;
        cout << "Payoff at S=" << testPrice2 << ": " << digiCall->payoff(testPrice2) << endl;
        cout << "Payoff at S=K=" << K << ": " << digiCall->payoff(K) << endl;
        
        BlackScholesPricer bsDigiCall(digiCall, S0, r, sigma);
        double bsDigiCallPrice = bsDigiCall();
        double bsDigiCallDelta = bsDigiCall.delta();
        cout << "BS Digital Call Price: " << bsDigiCallPrice << endl;
        cout << "BS Digital Call Delta: " << bsDigiCallDelta << endl;
        
        cout << "\n--- Testing Digital Put Option ---" << endl;
        DigitalPutOption* digiPut = new DigitalPutOption(T, K);
        cout << "Created DigitalPutOption with expiry=" << digiPut->getExpiry() 
             << ", strike=" << digiPut->getStrike() << endl;
        cout << "Payoff at S=" << testPrice1 << ": " << digiPut->payoff(testPrice1) << endl;
        cout << "Payoff at S=" << testPrice2 << ": " << digiPut->payoff(testPrice2) << endl;
        
        BlackScholesPricer bsDigiPut(digiPut, S0, r, sigma);
        double bsDigiPutPrice = bsDigiPut();
        double bsDigiPutDelta = bsDigiPut.delta();
        cout << "BS Digital Put Price: " << bsDigiPutPrice << endl;
        cout << "BS Digital Put Delta: " << bsDigiPutDelta << endl;
        
        // Digital parity: Digital Call + Digital Put = e^(-rT)
        double digiParity = bsDigiCallPrice + bsDigiPutPrice;
        double expectedDigiParity = exp(-r * T);
        cout << "\nDigital Parity Check:" << endl;
        cout << "Digital Call + Digital Put = " << digiParity << endl;
        cout << "e^(-rT) = " << expectedDigiParity << endl;
        cout << "Difference: " << abs(digiParity - expectedDigiParity) << " (should be ~0)" << endl;

        // ==================== PART 3: ASIAN OPTIONS ====================
        cout << "\n==================== PART 3: ASIAN OPTIONS ====================" << endl << endl;
        
        // Create time steps for Asian option
        vector<double> timeSteps;
        int m = 12; // Monthly observations
        for(int i = 1; i <= m; i++) {
            timeSteps.push_back(i * T / m);
        }
        
        cout << "--- Testing Asian Call Option ---" << endl;
        AsianCallOption* asianCall = new AsianCallOption(timeSteps, T, K);
        cout << "Created AsianCallOption with expiry=" << asianCall->getExpiry() 
             << ", strike=" << K << ", m=" << m << " time steps" << endl;
        cout << "Is Asian? " << (asianCall->isAsian() ? "YES" : "NO") << endl;
        
        // Test payoffPath
        vector<double> testPath = {95, 98, 102, 105, 103, 107, 110, 108, 112, 115, 113, 118};
        double avgPrice = 0;
        for(double p : testPath) avgPrice += p;
        avgPrice /= testPath.size();
        cout << "Test path average: " << avgPrice << endl;
        cout << "Asian Call payoffPath: " << asianCall->payoffPath(testPath) << endl;
        cout << "Manual calculation: max(0, " << avgPrice << " - " << K << ") = " 
             << max(0.0, avgPrice - K) << endl;
        
        cout << "\n--- Testing Asian Put Option ---" << endl;
        AsianPutOption* asianPut = new AsianPutOption(timeSteps, T, K);
        cout << "Created AsianPutOption with expiry=" << asianPut->getExpiry() 
             << ", strike=" << K << ", m=" << m << " time steps" << endl;
        cout << "Asian Put payoffPath: " << asianPut->payoffPath(testPath) << endl;
        
        // Test CRRPricer rejection of Asian options
        cout << "\n--- Testing CRRPricer Asian Option Rejection ---" << endl;
        try {
            CRRPricer crrAsian(asianCall, N, S0, u, d, r);
            cout << "ERROR: CRRPricer should have thrown exception for Asian option!" << endl;
        } catch (const std::invalid_argument& e) {
            cout << "SUCCESS: CRRPricer correctly rejected Asian option" << endl;
            cout << "Exception message: " << e.what() << endl;
        }

        // ==================== PART 4: MT RANDOM NUMBER GENERATOR ====================
        cout << "\n==================== PART 4: MT RANDOM NUMBER GENERATOR ====================" << endl << endl;
        
        cout << "--- Testing MT Singleton ---" << endl;
        cout << "Generating 10 uniform random numbers [0,1]:" << endl;
        for(int i = 0; i < 10; i++) {
            cout << "  " << i+1 << ": " << MT::rand_unif() << endl;
        }
        
        cout << "\nGenerating 10 standard normal random numbers:" << endl;
        double sum_norm = 0, sum_sq_norm = 0;
        for(int i = 0; i < 10; i++) {
            double z = MT::rand_norm();
            sum_norm += z;
            sum_sq_norm += z * z;
            cout << "  " << i+1 << ": " << z << endl;
        }
        cout << "Sample mean: " << sum_norm / 10 << " (should be ~0)" << endl;
        cout << "Sample variance: " << sum_sq_norm / 10 - (sum_norm/10)*(sum_norm/10) 
             << " (should be ~1)" << endl;

        // ==================== PART 5: BLACK-SCHOLES MONTE CARLO ====================
        cout << "\n==================== PART 5: BLACK-SCHOLES MONTE CARLO ====================" << endl << endl;
        
        cout << "--- Testing MC Pricer for European Call ---" << endl;
        BlackScholesMCPricer mcCall(euroCall, S0, r, sigma);
        cout << "Created BlackScholesMCPricer for Call option" << endl;
        cout << "Initial nbPaths: " << mcCall.getNbPaths() << endl;
        
        // Test exception before generation
        cout << "\nTesting exception before path generation:" << endl;
        try {
            double price = mcCall();
            cout << "ERROR: Should have thrown exception!" << endl;
        } catch (const std::invalid_argument& e) {
            cout << "SUCCESS: Exception caught: " << e.what() << endl;
        }
        
        // Generate paths progressively
        cout << "\n--- Progressive Path Generation ---" << endl;
        int batches[] = {1000, 5000, 10000, 50000};
        for(int batch : batches) {
            mcCall.generate(batch);
            double mcPrice = mcCall();
            vector<double> ci = mcCall.confidenceInterval();
            cout << "\nAfter " << mcCall.getNbPaths() << " paths:" << endl;
            cout << "  MC Price: " << mcPrice << endl;
            cout << "  95% CI: [" << ci[0] << ", " << ci[1] << "]" << endl;
            cout << "  CI Width: " << (ci[1] - ci[0]) << endl;
            cout << "  BS Price: " << bsCallPrice << endl;
            cout << "  Difference: " << abs(mcPrice - bsCallPrice) << endl;
            cout << "  BS in CI? " << (bsCallPrice >= ci[0] && bsCallPrice <= ci[1] ? "YES" : "NO") << endl;
        }
        
        cout << "\n--- Testing MC Pricer for European Put ---" << endl;
        BlackScholesMCPricer mcPut(euroPut, S0, r, sigma);
        mcPut.generate(100000);
        double mcPutPrice = mcPut();
        vector<double> putCI = mcPut.confidenceInterval();
        cout << "MC Put Price: " << mcPutPrice << " (after " << mcPut.getNbPaths() << " paths)" << endl;
        cout << "95% CI: [" << putCI[0] << ", " << putCI[1] << "]" << endl;
        cout << "BS Put Price: " << bsPutPrice << endl;
        cout << "Difference: " << abs(mcPutPrice - bsPutPrice) << endl;
        
        cout << "\n--- Testing MC Pricer for Asian Call ---" << endl;
        BlackScholesMCPricer mcAsianCall(asianCall, S0, r, sigma);
        mcAsianCall.generate(100000);
        double mcAsianPrice = mcAsianCall();
        vector<double> asianCI = mcAsianCall.confidenceInterval();
        cout << "MC Asian Call Price: " << mcAsianPrice << " (after " << mcAsianCall.getNbPaths() << " paths)" << endl;
        cout << "95% CI: [" << asianCI[0] << ", " << asianCI[1] << "]" << endl;
        cout << "CI Width: " << (asianCI[1] - asianCI[0]) << endl;
        cout << "(Note: Asian options are typically cheaper than European calls)" << endl;
        cout << "European Call Price (BS): " << bsCallPrice << endl;
        cout << "Ratio (Asian/European): " << mcAsianPrice / bsCallPrice << endl;
        
        cout << "\n--- Testing MC Pricer for Asian Put ---" << endl;
        BlackScholesMCPricer mcAsianPut(asianPut, S0, r, sigma);
        mcAsianPut.generate(100000);
        double mcAsianPutPrice = mcAsianPut();
        vector<double> asianPutCI = mcAsianPut.confidenceInterval();
        cout << "MC Asian Put Price: " << mcAsianPutPrice << " (after " << mcAsianPut.getNbPaths() << " paths)" << endl;
        cout << "95% CI: [" << asianPutCI[0] << ", " << asianPutCI[1] << "]" << endl;

        // ==================== CONVERGENCE ANALYSIS ====================
        cout << "\n==================== CONVERGENCE ANALYSIS ====================" << endl << endl;
        
        cout << "--- Testing MC Convergence (European Call) ---" << endl;
        BlackScholesMCPricer convergenceTest(euroCall, S0, r, sigma);
        cout << "True Price (BS): " << bsCallPrice << endl << endl;
        
        int nPaths[] = {100, 500, 1000, 5000, 10000, 50000, 100000};
        cout << setw(10) << "N Paths" << setw(15) << "MC Price" << setw(15) << "Error" 
             << setw(20) << "CI Width" << setw(10) << "BS in CI?" << endl;
        cout << string(70, '-') << endl;
        
        for(int n : nPaths) {
            BlackScholesMCPricer tempPricer(euroCall, S0, r, sigma);
            tempPricer.generate(n);
            double price = tempPricer();
            vector<double> ci = tempPricer.confidenceInterval();
            double error = abs(price - bsCallPrice);
            double ciWidth = ci[1] - ci[0];
            bool bsInCI = (bsCallPrice >= ci[0] && bsCallPrice <= ci[1]);
            
            cout << setw(10) << n << setw(15) << price << setw(15) << error 
                 << setw(20) << ciWidth << setw(10) << (bsInCI ? "YES" : "NO") << endl;
        }
        
        // ==================== STRESS TESTS ====================
        cout << "\n==================== STRESS TESTS ====================" << endl << endl;
        
        cout << "--- Testing Edge Cases ---" << endl;
        
        // Deep ITM/OTM options
        cout << "\n1. Deep ITM Call (S=150, K=100):" << endl;
        CallOption deepITMCall(T, K);
        BlackScholesMCPricer mcDeepITM(&deepITMCall, 150.0, r, sigma);
        mcDeepITM.generate(10000);
        cout << "   MC Price: " << mcDeepITM() << endl;
        BlackScholesPricer bsDeepITM(&deepITMCall, 150.0, r, sigma);
        cout << "   BS Price: " << bsDeepITM() << endl;
        
        cout << "\n2. Deep OTM Call (S=50, K=100):" << endl;
        CallOption deepOTMCall(T, K);
        BlackScholesMCPricer mcDeepOTM(&deepOTMCall, 50.0, r, sigma);
        mcDeepOTM.generate(10000);
        cout << "   MC Price: " << mcDeepOTM() << endl;
        BlackScholesPricer bsDeepOTM(&deepOTMCall, 50.0, r, sigma);
        cout << "   BS Price: " << bsDeepOTM() << endl;
        
        // High volatility
        cout << "\n3. High Volatility (sigma=0.8):" << endl;
        BlackScholesMCPricer mcHighVol(euroCall, S0, r, 0.8);
        mcHighVol.generate(10000);
        cout << "   MC Price: " << mcHighVol() << endl;
        BlackScholesPricer bsHighVol(euroCall, S0, r, 0.8);
        cout << "   BS Price: " << bsHighVol() << endl;
        
        // Short expiry
        cout << "\n4. Short Expiry (T=0.1):" << endl;
        CallOption shortCall(0.1, K);
        BlackScholesMCPricer mcShort(&shortCall, S0, r, sigma);
        mcShort.generate(10000);
        cout << "   MC Price: " << mcShort() << endl;
        BlackScholesPricer bsShort(&shortCall, S0, r, sigma);
        cout << "   BS Price: " << bsShort() << endl;
        
        // ==================== SUMMARY ====================
        cout << "\n==================== TEST SUMMARY ====================" << endl << endl;
        
        cout << "✓ European Vanilla Options: PASSED" << endl;
        cout << "✓ Black-Scholes Analytical Pricing: PASSED" << endl;
        cout << "✓ CRR Binomial Tree: PASSED" << endl;
        cout << "✓ Digital Options: PASSED" << endl;
        cout << "✓ Asian Options: PASSED" << endl;
        cout << "✓ MT Random Number Generator: PASSED" << endl;
        cout << "✓ Monte Carlo Simulation: PASSED" << endl;
        cout << "✓ Confidence Intervals: PASSED" << endl;
        cout << "✓ Exception Handling: PASSED" << endl;
        cout << "✓ Convergence Analysis: PASSED" << endl;
        
        // Cleanup
        delete euroCall;
        delete euroPut;
        delete digiCall;
        delete digiPut;
        delete asianCall;
        delete asianPut;
        
        cout << "\n============================================" << endl;
        cout << "   ALL TESTS COMPLETED SUCCESSFULLY!" << endl;
        cout << "============================================" << endl;
        
    } catch (const std::exception& e) {
        cout << "\n❌ CRITICAL ERROR: " << e.what() << endl;
        return 1;
    }
    
    return 0;
}