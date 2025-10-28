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

        virtual double payoffPath( std::vector<double> d)
        {
            double stm = payoff(d[d.size()-1]);
            return stm;
        }
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

//====================================================== PART 3 ======================================================

class AsianOption : public Option{
    private:
        std::vector<double> time_steps;
    public:
        AsianOption(std::vector<double> ts , double e) : Option(e), time_steps(ts){}

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
            return payoff(sum/St.size());
        }

        bool isAsian() override
        {
            return true;
        }

};

class AsianCallOption : public AsianOption{
    private:
        double _strike;
    public:
        AsianCallOption(double e, double s, vector<double> ts) : AsianOption(ts, e) , _strike(s) {}

        double payoff(double underlying)
        {
            if(underlying>_strike)
                return underlying - _strike;
            else
                return 0;
        }
};

class AsianPutOption : public AsianOption{
    private:
        double _strike;
    public:
        AsianPutOption(double e, double s, vector<double> ts) : AsianOption(ts, e) , _strike(s) {}

        double payoff(double underlying)
        {
            if(underlying<_strike)
                return   _strike - underlying;
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
        int nbPaths = 0;
        double initial_price;
        double interest_rate;
        Option* option;
    public:
        BlackScholesMCPricer(Option* o, double ip, double ir, double vol) : option(o), initial_price(ip), interest_rate(ir){}
        double getNbPaths()
        {
            return nbPaths;
        }
        void generate(int nb_paths)
        {
            
        }

};


