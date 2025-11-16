#ifndef MT_H
#define MT_H

#include <random>

class MT{
    private:
        static std::mt19937& get_generator();
    
    public:
        static double rand_unif();
        static double rand_norm();
};

#endif