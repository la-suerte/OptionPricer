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