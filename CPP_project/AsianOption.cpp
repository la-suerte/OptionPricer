#include "AsianOption.h"
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