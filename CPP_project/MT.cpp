#include "MT.h"

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