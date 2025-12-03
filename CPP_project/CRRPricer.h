#ifndef CRRPRICER_H
#define CRRPRICER_H

#include "Option.h"
#include "BinaryTree.h"

class CRRPricer{
    private:    
        Option* option;
        int N;
        double asset_price;
        double up;
        double down;
        double interest_rate;
        BinaryTree<double>* _tree;
        bool computed;
        bool closed_form;
        BinaryTree<bool>* exercised;
        double max(double a, double b);
    public:
        CRRPricer(Option* option, int depth, double asset_price, double up, double down, double interest_rate);
        CRRPricer(Option* option, int depth, double asset_price, double r, double volatility);
        ~CRRPricer();
        bool getExercise(int a, int b);
        void compute();
        double get(int a, int b);
        double operator()(bool closed_form = false);
        void display_tree();
};

#endif