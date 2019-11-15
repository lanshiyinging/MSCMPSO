//
// Created by lansy on 2019/11/12.
//

#include "function.h"

double Tablet(vector<double> x){
    double y;
    y = 10e6 * x[0] * x[0];
    for(int i = 1; i < dim; ++ i){
        y += x[i]*x[i];
    }
    return y;
}

double Quadric(vector<double> x){
    double y = 0;
    for(int i = 0; i < dim; ++ i){
        double tmp = 0;
        for(int j = 0; j < i; ++ j){
            tmp += x[j];
        }
        y += tmp*tmp;
    }
    return y;
}

double Rosenbrock(vector<double> x){
    double y = 0;
    for(int d = 0; d < dim-1; ++ d){
        y += 100*pow(pow(x[d+1]-x[d]*x[d], 2), 2) + pow((1-x[d]), 2);
    }
    return y;
}

double Griewank(vector<double> x){
    double y = 0;
    double y1 = 0, y2 = 0;
    for(int d = 0; d < dim; ++ d){
        y1 += x[d]*x[d];
        y2 *= cos(x[d]/sqrt(d));
    }
    y = y1/4000 - y2 + 1;
    return y;
}

double Rastrigin(vector<double> x){
    double y = 0;
    for(int d = 0; d < dim; ++ d){
        y += x[d] * x[d] - 10*cos(2*PI*x[d]) + 10;
    }
    return y;
}

double SchafferF7(vector<double> x){
    double y = 0;
    for(int d = 0; d < dim-1; ++ d){
        y += pow((x[d]*x[d]+x[d+1]*x[d+1]), 0.25)*(sin(50*pow((x[d]*x[d] + x[d+1]*x[d+1]), 0.1))+1.0);
    }
    return y;
}

