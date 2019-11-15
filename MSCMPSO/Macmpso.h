//
// Created by lansy on 2019/11/13.
//

#ifndef MSCMPSO_MACMPSO_H
#define MSCMPSO_MACMPSO_H


#include <vector>
#include <random>
using namespace std;

class Macmpso {
public:
    const int times = 10;
    const int generation = 6000;
    const int size = 20;    //种群规模
    const int dim = 30;     //维度
    const double c1 = 1.4;
    const double c2 = 1.4;
    const int K1 = 5;
    const int K2 = 10;
    const int M = 5;
    const double MAX_DOUBLE = numeric_limits<double>::max();
    default_random_engine e;

    vector<vector<double>> x;
    vector<vector<double>> v;
    vector<vector<double>> pbest;
    vector<double> pgbest;
    vector<double> fit;
    vector<int> G;
    vector<double> T;
    vector<double> sigma;
    double (*f)(vector<double>);
    double width;
    double w = 0.95;
    double bestFit;
    double Vmax;

    Macmpso(double w, double (*fc)(vector<double>));
    void init();
    void updataPbest();
    void updataPgbest();
    void updateV(double w);
    void escape();
    void updatePos();
    void updataSigma();
    void solution();
    vector<double> addToX(vector<double> x1, double value, int d);


};


#endif //MSCMPSO_MACMPSO_H
