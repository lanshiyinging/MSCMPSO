//
// Created by lansy on 2019/11/12.
//

#ifndef MSCMPSO_FUNCTION_H
#define MSCMPSO_FUNCTION_H

#include <cmath>
#include <vector>
using namespace std;

#define PI 3.14159265

const double TabletWidth = 100;
const double QuadricWidth = 100;
const double RosenbrockWidth = 50;
const double GriewankWidth = 300;
const double RastriginWidth = 5.12;
const double SchafferWidth = 100;
const int dim = 30;

double Tablet(vector<double> x);
double Quadric(vector<double> x);
double Rosenbrock(vector<double> x);
double Griewank(vector<double> x);
double Rastrigin(vector<double> x);
double SchafferF7(vector<double> x);



#endif //MSCMPSO_FUNCTION_H

