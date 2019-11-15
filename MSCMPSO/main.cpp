#include <iostream>
#include <cmath>
#include "function.h"
#include "Macmpso.h"

using namespace std;


int main() {


    Macmpso m1(TabletWidth, Tablet);
    cout << "The result of Tablet: " << endl;
    m1.solution("tablet.txt");


    Macmpso m2(QuadricWidth, Quadric);
    cout << "The result of Quadric: " << endl;
    m2.solution("quadric.txt");


    Macmpso m3(RosenbrockWidth, Rosenbrock);
    cout << "The result of Rosenbrock: " << endl;
    m3.solution("rosenbrock.txt");


    Macmpso m4(GriewankWidth, Griewank);
    cout << "The result of Griewank: " << endl;
    m4.solution("griewank.txt");


    Macmpso m5(RastriginWidth, Rastrigin);
    cout << "The result of Rastrigin: " << endl;
    m5.solution("rastrigin.txt");


    Macmpso m6(SchafferWidth, SchafferF7);
    cout << "The result of SchafferF7: " << endl;
    m6.solution("schaffer.txt");


    return 0;
}