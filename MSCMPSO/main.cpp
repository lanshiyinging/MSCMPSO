#include <iostream>
#include <cmath>
#include "function.h"
#include "Macmpso.h"

using namespace std;


int main() {

    Macmpso m1(TabletWidth, Tablet);
    cout << "The result of Tablet: " << endl;
    m1.solution();
    //cout << "The result of Tablet: " << m1.bestFit << endl;

    /*
    Macmpso m2(QuadricWidth, Quadric);
    m2.solution();
    cout << "The result of Quadric: " << m2.bestFit << endl;

    Macmpso m3(RosenbrockWidth, Rosenbrock);
    m3.solution();
    cout << "The result of Rosenbrock: " << m3.bestFit << endl;

    Macmpso m4(GriewankWidth, Griewank);
    m4.solution();
    cout << "The result of Griewank): " << m4.bestFit << endl;

    Macmpso m5(RastriginWidth, Rastrigin);
    m5.solution();
    cout << "The result of Rastrigin: " << m5.bestFit << endl;

    Macmpso m6(SchafferWidth, SchafferF7);
    m6.solution();
    cout << "The result of SchafferF7: " << m6.bestFit << endl;
     */

    return 0;
}