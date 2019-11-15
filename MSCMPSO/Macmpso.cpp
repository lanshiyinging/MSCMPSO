//
// Created by lansy on 2019/11/13.
//

#include <iostream>
#include "Macmpso.h"

Macmpso::Macmpso(double wid, double (*fc)(vector<double>)){
    f = fc;
    width = wid;
    e.seed(static_cast<unsigned int>(time(nullptr)));
}

void Macmpso::init() {
    x.clear();
    v.clear();
    pbest.clear();
    pgbest.clear();
    fit.clear();
    G.clear();
    T.clear();
    sigma.clear();
    uniform_real_distribution<double> u1(-width, width);
    uniform_real_distribution<double> u2(-1, 1);
    double minFit = MAX_DOUBLE;
    for(int i = 0; i < size; ++ i){
        vector<double> tempx;
        vector<double> tempv;
        double tmpFit;
        for(int d = 0; d < dim; ++ d){
            tempx.push_back(u1(e));
            tempv.push_back(u2(e));
        }
        tmpFit = f(tempx);
        fit.push_back(tmpFit);
        if(tmpFit < minFit){
            minFit = tmpFit;
            pgbest.assign(tempx.begin(), tempx.end());
        }
        x.push_back(tempx);
        v.push_back(tempv);
        G.push_back(0);
        T.push_back(0.5);
    }
    pbest.assign(x.begin(), x.end());
    bestFit = f(pgbest);
    for(int m = 0; m < M; ++ m){
        sigma.push_back(width*2);
    }
}

void Macmpso::updatePbest() {
    for(int i = 0; i < size; ++ i){
        if(fit[i] < f(pbest[i]))
            pbest[i].assign(x[i].begin(), x[i].end());
    }
}

void Macmpso::updatePgbest() {
    for(int i = 0; i < size; ++ i){
        if(f(pbest[i]) < f(pgbest))
            pgbest.assign(pbest[i].begin(), pbest[i].end());
    }
    bestFit = f(pgbest);
}


void Macmpso::updateV(double w) {
    uniform_real_distribution<double> u(0, 1);
    for(int i = 0; i < size; ++ i){
        for(int d = 0; d < dim; ++ d){
            double r1 = u(e);
            double r2 = u(e);
            v[i][d] = w*v[i][d] + c1*r1*(pbest[i][d]-x[i][d]) + c2*r2*(pgbest[i]-x[i][d]);
        }
    }
}

void Macmpso::escape() {
    normal_distribution<double> u;
    uniform_real_distribution<double> u2(0, 1);
    for(int d = 0; d < dim; ++ d){
        for(int i = 0; i < size; ++ i){
            //double r = u(e);
            if(v[i][d] < T[d]){
                double minf = MAX_DOUBLE;
                double randSigma = 0;
                for(int j = 0; j < M; ++ j){
                    double r = u(e);
                    double tmpf = f(addToX(x[i], r*sigma[j], d));
                    if(tmpf < minf){
                        minf = tmpf;
                        randSigma = r*sigma[j];
                    }
                }
                double r2 = u2(e);
                Vmax = width - abs(x[i][d]);
                if(minf < f(addToX(x[i], r2*Vmax, d)))
                    v[i][d] = randSigma;
                else
                    v[i][d] = r2*Vmax;
                G[d] += 1;
            }
        }
        if(G[d] > K1){
            G[d] = 0;
            T[d] /= K2;

        }
    }
}

void Macmpso::updatePos() {
    for(int i = 0; i < size; ++ i){
        for(int d = 0; d < dim; ++ d){
            x[i][d] += v[i][d];
        }
        fit[i] = f(x[i]);
    }
}

void Macmpso::updataSigma() {
    int pNum = size / M;
    vector<double> fits;
    fits.assign(fit.begin(), fit.end());
    sort(fits.begin(), fits.end());
    vector<double> fitX;
    double fitXsum = 0;
    double fitXmax = -(MAX_DOUBLE-1);
    double fitXmin = MAX_DOUBLE;
    for(int m = 0; m < M; ++ m){
        double fitsum = 0;
        for(int i = m*pNum; i < (m+1)*pNum; ++ i){
            fitsum += fits[i];
        }
        fitsum /= pNum;
        fitX.push_back(fitsum);
        fitXsum += fitsum;
        if(fitsum > fitXmax)
            fitXmax = fitsum;
        if(fitsum < fitXmin)
            fitXmin = fitsum;
    }
    double range = fitXmax - fitXmin;
    for(int m = 0; m < M; ++ m){
        sigma[m] *= exp((M*fitX[m]-fitXsum)/range);
        while(sigma[m] > width/4)
            sigma[m] = abs(width/4 - sigma[m]);
    }
}


void Macmpso::solution() {
    for(int i = 0; i < times; ++ i){
        init();
        double w;
        for(int j = 0; j < generation; ++ j){
            updatePbest();
            updatePgbest();
            w = wmax-(wmax-wmin)*j/6000;
            updateV(w);
            escape();
            updatePos();
            updataSigma();
            printf("iterator %d\tbest fitness: %lf\n", j, bestFit);
            //cout << "iterator " << j << "\tbest fitness: " << bestFit << endl;
            if(bestFit == 0){
                cout << "finished" << endl;
                break;
            }
        }
        cout << "train " << i << "\tresult: " << bestFit << endl;

    }

}

vector<double> Macmpso::addToX(vector<double> x1, double value, int d){
    vector<double> x2;
    x2.assign(x1.begin(), x1.end());
    x2[d] += value;
    return x2;
}


