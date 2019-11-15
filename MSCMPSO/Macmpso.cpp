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
    for(int m = 0; m < M; ++ m){
        sigma.push_back(width*2);
    }
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
    for(int i = 0; i < size; ++ i){
        for(int d = 0; d < dim; ++ d){
            //double r = u(e);
            if(fabs(v[i][d]) < T[d]){
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
            x[i][d] += v[i][d];
            if(f(x[i]) < f(pbest[i]))
                pbest[i].assign(x[i].begin(), x[i].end());
            if(f(pbest[i]) < f(pgbest))
                pgbest.assign(pbest[i].begin(), pbest[i].end());
        }
    }
}

void Macmpso::updataSigma() {
    int pNum = size / M;
    vector<double> fits;
    for(int i = 0; i < size; ++ i){
        fits.push_back(f(x[i]));
    }
    //calFitness(fits);
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
    double range = fitXmax - fitXmin + pow(10, -10);
    for(int m = 0; m < M; ++ m){
        sigma[m] *= exp((M*fitX[m]-fitXsum)/range);
        while(sigma[m] > width/4)
            sigma[m] -= width/4;
    }
}


void Macmpso::updateGT(){
    for(int d = 0; d < dim; ++ d){
        if(G[d] > K1){
            G[d] = 0;
            T[d] = T[d]/K2;
        }
    }
}

void Macmpso::solution(string filename) {
    ofstream out;
    out.open(filename, ios::out);
    vector<double> eval;
    for(int i = 0; i < times; ++ i){
        init();
        double w;
        for(int j = 0; j < generation; ++ j){
            w = wmax-(wmax-wmin)*j/6000;
            updateV(w);
            escape();
            updataSigma();
            updateGT();
            if(i == 0)
                out << i << " " << log(f(pgbest)) << endl;
        }
        eval.push_back(f(pgbest));
        cout << scientific << "train " << i << "\t best fitness: " << f(pgbest) << "\tbest position: ";
        printBestPosition();
    }
    out.close();
    out.open("eval.txt", ios::app);
    double min = MAX_DOUBLE, max = -1, mean = 0, sd = 0;
    for(int t = 0; t < times; ++ t){
        mean += eval[t];
        if(eval[t] > max)
            max = eval[t];
        if(eval[t] < min)
            min = eval[t];
    }
    mean /= times;
    for(int t = 0; t < times; ++ t){
        sd += (eval[t]-mean)*(eval[t]-mean);
    }
    out << scientific << min << "\t" << max << "\t" << mean << "\t" << sd << endl;
}

vector<double> Macmpso::addToX(vector<double> x1, double value, int d){
    vector<double> x2;
    x2.assign(x1.begin(), x1.end());
    x2[d] += value;
    return x2;
}

void Macmpso::printBestPosition() {
    cout << "(";
    for(int d = 0; d < dim; ++ d){
        cout << scientific << pgbest[d] << " ";
    }
    cout << ")" << endl;
}


