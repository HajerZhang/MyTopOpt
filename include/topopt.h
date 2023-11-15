#pragma once
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>

using namespace std;

class topopt
{
private:
    int nelx, nely, loop;
    double volfrac, penal, rmin, ft;
    double change;
    double ** x, ** xold, ** dc;
    double c;
    double KE[8][8];
    double **K, *F, *U;
    int * fixdofs, * alldofs, * freedofs;

public:
    topopt();
    ~topopt();
    topopt(int x, int y, double v, double pe, double r, int f);
    void sigmundopt();
    void check();
    void OC();
    void slovegauss(vector<vector<double>> &K, vector<double> &F, vector<double> &U);
    void elestiffmatix();
    void constrain();
    void output();
    bool solveLinearSystem(const vector<vector<double>>& A, const vector<double>& b, vector<double>& x);
    bool choleskyDecomposition(const vector<vector<double>>& A, vector<vector<double>>& L);
};