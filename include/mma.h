#pragma once
#include <vector>
#include <cmath>

using namespace std;

class mma 
{
public:

	mma();
	~mma();
	mma(int nn, int mm, vector<double> x);
    double Update(int iter, vector<double> &xval, const double &f0val,const vector<double> &fval, 
            const vector<double> &df0dx, const vector<vector<double>> &dfdx, 
            const vector<double> &df0dx2, const vector<vector<double>> &dfdx2);
    void Setasymptotes(int iter, const vector<double>& x);
    void GenSub(vector<double> &xval, const double &f0val,const vector<double> &fval, 
            const vector<double> &df0dx, const vector<vector<double>> &dfdx, 
            const vector<double> &df0dx2, const vector<vector<double>> &dfdx2);
    void Setbounds(const vector<double>& x);
    void SolveSub(vector<double> &xval);
    vector<double> SolveLinearSystem(vector<vector<double>> A, vector<double> b);
    
    
private:
    //MMA parameters
    int m;
    int n;
    vector<double> xmin;
    vector<double> xmax;
    vector<double> xold1;
    vector<double> xold2;
    vector<double> low;
    vector<double> upp;    
    double a0;
    vector<double> a;
    vector<double> c;
    vector<double> d;

    double feps;
    double epsimin;
    double asyinit;
    double asyincr;
    double asydecr;
    double albefa;

    // subproblem parameters
    vector<double> alpha;
    vector<double> beta;
    vector<double> p0;
    vector<double> q0;
    vector<vector<double>> P;
    vector<vector<double>> Q;
    vector<double> b;
};
