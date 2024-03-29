#pragma once
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <dirent.h>
#include <sys/stat.h>

using namespace std;

class topopt
{
private:

public:
    int nelx, nely;
    double volfrac;
    double penal;
    double rmin;
    int ft;
    int nele;
    int maxiter;
    double tolerance;

    double E;
    double nu;
    double Emax, Emin;

    vector<double> x;
    vector<double>  xold; 
    vector<double>  xPhys;


    double fx;
    vector<double> gx;
    vector<double> dfdx;
    vector<vector<double>> dgdx;
    vector<double> hessf;
    vector<vector<double>> hessg;

    // mmaparameters
    int n;
    int m;
    
    topopt();
    ~topopt();
    topopt(int nx, int ny, double v, double pe, double r, int f);
    void updateold(){
        for (int i = 0; i < n; i++)
        {
            xold[i] = x[i];
        }
    }
};

topopt::topopt()
{
}

topopt::~topopt()
{
}

topopt::topopt(int nx, int ny, double v, double pe, double r, int f)
{
    nelx = nx;
    nely = ny;
    volfrac=v;
    penal=pe;
    rmin=r;
    ft=f;
    nele = nelx*nely;

    x.resize(nele,volfrac);
    xold.resize(nele,volfrac);
    xPhys.resize(nele,volfrac);
    

    E = 1.0;
    nu = 0.3;

    Emin = 0.000000001;
    Emax = E;

    maxiter = 120;
    tolerance = 0.01;

    m = 1;
    n = nele;

    dfdx.resize(n, 0.0);
    hessf.resize(n, 0.0);
    gx.resize(m, 1.0);
    dgdx.resize(m, vector<double>(nele, 0));
    hessg.resize(m, vector<double>(nele, 0));
}

