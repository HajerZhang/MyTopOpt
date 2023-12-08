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
};