#pragma once
#include <vector>
#include <cmath>


using namespace std;


class filter
{
private:
    vector<vector<double>> H;
    vector<double> Hs;
    vector<double> xTilde;
    vector<double> dx;
    int ft;
    int beta; // heaviside parameter
    
public:
    int loopbeta; // loop number for heaviside parameter
    filter();
    ~filter();
    filter(vector<double>& x, vector<double>& xP, int nx, int ny, double r, int fit);
    void Filtering(vector<double>& x, vector<double>& dfx, vector<vector<double>>& dgdx);
    void UpdateFilter(vector<double>& xP, vector<double> x);
    void Heaviside(double& change, double tol);
};