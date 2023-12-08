#pragma once
#include <vector>
#include <cmath>


using namespace std;


class filter
{
private:
    vector<vector<double>> H;
    vector<double> Hs;
    int ft;
public:
    filter();
    ~filter();
    filter(int nx, int ny, double r, int fit);
    void Filtering(vector<double>& xP, vector<double>& dfx);
    void Filtering(vector<double>& xP, vector<vector<double>>& dgdx);
    void UpdateFilter(vector<double>& xP, vector<double> x);
};