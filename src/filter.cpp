#include "filter.h"
#include <iostream>

filter::filter()
{
}


filter::~filter()
{
}

filter::filter(int nx, int ny, double rmin, int fit)
{
    int nele = nx * ny;
    H.resize(nele, vector<double>(nele, 0.0));
    Hs.resize(nele, 0.0);
    ft = fit;
    for (int i1 = 1; i1 < nx + 1; i1++)
    {
        for (int j1 = 1; j1 < ny + 1; j1++)
        {
            int e1 = (i1 - 1) * ny + j1;
            int i2min = fmax(i1 - (ceil(rmin)-1), 1);
            int i2max = fmin(i1 + (ceil(rmin)-1), nx);
            int j2min = fmax(j1 - (ceil(rmin)-1), 1);
            int j2max = fmin(j1 + (ceil(rmin)-1), ny);
            for (int i2 = i2min; i2 < i2max + 1; i2++)
            {
                for (int j2 = j2min; j2 < j2max + 1; j2++)
                {
                    int e2 = (i2 - 1) * ny + j2;
                    double fac = rmin - sqrt((i1 - i2)*(i1 - i2) + (j1 - j2)*(j1 - j2));
                    H[e1 - 1][e2 - 1] = fmax(0.0, fac);
                }
            }
        }
    }
    for (int i = 0; i < nele; i++)
    {
        for (int j = 0; j < nele; j++)
        {
            Hs[i] += H[i][j];
        }
    }
}

void filter::Filtering(vector<double>& xP, vector<double>& dfx)
{
    if(ft == 1){
        vector<double> dfxtemp = dfx;
        for (int i = 0; i < dfx.size(); i++)
        {
            dfx[i] = 0.0;
            for (int j = 0; j < dfx.size(); j++)
            {
                dfx[i] += H[i][j] * dfxtemp[j] / Hs[j];
                
            }
            // cout << dfx[i] << endl;
        }
        
    }
    else if(ft == 2){

    }
}

void filter::Filtering(vector<double>& xP, vector<vector<double>>& dgdx)
{
    if(ft == 1){
        // volume fraction
        vector<double> dfxtemp = dgdx[0];
        for (int i = 0; i < dgdx[0].size(); i++)
        {
            dgdx[0][i] = 0.0;
            for (int j = 0; j < dgdx[0].size(); j++)
            {
                dgdx[0][i] += H[i][j] * dfxtemp[j] / Hs[j];
                
            }
            // cout << dgdx[0][i] << endl;
        }
        
    }
    else if(ft == 2){

    }
}
void filter::UpdateFilter(vector<double>& xP, vector<double> x)
{
    double sum;
    for(int i = 0; i < x.size(); i++){
        sum = 0;
        for(int j = 0; j < x.size(); j++){

            sum += H[i][j] * x[j];
        }
        xP[i] = sum / Hs[i];
    }
}