#pragma once
#include <vector>

using namespace std;

class FiniteElement
{
private:
    double E; 
    double nu;
    vector<vector<double>> KE;
    int nelx, nely, nelz;
    vector<vector<int>> edofMat;
    vector<int> iK, jK;
    vector<vector<double>> K;
    vector<double> F, U;
    vector<int> freedofs;
    
public:
    FiniteElement();
    ~FiniteElement();
    FiniteElement(double ym, double po, int nx, int ny);
    void GetEleStiffMat();
    void BoundaryCondition();
    void ReadyToAssemble();
    void AssembleGloStiffMat(vector<double> xphys, double Ema, double Emi, double pe);
    void Solve();
    void solveLinearSystem();
    void slovegauss(vector<vector<double>>& A, vector<double>& b, vector<double>& x);
    bool choleskyDecomposition(const vector<vector<double>>& A, vector<vector<double>>& L);
    bool solveLinearSystem(const vector<vector<double>>& A, const vector<double>& b, vector<double>& x);
    void GetUpdate(double& fx, vector<double>& gx, vector<double>& dfdx, vector<vector<double>>& dgdx,
            vector<double>&hessf, vector<vector<double>>&hessg, 
            vector<double>& xPhys, double Emin, double Emax, double penal, double volfrac, int itr);
};

