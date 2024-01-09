#include "FiniteElement.h"
#include <iostream>
#include <cmath>
FiniteElement::FiniteElement()
{
}

FiniteElement::~FiniteElement()
{
}

FiniteElement::FiniteElement(double ym, double po, int nx, int ny)
{
    E = ym;
    nu = po;
    nelx = nx;
    nely = ny;
	F.resize(2 * (nely + 1) * (nelx + 1), 0.0);
	U.resize(2 * (nely + 1) * (nelx + 1), 0.0);
}

void FiniteElement::GetEleStiffMat()
{

	KE.resize(8, vector<double>(8, 0.0));
	double E = 1.0 / (1.0 - nu * nu);
	double lk[8];
	lk[0] = 0.5 - nu / 6;
	lk[1] = 1.0 / 8 + nu / 8;
	lk[2] = -1.0 / 4 - nu / 12;
	lk[3] = -1.0 / 8 + 3.0 * nu / 8;
	lk[4] = -1.0 / 4 + nu / 12;
	lk[5] = -1.0 / 8 - nu / 8;
	lk[6] = nu / 6;
	lk[7] = 1.0 / 8 - 3.0 * nu / 8;
	KE[0][0] = E * lk[0]; KE[0][1] = E * lk[1]; KE[0][2] = E * lk[2]; KE[0][3] = E * lk[3];KE[0][4] = E * lk[4]; KE[0][5] = E * lk[5]; KE[0][6] = E * lk[6]; KE[0][7] = E * lk[7];
	KE[1][0] = E * lk[1]; KE[1][1] = E * lk[0]; KE[1][2] = E * lk[7]; KE[1][3] = E * lk[6];KE[1][4] = E * lk[5]; KE[1][5] = E * lk[4]; KE[1][6] = E * lk[3]; KE[1][7] = E * lk[2];
	KE[2][0] = E * lk[2]; KE[2][1] = E * lk[7]; KE[2][2] = E * lk[0]; KE[2][3] = E * lk[5];KE[2][4] = E * lk[6]; KE[2][5] = E * lk[3]; KE[2][6] = E * lk[4]; KE[2][7] = E * lk[1];
	KE[3][0] = E * lk[3]; KE[3][1] = E * lk[6]; KE[3][2] = E * lk[5]; KE[3][3] = E * lk[0];KE[3][4] = E * lk[7]; KE[3][5] = E * lk[2]; KE[3][6] = E * lk[1]; KE[3][7] = E * lk[4];
	KE[4][0] = E * lk[4]; KE[4][1] = E * lk[5]; KE[4][2] = E * lk[6]; KE[4][3] = E * lk[7];KE[4][4] = E * lk[0]; KE[4][5] = E * lk[1]; KE[4][6] = E * lk[2]; KE[4][7] = E * lk[3];
	KE[5][0] = E * lk[5]; KE[5][1] = E * lk[4]; KE[5][2] = E * lk[3]; KE[5][3] = E * lk[2];KE[5][4] = E * lk[1]; KE[5][5] = E * lk[0]; KE[5][6] = E * lk[7]; KE[5][7] = E * lk[6];
	KE[6][0] = E * lk[6]; KE[6][1] = E * lk[3]; KE[6][2] = E * lk[4]; KE[6][3] = E * lk[1];KE[6][4] = E * lk[2]; KE[6][5] = E * lk[7]; KE[6][6] = E * lk[0]; KE[6][7] = E * lk[5];
	KE[7][0] = E * lk[7]; KE[7][1] = E * lk[2]; KE[7][2] = E * lk[1]; KE[7][3] = E * lk[4];KE[7][4] = E * lk[3]; KE[7][5] = E * lk[6]; KE[7][6] = E * lk[5]; KE[7][7] = E * lk[0];

	ReadyToAssemble();
}

void FiniteElement::BoundaryCondition()
{
	F [1] = -1;

    vector<int> fixdofs(nely + 2, 0.0);
	for (int i = 0; i <= nely; i++)
	{
		fixdofs[i] = 2 * i;   //+1
	}
	fixdofs[nely + 1] = 2 * (nelx + 1) * (nely + 1) - 1;   //+1
	
	vector<int> alldofs(2 * (nely + 1) * (nelx + 1), 0.0);
	for (int i = 0; i < 2 * (nely + 1) * (nelx + 1); i++)
	{
		alldofs[i] = i;
	}

	freedofs.resize(2 * nely * nelx + 2 * nelx + nely, 0.0);

	int i, j, k;
	i = 0; j = 0; k = 0;
	while (i < nely + 2 && j < 2 * (nely + 1) * (nelx + 1))
	{
		if (fixdofs[i] < alldofs[j]) {
			freedofs[k] = fixdofs[i];
			k++; i++;
		}
		else if (fixdofs[i] > alldofs[j]) {
			freedofs[k] = alldofs[j];
			k++; j++;
		}
		else if (fixdofs[i] == alldofs[j]) {
			i++;
			j++;
		}
	}

}

void FiniteElement::ReadyToAssemble()
{
	vector<vector<int>> nodenrs(nely+1, vector<int>(nelx+1, 0.0));
	vector<int> edofVec(nelx*nely);
	int k = 0;
    for(int i = 0; i < nelx+1; i++){
		for(int j = 0; j < nely+1; j++){
			k++;
			nodenrs[j][i] = k;
		}
	}

	{	
		for(int k = 0; k < nelx*nely; k++){
			edofVec[k] = 2 * nodenrs[k % nely][k / nely] + 1;
		}
	}

	edofMat.resize(nelx*nely, vector<int>(8, 0.0));
	{
		for(int k = 0; k < nelx*nely; k++){
			edofMat[k][0] = edofVec[k];
			edofMat[k][1] = edofVec[k] + 1;
			edofMat[k][2] = edofVec[k] + 2*nely + 2;
			edofMat[k][3] = edofVec[k] + 2*nely + 3;
			edofMat[k][4] = edofVec[k] + 2*nely;
			edofMat[k][5] = edofVec[k] + 2*nely + 1;
			edofMat[k][6] = edofVec[k] - 2;
			edofMat[k][7] = edofVec[k] - 1;
		}
	}

	iK.resize(64*nelx*nely, 0.0);
	jK.resize(64*nelx*nely, 0.0);
	{
		int el = 0;
		for(int k = 0; k < nelx*nely; k++){
			for(int i = 0; i < 8; i++){
				for(int j = 0; j < 8; j++){
					iK[el] = edofMat[k][j];
					jK[el] = edofMat[k][i];
					el++;
				}
			}
		}
	}

}

void FiniteElement::AssembleGloStiffMat(vector<double> xphys, double Ema, double Emi, double pe)
{
    vector<double> sK(64 * nelx * nely);
    for (int i = 0; i < nelx; ++i) {
        for (int j = 0; j < nely; ++j) {
            for (int k = 0; k < 8; ++k) {
                for (int l = 0; l < 8; ++l) {
                    int index = 64 * (i * nely + j) + 8 * k + l;
                    sK[index] = KE[k][l] * (Emi + pow(xphys[i * nely + j], pe) * (Ema - Emi));
                }
            }
        }
    }

	K.resize(2 * (nely + 1) * (nelx + 1), vector<double>(2 * (nely + 1) * (nelx + 1), 0.0));
	// Assemble the global stiffness matrix
	// cout << iK.size() << ' ' << jK.size() << ' ' << sK.size() << endl;
	// int maxiK = 1, maxjK = 1, miniK = 1, minjK = 1; 
	// for (int i = 0; i < 64 * nelx * nely; ++i) {
	// 	maxiK = max(maxiK, iK[i]);
	// 	maxjK = max(maxjK, jK[i]);
	// 	miniK = min(miniK, iK[i]);
	// 	minjK = min(minjK, jK[i]);
	// }
	// cout << maxiK << ' ' << maxjK << ' ' << miniK << ' ' << minjK << endl;
	// cout << K.size() << endl;
	// cout << K[0][0] << endl;
	for (size_t i = 0; i < K.size(); ++i) {
        for (size_t j = 0; j < K[i].size(); ++j) {
            K[i][j] = 0;
        }
    }

	int e1, e2;
	for (int i = 0; i < 64 * nelx * nely; ++i) {
		e1 = iK[i]-1;
		e2 = jK[i]-1;
		K[e1][e2] += sK[i];
	}
	// cout << K[0][0] << endl;

}

void FiniteElement::Solve()
{
	int k = freedofs.size();
	vector<vector<double>> K_freedofs(k, vector<double>(k));
	vector<vector<double>> K_freedofs_copy(k, vector<double>(k));
	for (int i = 0; i < k; i++) {
		for (int j = 0; j < k; j++) {
			K_freedofs[i][j] = K[freedofs[i]][freedofs[j]];
		}
	}
	vector<double> F_freedofs(k);
	for (int i = 0; i < k; i++) {
		F_freedofs[i] = F[freedofs[i]];
	}
	vector<double> U_freedofs(k);
	for (int i = 0; i < k; i++) {
		for (int j = 0; j < k; j++) {
			K_freedofs_copy[i][j] = (K_freedofs[i][j] + K_freedofs[j][i]) / 2.0;
		}
	}

	solveLinearSystem(K_freedofs_copy, F_freedofs, U_freedofs);
	for (int i = 0; i < k; i++) {
		U[freedofs[i]] = U_freedofs[i];
	}
	// for(int i=0; i< 40; i++){
	// 	cout << U[i] << endl;
	// }


}

void FiniteElement::slovegauss(vector<vector<double>>& A, vector<double>& b, vector<double>& x) 
{
	int n = b.size();

	for (int k = 0; k < n - 1; k++) {
		for (int i = k + 1; i < n; i++){
			double factor = A[i][k] / A[k][k];
			for (int j = k; j < n; j++) {
				A[i][j] -= factor * A[k][j];
			}
			b[i] -= factor * b[k];
		}
	}
	x[n - 1] = b[n - 1] / A[n - 1][n - 1];
	for (int i = n - 2; i >= 0; i--) {
		double sum = 0.0;
		for (int j = i + 1; j < n; j++) {
			sum += A[i][j] * x[j];
		}
		x[i] = (b[i] - sum) / A[i][i];
	}
}

bool FiniteElement::choleskyDecomposition(const vector<vector<double>>& A, vector<vector<double>>& L) {
    int n = A.size();
    L.resize(n, vector<double>(n, 0.0));

    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            double sum = 0.0;
            if (j == i) {
                for (int k = 0; k < j; k++) {
                    sum += L[j][k] * L[j][k];
                }
                L[j][j] = sqrt(A[j][j] - sum);
            } else {
                for (int k = 0; k < j; k++) {
                    sum += L[i][k] * L[j][k];
                }
                L[i][j] = (A[i][j] - sum) / L[j][j];
            }
        }
        if (L[i][i] <= 0.0) {

            return false;
        }
    }

    return true;
}

// Cholesky
bool FiniteElement::solveLinearSystem(const vector<vector<double>>& A, const vector<double>& b, vector<double>& x) {
    int n = A.size();
    vector<vector<double>> L;

    if (!choleskyDecomposition(A, L)) {
        cout << "Matrix is not symmetric positive definite." << endl;
        return false;
    }

    // Ly = b
    vector<double> y(n, 0.0);
    for (int i = 0; i < n; i++) {
        double sum = 0.0;
        for (int j = 0; j < i; j++) {
            sum += L[i][j] * y[j];
        }
        y[i] = (b[i] - sum) / L[i][i];
    }

    // L^T x = y
    x.resize(n, 0.0);
    for (int i = n - 1; i >= 0; i--) {
        double sum = 0.0;
        for (int j = i + 1; j < n; j++) {
            sum += L[j][i] * x[j];
        }
        x[i] = (y[i] - sum) / L[i][i];
    }

    return true;
}

void FiniteElement::GetUpdate(double& fx, vector<double>& gx, vector<double>& dfdx, vector<vector<double>>& dgdx,
		vector<double>& hessf, vector<vector<double>>& hessg, 
		vector<double>& xPhys, double Emin, double Emax, double penal, double volfrac, int itr)
{	
	fx = 0.0;
	for(int i = 0; i < nelx; i++){
		for(int j = 0; j < nely; j++){
			int e = i * nely + j;
			double ue[8];
			for(int k = 0; k < 8; k++){
				ue[k] = U[edofMat[e][k]-1];
			}
			double sum = 0.0;
			for(int k = 0; k < 8; k++){
				for (int l = 0; l < 8; l++){
					sum += ue[k] * KE[k][l] * ue[l];
				}
			}
			dfdx[e] = -penal * pow(xPhys[e], penal - 1) * (Emax - Emin) * sum;
			hessf[e] = 2 * pow((penal * (Emax - Emin) * pow(xPhys[e], penal - 1)), 2) / (Emax + (Emax-Emin)*pow(xPhys[e],penal))*sum;
			// cout << hessf[e] << endl;
			fx += pow(xPhys[e], penal) * (Emax - Emin) * sum;
		}
	}
	// volume constraint
	gx[0] = 0.0;
	for(int i = 0; i < nelx; i++){
		for(int j = 0; j < nely; j++){
			int e = i * nely + j;
			gx[0] += xPhys[e];
		}
	}
	gx[0] = gx[0] / (nelx * nely * volfrac) -1;
    for (int i = 0; i < nelx * nely; i++) {
        dgdx[0][i] = 1.0 ;
		hessg[0][i] = 0.0;
    }
	
}