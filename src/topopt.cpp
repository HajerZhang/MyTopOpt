#include "topopt.h"
#include <fstream>
#include <sstream>
#include <filesystem>

using namespace std;

topopt::topopt() {}
topopt::~topopt() {}

topopt::topopt(int nx,int ny, double v, double pe, double r,int f)
{
    nelx = nx;
    nely = ny;
    volfrac = v;
    penal = pe;
    rmin = r;
    ft = f;
    
    x = new double*[nely];
    xold = new double*[nely];
    dc = new double*[nely];
    for (int i = 0; i < nely; i++)
    {
        x[i] = new double[nelx];
        xold[i] = new double[nelx];
        dc[i] = new double[nelx];
        for (int j = 0; j < nelx; j++)
        {
            x[i][j] = volfrac;
            xold[i][j] = volfrac;
            dc[i][j] = 0.0;
        }
    }
    K = new double*[2 * (nelx + 1)*(nely + 1)];
    for (int i = 0; i < 2 * (nelx + 1)*(nely + 1); i++)
    {
        K[i] = new double[2 * (nelx + 1)*(nely + 1)];
        for (int j = 0; j < 2 * (nelx + 1)*(nely + 1); j++)
        {
            K[i][j] = 0.0;
        }
    }
    F = new double[2 * (nelx + 1)*(nely + 1)];
    for (int i = 0; i < 2 * (nelx + 1)*(nely + 1); i++)
    {
        F[i] = 0.0;
    }
    U = new double[2 * (nelx + 1)*(nely + 1)];
    for (int i = 0; i < 2 * (nelx + 1)*(nely + 1); i++)
    {
        U[i] = 0.0;
    }
    loop = 0;
    change = 1.0;
    elestiffmatix();
    constrain();
}

void topopt::elestiffmatix()
{
    double nu = 0.3;
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
}

void topopt::constrain()
{
    F [1] = -1;
    fixdofs = new int[nely + 2];
	for (int i = 0; i <= nely; i++)
	{
		fixdofs[i] = 2 * i;   //+1
	}
	fixdofs[nely + 1] = 2 * (nelx + 1) * (nely + 1) - 1;   //+1
	
	alldofs = new int[2 * (nely + 1) * (nelx + 1)];
	for (int i = 0; i < 2 * (nely + 1) * (nelx + 1); i++)
	{
		alldofs[i] = i;
	}

	freedofs = new int[2 * nelx * nely + 2 * nelx + nely];

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

void topopt::sigmundopt()
{
    while(change > 0.01)
    {
        loop++;
        for(int i=0; i < nely; i++)
        {
            for(int j=0; j < nelx; j++)
                xold[i][j] = x[i][j];
            
        }
        int n1, n2;
		int edof[8];
		
		for (int i = 0; i < 2 * (nelx + 1) * (nely + 1); i++)
		{
			for (int j = 0; j < 2 * (nelx + 1) * (nely + 1); j++)
				K[i][j] = 0.0;
		}
		for (int elx = 1; elx <= nelx; elx++)
		{
			for (int ely = 1; ely <= nely; ely++)
			{
				n1 = (nely + 1) * (elx - 1) + ely;
				n2 = (nely + 1) * elx + ely;
				
				edof[0] = 2 * n1 - 1 - 1; edof[1] = 2 * n1 - 1;
				edof[2] = 2 * n2 - 1 - 1; edof[3] = 2 * n2 - 1;
				edof[4] = 2 * n2 + 1 - 1; edof[5] = 2 * n2 + 2 - 1;
				edof[6] = 2 * n1 + 1 - 1; edof[7] = 2 * n1 + 2 - 1;
	
				for (int i = 0; i < 8; i++) {
					for (int j = 0; j < 8; j++) {
						K[edof[i]][edof[j]] += pow(x[ely - 1][elx - 1], penal) * KE[i][j];
					}
				}
			}
		}

        int k = 2 * nelx * nely + 2 * nelx + nely;
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
        // slovegauss(K_freedofs, F_freedofs, U_freedofs);

        for (int i = 0; i < k; i++) {
			U[freedofs[i]] = U_freedofs[i];
		}
		//fixeddofs = new int[nely + 2];
		for (int i = 0; i < nely + 2; i++) {
			U[fixdofs[i]] = 0;
		}
		//for (int i = 0; i < 2 * (nelx + 1) * (nely + 1); i++)
		//{
		//	std::cout << U_2D[i] << "\t";
		//}
		c = 0.0;

		double Ue[8];
		for (int elx = 1; elx <= nelx; elx++)
		{
			for (int ely = 1; ely <= nely; ely++)
			{
				n1 = (nely + 1) * (elx - 1) + ely;
				n2 = (nely + 1) * elx + ely;
				Ue[0] = U[2 * n1 - 2]; Ue[1] = U[2 * n1 - 1];
				Ue[2] = U[2 * n2 - 2]; Ue[3] = U[2 * n2 - 1];
				Ue[4] = U[2 * n2]; Ue[5] = U[2 * n2 + 1];
				Ue[6] = U[2 * n1]; Ue[7] = U[2 * n1 + 1];

				double sum1, sum2;
				sum2 = 0.0;
				for (int i = 0; i < 8; i++) {
					sum1 = 0.0;
					for (int j = 0; j < 8; j++) {
						sum1 += Ue[j] * KE[i][j];
					}
					sum2 += sum1 * Ue[i];
				}
				c += pow(x[ely - 1][elx - 1], penal) * sum2;
				dc[ely - 1][elx - 1] = -penal * pow(x[ely - 1][elx - 1], penal - 1) * sum2;
			}
		}

		check();
		OC();

		double mx = 0.0;
		for (int i = 0; i < nely; i++){
			for (int j = 0; j < nelx; j++) {
				mx = max(mx, fabs(x[i][j] - xold[i][j]));
			}
		}
		change = mx;
        
        cout << "loop:" << loop <<" change:" << change << endl;
        if(ft==1) output();
    }
	
}

void topopt::slovegauss(vector<vector<double>>& A, vector<double>& b, vector<double>& x) 
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

void topopt::check()
{
	double** dcn;
	dcn = new double* [nely];
	for (int i = 0; i < nely; i++) {
		dcn[i] = new double[nelx];
		for (int j = 0; j < nelx; j++){
			dcn[i][j] = 0.0;
		}
	}
	double sum, fac;
	for (int i = 1; i <= nelx; i++){
		for (int j = 1; j <= nely; j++){
			sum = 0.0;
			for (int k = max(i - (int)floor(rmin), 1); k <= min(i + (int)floor(rmin), nelx); k++)
			{
				for (int m = max(j - (int)floor(rmin), 1); m <= min(j + (int)floor(rmin), nely); m++)
				{
					fac = rmin - sqrt(pow(i - k, 2) + pow(j - m, 2));
					sum = sum + max(0.0, fac);

					dcn[j - 1][i - 1] += max(0.0, fac) * x[m - 1][k - 1] * dc[m - 1][k - 1];
				}
			}
			dcn[j - 1][i - 1] = dcn[j - 1][i - 1] / (x[j - 1][i - 1] * sum);
		}
	}
	for (int i = 0; i < nely; i++){
		for (int j = 0; j < nelx; j++){
			dc[i][j] = dcn[i][j];
		}
	}
	for (int i = 0; i < nely; i++){
		delete[] dcn[i];
	}
	delete[] dcn;
}

void topopt::OC()
{
	double mm, nn, move, lmid, sum;
	double** xnew;
	xnew = new double* [nely];
	for (int i = 0; i < nely; i++) {
		xnew[i] = new double[nelx];
		for (int j = 0; j < nelx; j++) {
			xnew[i][j] = 0.0;
		}
	}

	mm = 0.0; nn = 100000; move = 0.2;
	while (nn - mm > 0.0001)
	{
		lmid = 0.5 * (mm + nn);
		sum = 0.0;
		for (int i = 0; i < nely; i++) {
			for (int j = 0; j < nelx; j++) {
				xnew[i][j] = max(0.001, max(x[i][j] - move, min(1.0, min(x[i][j] + move, x[i][j] * sqrt(-dc[i][j] / lmid)))));
				sum += xnew[i][j];
			}
		}
		if (sum - volfrac * nelx * nely > 0) {
			mm = lmid;
		}
		else{
			nn = lmid;
		}
	}
	for (int i = 0; i < nely; i++) {
		for (int j = 0; j < nelx; j++) {
			x[i][j] = xnew[i][j];
		}
	}
	for (int i = 0; i < nely; i++) {
		delete[] xnew[i];
	}
	delete[] xnew;
}

void topopt::output()
{

	string filename = "./data/topopt_"+to_string(loop)+".txt";
	ofstream outputFile(filename);

    if (!outputFile.is_open()) {
        cerr << "Error opening file: " << filename << std::endl;
        return;
    }
	
    for (int i = 0; i < nely; ++i) {
        for (int j = 0; j < nelx; ++j) {
            outputFile << x[i][j] << " ";
        }
        outputFile << endl; 
    }

    outputFile.close();
}

bool topopt::choleskyDecomposition(const vector<vector<double>>& A, vector<vector<double>>& L) {
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
bool topopt::solveLinearSystem(const vector<vector<double>>& A, const vector<double>& b, vector<double>& x) {
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