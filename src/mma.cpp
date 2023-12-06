#include "mma.h"
#include <iostream>
mma::mma()
{
}

mma::~mma()
{
}

mma::mma(int nn, int mm, vector<double> x)
{   
    double ai = 0.0;
    double ci = 1000.0;
    double di = 0.0;

    n = nn;
    m = mm;

    xmin.resize(n,0.0);
    xmax.resize(n,1.0);
    xold1 = x;
    xold2 = x;
    upp.resize(n,1.0);
    low.resize(n,1.0);
    a0 = 1.0;
    a.resize(m, ai);
    c.resize(m, ci);
    d.resize(m, di);

    epsimin = sqrt(n + m) * 0.000000001;
    feps = 0.000001;
    albefa = 0.1;
    asyinit = 0.5;
    asyincr = 1.2;

    alpha.resize(n, 0.0);
    beta.resize(n, 0.0);
    p0.resize(n, 0.0);
    q0.resize(n, 0.0);
    P.resize(m, vector<double>(n, 0.0));
    Q.resize(m, vector<double>(n, 0.0));
    b.resize(m, 0.0);
}

double mma::Update(int iter, vector<double> &xval, const double &f0val,const vector<double> &fval, 
            const vector<double> &df0dx, const vector<vector<double>> &dfdx, 
            const vector<double> &df0dx2, const vector<vector<double>> &dfdx2)
{   
    double change = 0.0;

    Setasymptotes(iter, xval);
    Setbounds(xval);
    GenSub(xval, f0val, fval, df0dx, dfdx, df0dx2, dfdx2);

    for(int i = 0; i < n; i++){
        xold2[i] = xold1[i];
        xold1[i] = xval[i];
    }

    SolveSub(xval);

    change = fabs(xval[0] - xold1[0]);
    for(int i = 0; i < n; i++){
        change = fmax(change, fabs(xval[i] - xold1[i]));
    }

    return change;
}

void mma::Setasymptotes(int iter, const vector<double>& x)
{
    if(iter <= 2){
        for(int i = 0; i < n; i++){
            low[i] = x[i] - asyinit * (xmax[i] - xmin[i]);
            upp[i] = x[i] + asyinit * (xmax[i] - xmin[i]);
        }
    }
    else if(iter > 3 ){
        for(int i = 0; i < n; i++){
            double factor = (x[i] - xold1[i]) * (xold1[i] - xold2[i]);
            if(factor > 0.0){
                factor = asyincr;
            }
            else if(factor < 0.0){
                factor = asydecr;
            }
            else{
                factor = 1.0;
            }
            low[i] = x[i] - factor * (xold1[i] - low[i]);
            upp[i] = x[i] + factor * (upp[i] - xold1[i]);
        }
    }
}

void mma::Setbounds(const vector<double>& x){
    double xxx;
    for(int i = 0; i < n; i++){
        xxx = low[i] + albefa * (x[i] - low[i]);
        alpha[i] = fmax(xxx, xmin[i]);
        xxx = upp[i] - albefa * (upp[i] - x[i]);
        beta[i] = fmin(xxx, xmax[i]);
    }
}

void mma::GenSub(vector<double> &xval, const double &f0val,const vector<double> &fval, 
            const vector<double> &df0dx, const vector<vector<double>> &dfdx, 
            const vector<double> &df0dx2, const vector<vector<double>> &dfdx2)
{   
    vector<double> ux1(n, 0.0);
    vector<double> ux2(n, 0.0);
    vector<double> ux3(n, 0.0);
    vector<double> xl1(n, 0.0);
    vector<double> xl2(n, 0.0);
    vector<double> xl3(n, 0.0);
    vector<double> ul1(n, 0.0);
    vector<double> ulinv1(n, 0.0);
    vector<double> uxinv1(n, 0.0);
    vector<double> xlinv1(n, 0.0);
    vector<double> uxinv3(n, 0.0);
    vector<double> xlinv3(n, 0.0);
    vector<double> diap(n, 0.0);
    vector<double> diaq(n, 0.0);

    for(int i = 0; i < n; i++){
        ux1[i] = upp[i] - xval[i];
        ux2[i] = ux1[i] * ux1[i];
        ux3[i] = ux2[i] * ux1[i];
        xl1[i] = xval[i] - low[i];
        xl2[i] = xl1[i] * xl1[i];
        xl3[i] = xl2[i] * xl1[i];
        ul1[i] = upp[i] - low[i];
        ulinv1[i] = 1.0 / ul1[i];
        uxinv1[i] = 1.0 / ux1[i];
        xlinv1[i] = 1.0 / xl1[i];
        uxinv3[i] = 1.0 / ux3[i];
        xlinv3[i] = 1.0 / xl3[i];
        diap[i] = (ux3[i] * xl1[i]) / (2 * ul1[i]);
        diaq[i] = (ux1[i] * xl3[i]) / (2 * ul1[i]);
    }

    vector<double> dg0dx2(n, 0.0);
    double del0 = 0.0;
    for(int i = 0; i < n; i++){
        p0[i] = 0.0;
        q0[i] = 0.0;
        if(df0dx[i] > 0.0){
            p0[i] = df0dx[i];
        }
        else if(df0dx[i] < 0.0){
            q0[i] = -df0dx[i];
        }
        p0[i] = p0[i] + 0.001 * abs(df0dx[i]) + feps * ulinv1[i];
        p0[i] = p0[i] * ux2[i];
        q0[i] = q0[i] + 0.001 * abs(df0dx[i]) + feps * ulinv1[i];
        q0[i] = q0[i] * xl2[i];
        dg0dx2[i] = 2*(p0[i] / ux3[i] + q0[i] / xl3[i]);
        del0 = df0dx2[i] - dg0dx2[i];
        if(del0 > 0.0){
            p0[i] = p0[i] + del0 * diap[i];
            q0[i] = q0[i] + del0 * diaq[i];
        }

    }

    vector<vector<double>> dgdx2(m, vector<double>(n, 0.0));
    double del = 0.0;
    for(int i = 0; i < m; i++){
        for(int j = 0; j < n; j++){
            P[i][j] = 0.0;
            Q[i][j] = 0.0;
            if(dfdx[i][j] > 0.0){
                P[i][j] = dfdx[i][j];
            }
            else if(dfdx[i][j] < 0.0){
                Q[i][j] = -dfdx[i][j];
            }
            P[i][j] = P[i][j] * ux2[j];
            Q[i][j] = Q[i][j] * xl2[j];
            dgdx2[i][j] = 2 * (P[i][j] / uxinv3[j] + Q[i][j] / xlinv3[j]);
            del = dfdx2[i][j] - dgdx2[i][j];
            if(del > 0.0){
                P[i][j] = P[i][j] + del * diap[j];
                Q[i][j] = Q[i][j] + del * diaq[j];
            }
        }
    }

    for(int i = 0; i < m; i++){
        b[i] = 0.0;
        for(int j = 0; j < n; j++){
            b[i] += P[i][j] * uxinv1[j] - Q[i][j] / xlinv1[j] - fval[i];
        }
    }
}

void mma::SolveSub(vector<double>& xval) 
{
	// double tidnewton = 0;
    double epsi = 1.0;
	vector<double> epsvecn(n, epsi * 1.0);
	vector<double> epsvecm(m, epsi * 1.0);
    vector<double> x(n, 0.0);
    vector<double> y(m, 1.0);
    double z = 1.0;
    vector<double> lam(m, 1.0);
    vector<double> xsi(n, 1.0);
    vector<double> eta(n, 1.0);
    vector<double> mu(m, 1.0);
    double zet = 1.0;
    vector<double> s(m, 1.0);
    vector<double> ux1(n, 0.0);
	vector<double> xl1(n, 0.0);
	vector<double> ux2(n, 0.0);
	vector<double> xl2(n, 0.0);
	vector<double> uxinv1(n, 0.0);
	vector<double> xlinv1(n, 0.0);
	vector<double> plam(n, 0.0);
	vector<double> qlam(n, 0.0);
    vector<double> gvec(m, 0.0);
    vector<double> dpsidx(n, 0.0);
    vector<double> rex(n, 0.0);
    vector<double> rey(m, 0.0);
    double rez;
    vector<double> relam(m, 0.0);
    vector<double> rexsi(n, 0.0);
    vector<double> reeta(n, 0.0);
    vector<double> remu(m, 0.0);
    double rezet;
    vector<double> res(m, 0.0);
    vector<double> residu1(m + n + 1, 0.0);
    vector<double> residu2(3 * m + 2 * n + 1, 0.0);
    vector<double> residu(4 * m + 3 * n + 2, 0.0);
    double residunorm;
    double residumax;

    vector<double> ux3(n, 0.0);
    vector<double> xl3(n, 0.0);
    vector<double> uxinv2(n, 0.0);
    vector<double> xlinv2(n, 0.0);
    vector<vector<double>> GG(m, vector<double>(n, 0.0));

    vector<double> delx(n, 0.0);
    vector<double> dely(m, 0.0);
    double delz = 0.0;
    vector<double> dellam(m, 0.0);
    vector<double> diagx(n, 0.0);
    vector<double> diagxinv(n, 0.0);
    vector<double> diagy(m, 0.0);
    vector<double> diagyinv(m, 0.0);
    vector<double> diaglam(m, 0.0);
    vector<double> diaglamyi(m, 0.0);
    vector<double> diaglamyiinv(m, 0.0);
    vector<double> dellamyi(m, 0.0);

    vector<vector<double>> Axx(n, vector<double>(n , 0.0));
    vector<vector<double>> azz(m, vector<double>(m , 0.0));
    vector<vector<double>> axz(n, vector<double>(m , 0.0));
    vector<double> bx(n, 0.0);
    vector<double> bz(m, 0.0);

    vector<vector<double>> AA(m + n, vector<double>(m + n , 0.0));
    vector<double> bb(m + n, 0.0);

    vector<double> result(m + n, 0.0);

    vector<double> dx(n, 0.0);
    double dz = 0.0;
    vector<double> dy(m, 0.0);
    vector<double> dlam(m, 0.0);
    vector<double> dxsi(n, 0.0);
    vector<double> deta(n, 0.0);
    vector<double> dmu(m, 0.0);
    double dzet = 0.0;
    vector<double> ds(m, 0.0);
    vector<double> xx(5 * m + 2 * n + 1, 0.0);
    vector<double> dxx(5 * m + 2 * n + 1, 0.0);

    vector<double> stepxx(5 * m + 2 * n + 1, 0.0);
    double stmxx;
    vector<double> stepalpha(n , 0.0);
    double stmalpha;
    vector<double> stepbeta(n , 0.0);
    double stmbeta;
    double stmalbe;
    double stmalbexx;
    double stminv;
    double steg;

    vector<double> xold(n, 0.0);
    vector<double> yold(m, 0.0);
    double zold = 0.0;
    vector<double> lamold(m, 0.0);
    vector<double> xsiold(n, 0.0);
    vector<double> etaold(n, 0.0);
    vector<double> muold(m, 0.0);
    double zetold = 0.0;
    vector<double> sold(m, 0.0);

    int itto;
    double resinew;


    for (int i = 0; i < n; ++i) {
        x[i] = 0.5 * (alpha[i] + beta[i]);
        xsi[i] = fmax(1.0, 1.0/ (x[i] - alpha[i]));
        eta[i] = fmax(1.0, 1.0/ (beta[i] - x[i]));
    }
    for(int i = 0; i < m; ++i){
        mu[i] = fmax(1.0, 0.5 * c[i]);
    }
    int itera = 0;
	
    while (epsi > epsimin){
		for(int i = 0; i < n; ++i) epsvecn[i] = epsi * 1.0;
		for(int i = 0; i < m; ++i) epsvecm[i] = epsi * 1.0;
		for(int i = 0; i < n; ++i){
			ux1[i] = upp[i] - x[i];
			xl1[i] = x[i] - low[i];
			ux2[i] = ux1[i] * ux1[i];
			xl2[i] = xl1[i] * xl1[i];
			uxinv1[i] = 1.0 / ux1[i];
			xlinv1[i] = 1.0 / xl1[i];
		}
		for(int i = 0; i < n; ++i){
			double sum1 = 0.0, sum2 = 0.0;
			for(int j = 0; j < m; ++j){
				sum1 += P[j][i] * lam[j];
				sum2 += Q[j][i] * lam[j];
			}
			plam[i] =  p0[i] + sum1;
			qlam[i] =  q0[i] + sum2;
    	}

		for(int i = 0; i < m; i++){
            gvec[i] = 0.0;
            for(int j = 0; j < n; j++){
                gvec[i] += P[i][j] * uxinv1[j] + Q[i][j] / xlinv1[j];
            }
        }
        for(int i = 0; i < n; i++){
            dpsidx[i] = plam[i] / ux2[i] - qlam[i] / xl2[i];
        }
        for(int i = 0; i < n; i++) rex[i] = dpsidx[i] - xsi[i] + eta[i];
        for(int i = 0; i < m; i++) rey[i] = c[i] + d[i] * y[i] - mu[i] - lam[i];
        rez = a0 - zet;
        for(int i = 0; i < m; i++) rez += a[i] * lam[i];
        for(int i = 0; i < m; i++){
            relam[i] = gvec[i] - z * a[i] - y[i] + s[i] - b[i];
        }
        for(int i = 0; i < n; i++){
            rexsi[i] = xsi[i] * (x[i] - alpha[i]) - epsvecn[i];
            reeta[i] = eta[i] * (beta[i] - x[i]) - epsvecn[i];
        }
        for(int i = 0; i < m; i++){
            remu[i] = mu[i] * y[i] - epsvecm[i];
        }
        rezet = zet * z - epsi;
        for(int i = 0; i < m; i++){
            res[i] = lam[i] * s[i] - epsvecm[i];
        }
        
        for(int i = 0; i < n; i++){
            residu1[i] = rex[i];
        }
        for(int i = 0; i < m; i++){
            residu1[n + i] = rey[i];
        }
        residu1[n + m] = rez;

        for(int i = 0; i < m; i++){
            residu2[i] = relam[i];
        }
        for(int i = 0; i < n; i++){
            residu2[m + i] = rexsi[i];
        }
        for(int i = 0; i < n; i++){
            residu2[m + n + i] = reeta[i];
        }
        for(int i = 0; i < m; i++){
            residu2[2 * n + m + i] = remu[i];
        }
        residu2[2 * n + 2 * m] = rezet;
        for(int i = 0; i < m; i++){
            residu2[2 * n + 2 * m + 1 + i] = res[i];
        }

        for(int i = 0; i < n + m + 1; i++){
            residu[i] = residu1[i];
        }
        for(int i = 0; i < 3 * m + 2 * n + 1; i++){
            residu[i + n + m + 1] = residu2[i];
        }
        residunorm = 0.0;
        for(int i = 0; i < 4 * m + 3 * n + 2; i++){
            residunorm += residu[i] * residu[i];
        }
        residunorm = sqrt(residunorm);
        residumax = residu[0];
        for(int i = 0; i < residu.size(); i++){
            residumax = fmax(residumax, fabs(residu[i]));
        }
        int ittt = 0;
        while (residumax > 0.9 * epsi && ittt < 100){
            ittt++;
            itera++;
            for(int i = 0; i < n; i++){
                ux1[i] = upp[i] - x[i];
                xl1[i] = x[i] - low[i];
                ux2[i] = ux1[i] * ux1[i];
                xl2[i] = xl1[i] * xl1[i];
                ux3[i] = ux2[i] * ux1[i];
                xl3[i] = xl2[i] * xl1[i];
                uxinv1[i] = 1.0 / ux1[i];
                xlinv1[i] = 1.0 / xl1[i];
                uxinv2[i] = 1.0 / ux2[i];
                xlinv2[i] = 1.0 / xl2[i];
            }
            for(int i = 0; i < n; ++i){
                double sum1 = 0.0, sum2 = 0.0;
                for(int j = 0; j < m; ++j){
                    sum1 += P[j][i] * lam[j];
                    sum2 += Q[j][i] * lam[j];
                }
                plam[i] =  p0[i] + sum1;
                qlam[i] =  q0[i] + sum2;
    	    }
            for(int i = 0; i < m; i++){
                gvec[i] = 0.0;
                for(int j = 0; j < n; j++){
                    gvec[i] += P[i][j] * uxinv1[j] + Q[i][j] / xlinv1[j];
                }
            }
            for(int i = 0; i < m; i++){
                for(int j = 0; j < n; j++){
                    GG[i][j] = P[i][j] * uxinv2[j] + Q[i][j] * xlinv2[j];
                }
            }
            for(int i = 0; i < n; i++){
                dpsidx[i] = plam[i] / ux2[i] - qlam[i] / xl2[i];
            }
            for(int i = 0; i < n; i++) 
            delx[i] = dpsidx[i] - epsvecn[i] / (x[i] - alpha[i]) + epsvecn[i] / (beta[i] - x[i]);
            for(int i = 0; i < m; i++)
            dely[i] = c[i] + d[i] * y[i] - lam[i] - epsvecm[i] / y[i];
            delz = a0 - epsi/z;
            for(int i = 0; i < m; i++)
            delz += a[i] * lam[i];
            for(int i = 0; i < m; i++){
                dellam[i] = gvec[i] - a[i] * z - y[i] - b[i] + epsvecm[i] / lam[i];
            }
            for(int i = 0; i < n; i++){
                diagx[i] = plam[i] / ux3[i] + qlam[i] / xl3[i];
                diagx[i] = 2 * diagx[i] + xsi[i] / (x[i] - alpha[i]) + eta[i] / (beta[i] - x[i]);
                diagxinv[i] = 1.0 / diagx[i];
            }
            for(int i = 0; i < m; i++){
                diagy[i] = d[i] + mu[i] / y[i];
                diagyinv[i] = 1.0 / diagy[i];
                diaglam[i] = s[i] / lam[i];
                diaglamyi[i] = diaglam[i] + diagyinv[i];
                diaglamyiinv[i] = 1.0 / diaglamyi[i];
            }
            for(int i = 0; i < m; i++){
                dellamyi[i] = dellam[i] + dely[i] / diagy[i];
            }
            for(int i = 0; i < n; i++){
                Axx[i][i] = diagx[i];
                for(int j = 0; j < n; j++){
                    for(int k = 0; k < m; k++){
                        Axx[j][i] += GG[k][i] * diaglamyiinv[k] * GG[k][j];
                    }    
                }
            }
            for(int i = 0; i < m; i++){
                azz[i][i] = zet / z;
                for(int j = 0; j < m; j++){
                    azz[j][i] += a[j] * a[j] /diaglamyiinv[j];
                }
            }
            for(int i = 0; i < n; i++){
                for(int j = 0; j < m; j++){
                    axz[i][j] = -GG[j][i] * a[j] / diaglamyi[j];
                }
            }
            for(int i = 0; i < n; i++){
                bx[i] = delx[i];
                for(int j = 0; j < m; j++){
                    bx[i] += GG[j][i] * dellamyi[j] / diaglamyi[j];
                }
            }
            for(int i = 0; i < m; i++){
                bz[i] = delz;
                for(int j = 0; j < m; j++){
                    bz[i] -= a[j] * dellamyi[j] / diaglamyi[j];
                }
            }
            for(int i = 0; i < n; i++){
                for(int j = 0; j < n; j++){
                    AA[i][j] = Axx[i][j];
                }
            }
            for(int i = 0; i < m; i++){
                for(int j = 0; j < n; j++){
                    AA[n + i][j] = axz[j][i];
                    AA[j][n + i] = axz[j][i];
                }
            }
            for(int i = 0; i < m; i++){
                for(int j = 0; j < m; j++){
                    AA[n + i][n + j] = azz[i][j];
                }
            }
            for(int i = 0; i < n; i++){
                bb[i] = -bx[i];
            }
            for(int i = 0; i < m; i++){
                bb[n + i] = -bz[i];
            }

            result = SolveLinearSystem(AA, bb);

            for(int i = 0; i < n; i++){
                dx[i] = result[i];
            }
            dz = result[n];
            for(int i = 0; i < m; i++){
                double sum = 0.0;
                for(int j = 0; j < n; j++){
                    sum += GG[i][j] * dx[j];
                }
                dlam[i] = sum / diaglamyi[i] 
                - dz * a[i] / diaglamyi[i] 
                + dellamyi[i] / diaglamyi[i];
            }
            for(int i = 0; i < m; i++){
                dy[i] = -dely[i] / diagy[i] + dlam[i] / diagy[i];
            }
            for(int i = 0; i < n; i++){
                dxsi[i] = -xsi[i] + epsvecn[i] / (x[i] - alpha[i]) - xsi[i] * dx[i] / (x[i] - alpha[i]);
                deta[i] = -eta[i] + epsvecn[i] / (beta[i] - x[i]) + eta[i] * dx[i] / (beta[i] - x[i]);
            }
            for (int i = 0; i < m; i++){
                dmu[i] = -mu[i] + epsvecm[i] / y[i] - mu[i] * dy[i] / y[i];
            }
            dzet = -zet + epsi / z - zet * dz / z;
            for(int i = 0; i < m; i++){
                ds[i] = -s[i] + epsvecm[i] / lam[i] - s[i] * dlam[i] / lam[i];
            }

            for(int i = 0; i < m; i++){
                xx[i] = y[i];
            }
            xx[m] = z;
            for(int i = 0; i < m; i++){
                xx[m + 1 + i] = lam[i];
            }
            for(int i = 0; i < n; i++){
                xx[2 * m + 1 + i] = xsi[i];
            }
            for(int i = 0; i < n; i++){
                xx[2 * m +  n + 1 + i] = eta[i];
            }
            for(int i = 0; i < m; i++){
                xx[2 * m + 2 * n + 1 + i] = mu[i];
            }
            xx[3 * m + 2 * n + 1] = zet;
            for(int i = 0; i < m; i++){
                xx[3 * m + 2 * n + 2 + i] = s[i];
            }
            
            for(int i = 0; i < m; i++){
                dxx[i] = dy[i];
            }
            dxx[m] = dz;
            for(int i = 0; i < m; i++){
                dxx[m + 1 + i] = dlam[i];
            }
            for(int i = 0; i < n; i++){
                dxx[2 * m + 1 + i] = dxsi[i];
            }
            for(int i = 0; i < n; i++){
                dxx[2 * m +  n + 1 + i] = deta[i];
            }
            for(int i = 0; i < m; i++){
                dxx[2 * m + 2 * n + 1 + i] = dmu[i];
            }
            dxx[3 * m + 2 * n + 1] = dzet;
            for(int i = 0; i < m; i++){
                dxx[3 * m + 2 * n + 2 + i] = ds[i];
            }
            
            for(int i = 0; i < 5 * m + 2 * n + 1; i++){
                stepxx[i] = -1.01 * dxx[i] / xx[i];
            }
            stmxx = stepxx[0];
            for(int i = 0; i < 5 * m + 2 * n + 1; i++){
                stmxx = fmax(stmxx, stepxx[i]);
            }

            for(int i = 0; i < n; i++){
                stepalpha[i] = -1.01 * dx[i] / (x[i] - alpha[i]);
            }
            stmalpha = stepalpha[0];
            for(int i = 0; i < n; i++){
                stmalpha = fmax(stmalpha, stepalpha[i]);
            }

            for(int i = 0; i < n; i++){
                stepbeta[i] = 1.01 * dx[i] / (beta[i] - x[i]);
            }
            stmbeta = stepbeta[0];
            for(int i = 0; i < n; i++){
                stmbeta = fmax(stmbeta, stepbeta[i]);
            }

            stmalbe = fmax(stmalpha, stmbeta);
            stmalbexx = fmax(stmalbe, stmxx);
            stminv = fmax(1.0, stmalbexx);
            steg = 1.0 / stminv;

            for(int i = 0; i < n; i++){
                xold[i] = x[i];
            }
            for(int i = 0; i < m; i++){
                yold[i] = y[i];
            }
            zold = z;
            for(int i = 0; i < m; i++){
                lamold[i] = lam[i];
            }
            for(int i = 0; i < n; i++){
                xsiold[i] = xsi[i];
            }
            for(int i = 0; i < n; i++){
                etaold[i] = eta[i];
            }
            for(int i = 0; i < m; i++){
                muold[i] = mu[i];
            }
            zetold = zet;
            for(int i = 0; i < m; i++){
                sold[i] = s[i];
            }
            
            itto = 0;
            resinew = 2 * residunorm;
            while( resinew > residunorm && itto < 50){
                itto++;

                for(int i = 0; i < n; i++){
                    x[i] = xold[i] + steg * dx[i];
                }
                for(int i = 0; i < m; i++){
                    y[i] = yold[i] + steg * dy[i];
                }
                z = zold + steg * dz;
                for(int i = 0; i < m; i++){
                    lam[i] = lamold[i] + steg * dlam[i];
                }
                for(int i = 0; i < n; i++){
                    xsi[i] = xsiold[i] + steg * dxsi[i];
                }
                for(int i = 0; i < n; i++){
                    eta[i] = etaold[i] + steg * deta[i];
                }
                for(int i = 0; i < m; i++){
                    mu[i] = muold[i] + steg * dmu[i];
                }
                zet = zetold + steg * dzet;
                for(int i = 0; i < m; i++){
                    s[i] = sold[i] + steg * ds[i];
                }
                for(int i = 0; i < n; i++){
                    ux1[i] = upp[i] - x[i];
                    xl1[i] = x[i] - low[i];
                    ux2[i] = ux1[i] * ux1[i];
                    xl2[i] = xl1[i] * xl1[i];
                    uxinv1[i] = 1.0 / ux1[i];
                    xlinv1[i] = 1.0 / xl1[i];
                }
                for(int i = 0; i < n; ++i){
                    double sum1 = 0.0, sum2 = 0.0;
                    for(int j = 0; j < m; ++j){
                        sum1 += P[j][i] * lam[j];
                        sum2 += Q[j][i] * lam[j];
                    }
                    plam[i] =  p0[i] + sum1;
                    qlam[i] =  q0[i] + sum2;
                }
                for(int i = 0; i < m; i++){
                    gvec[i] = 0.0;
                    for(int j = 0; j < n; j++){
                        gvec[i] += P[i][j] * uxinv1[j] + Q[i][j] / xlinv1[j];
                    }
                }
                for(int i = 0; i < n; i++){
                    dpsidx[i] = plam[i] / ux2[i] - qlam[i] / xl2[i];
                }

                for(int i = 0; i < n; i++) rex[i] = dpsidx[i] - xsi[i] + eta[i];
                for(int i = 0; i < m; i++) rey[i] = c[i] + d[i] * y[i] - mu[i] - lam[i];
                rez = a0 - zet;
                for(int i = 0; i < m; i++) rez += a[i] * lam[i];
                for(int i = 0; i < m; i++){
                    relam[i] = gvec[i] - z * a[i] - y[i] + s[i] - b[i];
                }
                for(int i = 0; i < n; i++){
                    rexsi[i] = xsi[i] * (x[i] - alpha[i]) - epsvecn[i];
                    reeta[i] = eta[i] * (beta[i] - x[i]) - epsvecn[i];
                }
                for(int i = 0; i < m; i++){
                    remu[i] = mu[i] * y[i] - epsvecm[i];
                }
                rezet = zet * z - epsi;
                for(int i = 0; i < m; i++){
                    res[i] = lam[i] * s[i] - epsvecm[i];
                }
                
                for(int i = 0; i < n; i++){
                    residu1[i] = rex[i];
                }
                for(int i = 0; i < m; i++){
                    residu1[n + i] = rey[i];
                }
                residu1[n + m] = rez;

                for(int i = 0; i < m; i++){
                    residu2[i] = relam[i];
                }
                for(int i = 0; i < n; i++){
                    residu2[m + i] = rexsi[i];
                }
                for(int i = 0; i < n; i++){
                    residu2[m + n + i] = reeta[i];
                }
                for(int i = 0; i < m; i++){
                    residu2[2 * n + m + i] = remu[i];
                }
                residu2[2 * n + 2 * m] = rezet;
                for(int i = 0; i < m; i++){
                    residu2[2 * n + 2 * m + 1 + i] = res[i];
                }

                for(int i = 0; i < n + m + 1; i++){
                    residu[i] = residu1[i];
                }
                for(int i = 0; i < 3 * m + 2 * n + 1; i++){
                    residu[i + n + m + 1] = residu2[i];
                }
                resinew = 0.0;
                for(int i = 0; i < 4 * m + 3 * n + 2; i++){
                    resinew += residu[i] * residu[i];
                }
                resinew = sqrt(resinew);
                steg = steg / 2.0;
            }

            residunorm = resinew;
            residumax = residu[0];
            for(int i = 0; i < 4 * m + 3 * n + 2; i++){
                residumax = fmax(residumax, fabs(residu[i]));
            }

            steg = 2 * steg;
        }
        epsi = 0.1 * epsi;
	}
    for(int i = 0; i < n; i++){
        xval[i] = x[i];
    }
    cout << 'o';
}

vector<double> mma::SolveLinearSystem(vector<vector<double>> A, vector<double> b) 
{
    vector<double> x(n);
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
    return x;
}