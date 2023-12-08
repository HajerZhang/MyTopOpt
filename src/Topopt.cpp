#include "Topopt.h"

topopt::topopt()
{
}

topopt::~topopt()
{
}

topopt::topopt(int nx, int ny, double v, double pe, double r, int f)
    : nelx(nx)
    , nely(ny)
    , volfrac(v)
    , penal(pe)
    , rmin(r)
    , ft(f)
{
    nele = nelx * nely;

    x.resize(nele, volfrac);
    xold.resize(nele, volfrac);
    xPhys.resize(nele, volfrac);

    E = 1.0;
    nu = 0.3;

    Emin = 0.000000001;
    Emax = E;

    maxiter = 200;
    tolerance = 0.01;

    m = 1;
    n = nele;

    dfdx.resize(n, 0.0);
    hessf.resize(n, 0.0);
    gx.resize(m, 1.0);
    dgdx.resize(m, vector<double>(nele, 0));
    hessg.resize(m, vector<double>(nele, 0));
}
