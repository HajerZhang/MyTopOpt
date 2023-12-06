#include <iostream>
#include <dirent.h>
#include <string>
#include "Topopt.h"
#include "FiniteElement.h"
#include "filter.h"
#include "mma.h"

void deletePreviousImages(const std::string& directoryPath) {
    DIR* dir = opendir(directoryPath.c_str());
    if (dir) {
        struct dirent* entry;
        while ((entry = readdir(dir)) != nullptr) {
            std::string fileName = entry->d_name;
            // Delete files with names starting with 'topopt_' and ending with '.txt'.
            if (fileName.length() > 4 &&
                fileName.substr(0, 7) == "topopt_" &&
                fileName.substr(fileName.length() - 4) == ".txt") {
                std::string filePath = directoryPath + "/" + fileName;
                std::remove(filePath.c_str());
            }
        }
        closedir(dir);
    }
}

void topmain();
void output(int loop, vector<double> x, int nelx, int nely);

int main()
{

    string directoryPath = "data";
    deletePreviousImages(directoryPath);
    // int dim;
    // while(true) {
    //     cout << "Enter the dim" << endl;
    //     cout << "('2' for 2D, '3' for 3D or '-1' for exit)" << endl;
    //     cin >> dim;
    //     if (dim == 2) {
    //         topmain(dim);
    //         break;
    //     }
    //     else if (dim == 3) {
    //         topmain(dim);
    //         break;
    //     }
    //     else if (dim == -1) {
    //         return 0;
    //     }
    //     else {
    //         cout << "Invalid input. Try again." << endl;
    //     }
    // }
    topmain();

    return 0;
}

void topmain()
{   
    int nelx, nely;
    double volfrac, penal, rmin;
    int ft;

    cout << "Enter 'nelx', 'nely', 'volfrac', 'penal', 'rmin' and 'filter' in sequence" << endl;
    cout << "(Example: 3 3 0.5 3.0 1.5 1)" << endl;
    // cin >> nelx >> nely >> volfrac >> penal >> rmin >> ft;
    nelx = 10;
    nely = 10;
    volfrac = 0.5;
    penal = 3.0;
    rmin = 1.5;
    ft = 1;

    // topopt model intialization
    topopt* opt = new topopt(nelx, nely, volfrac, penal, rmin, ft);
    // finite element intialization
    FiniteElement fe(opt->E, opt->nu,  opt->nelx, opt->nely);
    fe.GetEleStiffMat();
    fe.BoundaryCondition();
    // filter intialization
    filter filt(opt->nelx, opt->nely, opt->rmin, opt->ft);
    // mma intialization
    mma *mmasolver = new mma(opt->n, opt->m, opt->x);

    int itr = 0;
    double ch = 1.0;

    while (itr < opt->maxiter && ch > opt->tolerance){
        itr++;

        fe.AssembleGloStiffMat(opt->xPhys, opt->Emax, opt->Emin, opt->penal); 
        fe.Solve();
        fe.GetUpdate(opt->fx, opt->gx, opt->dfdx, opt->dgdx,
            opt->hessf, opt->hessg,
            opt->xPhys, opt->Emin, opt->Emax, opt->penal, opt->volfrac, itr);
        
        filt.Filtering(opt->xPhys, opt->dfdx);
        filt.Filtering(opt->xPhys, opt->hessf);
        filt.Filtering(opt->x, opt->dgdx);
        filt.Filtering(opt->x, opt->hessg);

        // Call the Update function from the mma class
        ch = mmasolver->Update(itr, opt->x, opt->fx, opt->gx, opt->dfdx, opt->dgdx, opt->hessf, opt->hessg);
        cout << "iteration: " << itr << "  compliance: " << opt->fx << "  change: " << ch << endl;
        output(itr, opt->x, opt->nelx, opt->nely);
    }
}

void output(int loop, vector<double> x, int nelx, int nely)
{

    const char* folderPath = "./data";

    // Check if the "data" folder exists, create it if not
    struct stat info;
    if (stat(folderPath, &info) != 0 || !(info.st_mode & S_IFDIR)) {
        // Folder doesn't exist, create it
        if (mkdir(folderPath, 0777) != 0) {
            cerr << "Error creating directory: " << folderPath << std::endl;
            return;
        }
    }

    string filename = "./data/topopt_" + to_string(loop) + ".txt";
    ofstream outputFile(filename);

    if (!outputFile.is_open()) {
        cerr << "Error opening file: " << filename << std::endl;
        return;
    }
    int e=0;
    for (int i = 0; i < nely; ++i) {
        for (int j = 0; j < nelx; ++j) {
            outputFile << x[e] << " ";
            e++;
        }
        outputFile << endl; 
    }

    outputFile.close();
}