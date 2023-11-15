#include <iostream>
#include <dirent.h>
#include <string>
#include "topopt.h"

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

int main()
{
    using namespace std;

    string directoryPath = "data";
    deletePreviousImages(directoryPath);
    
    int nelx, nely;
    double volfrac, penal, rmin, f;
    
    cout << "Enter 'nelx', 'nely', 'volfrac', 'penal', 'rmin' in sequence" << endl;
    cout << "(Example: 30 10 0.5 3.0 1.5)" << endl;
    cin >> nelx >> nely >> volfrac >> penal >> rmin;

    // nelx = 30;
    // nely = 10;
    // volfrac = 0.5;
    // penal = 3.0;
    // rmin = 1.5;
    f = 1.0;

    topopt *top = new topopt(nelx, nely, volfrac, penal, rmin, f);
    top->sigmundopt();
    
    return 0;
}
