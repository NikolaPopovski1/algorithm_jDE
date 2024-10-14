#include <iostream>
#include <string>
#include <cstdlib>
#include <vector>
#include <utility>

// Constants
const double PI = 3.141592653589793;
const double XjL = -PI;
const double XjU = PI;
const double epsilon = 1e-6;

struct Solution {
    std::vector<double> xBetas;
    std::vector<double> xFis;
    double energy;
    double F;
    double Cr;
};

double energy1Calculation(double fi) {

}



int main(int argc, char* argv[]) {
    std::string aminoacids = "";
    unsigned int seed = 0;
    double target = 0.0;
    unsigned int nfesLmt = 0;
    unsigned int runtimeLmt = 0;
    unsigned int np = 0;

    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];

        if ((arg[0] == 'A' || arg[0] == 'B')) {
            aminoacids = argv[i];
        }
        else if (arg == "-seed" && i + 1 < argc) {
            seed = std::stoi(argv[++i]);
        }
        else if (arg == "-target" && i + 1 < argc) {
            target = std::stod(argv[++i]);
        }
        else if (arg == "-nfesLmt" && i + 1 < argc) {
            nfesLmt = std::stoi(argv[++i]);
        }
        else if (arg == "-runtimeLmt" && i + 1 < argc) {
            runtimeLmt = std::stoi(argv[++i]);
        }
        else if (arg == "-Np" && i + 1 < argc) {
            np = std::stoi(argv[++i]);
        }
    }

    // Print parsed values to check
    std::cout << "Input values:" << std::endl;
    std::cout << "\tAminoacids: " << aminoacids << std::endl;
    std::cout << "\tSeed: " << seed << std::endl;
    std::cout << "\tTarget: " << target << std::endl;
    std::cout << "\tNfesLmt: " << nfesLmt << std::endl;
    std::cout << "\tRuntimeLmt: " << runtimeLmt << std::endl;
    std::cout << "\tNp: " << np << std::endl;

    // Example usage of seed (set the random seed)
    srand(seed);

    // Initialisation
    std::vector<Solution> population;
    for (int i = 0; i < np; i++) {
        Solution sol;
        sol.F = 0.5;
        sol.Cr = 0.9;
        sol.xFis.push_back(rand() * (XjU - XjL) + XjL);

        

        for (int j = 3; j < aminoacids.size(); j++) {
            sol.xBetas.push_back(rand() * (XjU - XjL) + XjL);
            sol.xFis.push_back(rand() * (XjU - XjL) + XjL);
        }
        population.push_back(sol);
    }

    // 


    return 0;
}

