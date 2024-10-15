#include <iostream>
#include <string>
#include <vector>
#include <utility>
#include <chrono>
#include <cmath>
#include <tuple>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Constants
const double XjL = -M_PI;
const double XjU = M_PI;
const double epsilon = 1e-6;

struct Solution {
    std::vector<double> xThetas;
    std::vector<double> xBetas;
    double energy;
    double F; // ------------------------------- not sure if I should work with float or double values, float -> faster but less precise, double -> slower but more precise ----------------------------------
    double Cr;
};

double energyCalculation(std::string aminoacids, Solution element) {
    return energyOne(element.xThetas) + energyTwo(aminoacids, element);
}
double energyOne(std::vector<double> elementThetas) {
    double result = 0.0;
    for (int i = 0; i < elementThetas.size(); i++) {
        result += (1 - cos(elementThetas[i]));
    }
    return (1 / 4) * result;
}
double energyTwo(std::string aminoacids, Solution element) {
    double result = 0.0;
    for (int i = 0; i < element.xThetas.size(); i++) {
        for (int j = i + 2; j < element.xThetas.size() + 2; j++) {
            result = 
                powerOf(1 / distance(coordinate(element, i), coordinate(element, j)), 12)
                - c(aminoacids[i], aminoacids[j])
                * powerOf(1 / distance(coordinate(element, i), coordinate(element, j)), 6)
                ; // ------ in debug, check if I get a proper value -----------
        }
    }
    return (1 / 4) * result;
}

double distance(std::tuple<double, double, double> ri, std::tuple<double, double, double> rj) {
    return sqrt(
        (std::get<0>(rj) - std::get<0>(ri)) * (std::get<0>(rj) - std::get<0>(ri)) +
        (std::get<1>(rj) - std::get<1>(ri)) * (std::get<1>(rj) - std::get<1>(ri)) +
        (std::get<2>(rj) - std::get<2>(ri)) * (std::get<2>(rj) - std::get<2>(ri))
    );
}

std::tuple<double, double, double> coordinate(Solution element, int i) { // ----------------------- if there is time, optimise so that Solution doesen't get copied as much --------------------------------------
    if (i == 1) return std::make_tuple(0.0, 0.0, 0.0);
    else if (i == 2) return std::make_tuple(0.0, 1.0, 0.0);
    else if (i == 3) return std::make_tuple(cos(element.xThetas[0]), 1.0 + sin(element.xThetas[0]), 0.0);
    else {
        return std::make_tuple(
            std::get<0>(coordinate(element, i - 1)) + cos(element.xThetas[i - 2]) * cos(element.xBetas[i - 3]),
            std::get<1>(coordinate(element, i - 1)) + sin(element.xThetas[i - 2]) * cos(element.xBetas[i - 3]),
            std::get<2>(coordinate(element, i - 1)) + sin(element.xBetas[i - 3])
        );
    }
}

double c(char aminoacidI, char aminoacidJ) {
    if (aminoacidI == 'A' && aminoacidJ == 'A') return 1.0;
    if (aminoacidI == 'B' && aminoacidJ == 'B') return 0.5;
    if (aminoacidI != aminoacidJ) return -0.5;

    throw std::runtime_error("Invalid amino acid pair: " + std::string(1, aminoacidI) + " and " + std::string(1, aminoacidJ));
}



double powerOf(double num, int pow) {
    double result = 1.0;
    for (int i = 0; i < pow; i++) {
        result *= num;
    }
    return result;
}

int returnRandomIndex(int size, int i) {
    int returnIndex = i;
    while (i == returnIndex) returnIndex = rand() % size;
    return returnIndex;
}

int returnRandomIndex(int size, int i, int r1) {
    int returnIndex = i;
    while (i == returnIndex || i == r1) returnIndex = rand() % size;
    return returnIndex;
}

int returnRandomIndex(int size, int i, int r1, int r2) {
    int returnIndex = i;
    while (i == returnIndex || i == r1 || i == r2) returnIndex = rand() % size;
    return returnIndex;
}


float randValBetween0And1() {
    return (rand() % 11) * 0.1; // -------------------------------------------------------------- High probability that these are not done correctly and should be of higher decimal values --------------------------------------------------------------

}


int main(int argc, char* argv[]) {
    std::string aminoacids = "";
    unsigned int seed = 0;
    double target = 0.0;
    unsigned int nfesLmt = 0;
    unsigned int nfesCounter = 0;
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
    std::cout << "\tSeed: " << seed << std::endl; // Seed for the rand() function
    std::cout << "\tTarget: " << target << std::endl; // Target energy
    std::cout << "\tNfesLmt: " << nfesLmt << std::endl; // How many times the function for calculating energy can be called
    std::cout << "\tRuntimeLmt (ms): " << runtimeLmt << std::endl; // How much time the algorithm can run before it shuts down (900000 -> 15min)
    std::cout << "\tNp: " << np << std::endl;

    // Example usage of seed (set the random seed)
    srand(seed);

    // Initialisation
    std::vector<Solution> population;
    for (int i = 0; i < np; i++) {
        Solution sol;
        sol.F = 0.5;
        sol.Cr = 0.9;

        sol.xThetas.push_back(rand() * (XjU - XjL) + XjL);
        for (int j = 3; j < aminoacids.size(); j++) {
            sol.xThetas.push_back(rand() * (XjU - XjL) + XjL);
            sol.xBetas.push_back(rand() * (XjU - XjL) + XjL);
        }
        population.push_back(sol);
    }

    // Initialisation is properly done (probably, maybe there is a few semantical problems)

    // std::cout << population[3].xBetas.size() << "\n" << population[3].xThetas.size();

    auto start = std::chrono::high_resolution_clock::now();
    auto now = std::chrono::high_resolution_clock::now();
    auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(now - start).count();
    bool finished = true;
    while (elapsed_ms >= runtimeLmt || nfesCounter >= nfesLmt || finished) {
        Solution u = Solution();
        now = std::chrono::high_resolution_clock::now();
        elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(now - start).count();
        for (int i = 0; i < np; i++) {
            int randomValue = randValBetween0And1();
            int F = 0.5;
            int Cr = 0.9;

            // Mutations
            if (randomValue < 0.1) population[i].F = 0.1 + 0.9 * randomValue;
            else population[i].F = F;
            
            randomValue = randValBetween0And1();
            if (randomValue < 0.1) population[i].Cr = randomValue;
            else population[i].Cr = Cr;


            // indexes are done properly
            int r1 = returnRandomIndex(np, i);
            int r2 = returnRandomIndex(np, i, r1);
            int r3 = returnRandomIndex(np, i, r1, r2); // -------------------------------------------------------------- možno, da bi moral tu vzeti best namesto r3--------------------------------------------------------------


            int jRand = rand() % aminoacids.size(); // -------------------------------------------------------------- možno, da D ni pravilen --------------------------------------------------------------
        
            // first go through for all thetas
            for (int j = 0; j < aminoacids.size() - 2; i++) { // -------------------------------------------------------------- možno, da D ni pravilen --------------------------------------------------------------
                if (randValBetween0And1() < Cr || j == jRand) {
                    u.xThetas[j] = population[r3].xThetas[j] + F * (population[r1].xThetas[j] - population[r2].xThetas[j]);
                    if (u.xThetas[j] <= -M_PI) u.xThetas[j] = 2 * M_PI + u.xThetas[j];
                    if (u.xThetas[j] > M_PI) u.xThetas[j] = 2 * (-M_PI) + u.xThetas[j];
                }
                else u.xThetas[j] = population[i].xThetas[j];
            }
            // then go through for all betas
            for (int j = 0; j < aminoacids.size() - 3; i++) { // -------------------------------------------------------------- možno, da D ni pravilen --------------------------------------------------------------
                if (randValBetween0And1() < Cr || j == jRand) {
                    u.xBetas[j] = population[r3].xBetas[j] + F * (population[r1].xBetas[j] - population[r2].xBetas[j]);
                    if (u.xBetas[j] <= -M_PI) u.xBetas[j] = 2 * M_PI + u.xBetas[j];
                    if (u.xBetas[j] > M_PI) u.xBetas[j] = 2 * (-M_PI) + u.xBetas[j];
                }
                else u.xBetas[j] = population[i].xBetas[j];
            }

            u.energy = energyCalculation(aminoacids, u);
            std::cout << "PLKS WORK BRO: " << u.energy << " ms" << std::endl;
        }
    }
    auto end = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Time taken: " << duration << " ms" << std::endl;

    return 0;
}