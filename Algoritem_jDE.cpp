#include <iostream>
#include <string>
#include <vector>
#include <utility>
#include <chrono>
#include <cmath>
#include <tuple>

#ifndef M_PI
#define M_PI 3.1415926f
#endif

float powerOf(float num, int pow);

// Constants
const float XjL = -(float)M_PI;
const float XjU = (float)M_PI;
const float epsilon = (float)1e-6;

struct Solution {
    std::vector<float> xThetas;
    std::vector<float> xBetas;
    float energy;
    float F; // ------------------------------- not sure if I should work with float or float values, float -> faster but less precise, double -> slower but more precise ----------------------------------
    float Cr;
};

struct Coordinates {
    float x, y, z;

    Coordinates() {}
    Coordinates(float xValue, float yValue, float zValue) : x(xValue), y(yValue), z(zValue) {}
};

Coordinates coordinate(Solution& element, int i) { // ----------------------- if there is time, optimise so that Solution doesen't get copied as much --------------------------------------
    if (i == 0) return Coordinates(0.0f, 0.0f, 0.0f);
    else if (i == 1) return Coordinates(0.0f, 1.0f, 0.0f);
    else if (i == 2) return Coordinates(cos(element.xThetas[0]), 1.0f + sin(element.xThetas[0]), 0.0f);
    else if (3 <= i <= element.xThetas.size() + 3) {
        return Coordinates(
            coordinate(element, i - 1).x + cos(element.xThetas[i - 2]) * cos(element.xBetas[i - 3]),
            coordinate(element, i - 1).y + sin(element.xThetas[i - 2]) * cos(element.xBetas[i - 3]),
            coordinate(element, i - 1).z + sin(element.xBetas[i - 3])
        );
    }
    else {
        throw std::runtime_error("Index for calculating coordinate not valid...");
        return Coordinates();
    }
}

float distance(Coordinates& ri, Coordinates& rj) {
    float first = rj.x - ri.x;
    float second = rj.y - ri.y;
    float third = rj.z - ri.z;
    return sqrt(first * first + second * second + third * third);
}

float c(char& aminoacidI, char& aminoacidJ) {
    if (aminoacidI == 'A' && aminoacidJ == 'A') return 1.0f;
    if (aminoacidI == 'B' && aminoacidJ == 'B') return 0.5f;
    if (aminoacidI != aminoacidJ) return -0.5f;

    throw std::runtime_error("Invalid amino acid pair: " + std::string(1, aminoacidI) + " and " + std::string(1, aminoacidJ));
}

float energyCalculation(std::string& aminoacids, Solution& element) {
    float resultOne = 0.0;
    for (int i = 0; i < element.xThetas.size(); i++) {
        resultOne += (1 - cos(element.xThetas[i]));
    }
    resultOne /= 4;

    float resultTwo = 1.0;
    for (int i = 0; i < element.xThetas.size(); i++) {
        for (int j = i + 2; j < element.xThetas.size() + 2; j++) {
            //std::cout << "-----" << element.xThetas.size() << "-----";
            Coordinates coordI = coordinate(element, i);
            Coordinates coordJ = coordinate(element, j);
            float inv = distance(coordI, coordJ);
            float inv2 = inv * inv;
            float inv3 = inv * inv2;
            float inv6 = 1 / (inv3 * inv3);
            float inv12 = inv6 * inv6;
            //std::cout << "i:" << i << ", j:" << j << ";";
            resultTwo += inv12 - c(aminoacids[i], aminoacids[j]) * inv6; // ------ in debug, check if I get a proper value -----------
        }
    }
    return resultTwo * 4;
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
    return (float)rand() / (RAND_MAX + 1) + (rand() % 1);; // -------------------------------------------------------------- High probability that these are not done correctly and should be of higher decimal values --------------------------------------------------------------
}


int main(int argc, char* argv[]) {
    std::string aminoacids;
    unsigned int seed;
    float target;
    unsigned int nfesLmt;
    unsigned int nfesCounter = 0;
    unsigned int runtimeLmt;
    unsigned int np;

    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];

        if ((arg[0] == 'A' || arg[0] == 'B')) {
            aminoacids = argv[i];
        }
        else if (arg == "-seed" && i + 1 < argc) {
            seed = std::stoi(argv[++i]);
        }
        else if (arg == "-target" && i + 1 < argc) {
            target = std::stof(argv[++i]);
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
    std::cout << "Input values:" << '\n';
    std::cout << "\tAminoacids: " << aminoacids << '\n';
    std::cout << "\tSeed: " << seed << '\n'; // Seed for the rand() function
    std::cout << "\tTarget: " << target << '\n'; // Target energy
    std::cout << "\tNfesLmt: " << nfesLmt << '\n'; // How many times the function for calculating energy can be called
    std::cout << "\tRuntimeLmt (ms): " << runtimeLmt << '\n'; // How much time the algorithm can run before it shuts down (900000 -> 15min)
    std::cout << "\tNp: " << np << '\n';

    // Example usage of seed (set the random seed)
    srand(seed);

    // Initialisation
    std::cout << "Init:";
    std::vector<Solution*> populationPrevGen;
    std::vector<Solution*> populationCurrGen;
    std::vector<Solution*> populationNextGen;
    for (unsigned int i = 0; i < np; i++) {
        Solution* sol = new Solution();
        sol->F = 0.5f;
        sol->Cr = 0.9f;

        sol->xThetas.push_back(randValBetween0And1() * (XjU - XjL) + XjL);
        for (int j = 3; j < aminoacids.size(); j++) {
            sol->xThetas.push_back(randValBetween0And1() * (XjU - XjL) + XjL);
            sol->xBetas.push_back(randValBetween0And1() * (XjU - XjL) + XjL);
        }
        sol->energy = energyCalculation(aminoacids, *sol);
        populationPrevGen.push_back(sol);
        populationCurrGen.push_back(sol);

        std::cout << ".";
    }

    std::cout << "\nInit finished.";

    // Initialisation is properly done (probably, maybe there is a few semantical problems)

    // std::cout << populationCurrGen[3].xBetas.size() << "\n" << populationCurrGen[3].xThetas.size();

    std::cout << "\nStarting algorithm:";

    auto start = std::chrono::high_resolution_clock::now();
    auto now = std::chrono::high_resolution_clock::now();
    auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(now - start).count();

    float bestEnergy = 0.0f;
    bool finished = true;
    float F = 0.5f;
    float Cr = 0.9f;
    while (elapsed_ms <= runtimeLmt || nfesCounter <= nfesLmt || finished) {
        now = std::chrono::high_resolution_clock::now();
        elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(now - start).count();
        std::cout << "\n";

        for (unsigned int i = 0; i < np; i++) {
            Solution* u = new Solution();
            float randomValue = randValBetween0And1();

            // Mutation
            if (randomValue < 0.1) populationCurrGen[i]->F = 0.1f + 0.9f * randomValue;
            else populationCurrGen[i]->F = F;

            // Crossing
            randomValue = randValBetween0And1();
            if (randomValue < 0.1) populationCurrGen[i]->Cr = randomValue;
            else populationCurrGen[i]->Cr = Cr;


            // Part of mutation
            // indexes are done properly
            int r1 = returnRandomIndex(np, i);
            int r2 = returnRandomIndex(np, i, r1);
            // int r3 = returnRandomIndex(np, i, r1, r2); // -------------------------------------------------------------- možno, da bi moral tu vzeti best namesto r3--------------------------------------------------------------
            int rBest = 0;
            float tmp = populationCurrGen[0]->energy;
            for (int i = 1; i < populationCurrGen.size(); i++) {
                if (tmp > populationCurrGen[i]->energy) {
                    rBest = i;
                    tmp = populationCurrGen[i]->energy;
                }
            }

            // Part of crossing
            int jRand = rand() % aminoacids.size(); // -------------------------------------------------------------- možno, da D ni pravilen --------------------------------------------------------------

            // Part of crossing
            // first go through for all thetas
            for (int j = 0; j < aminoacids.size() - 2; j++) { // -------------------------------------------------------------- možno, da D ni pravilen, mozno tud da implementacija ni pravilna in da bi moral implementirat  --------------------------------------------------------------
                if (randValBetween0And1() < Cr || j == jRand) {
                    u->xThetas.push_back(populationCurrGen[rBest]->xThetas[j] + F * (populationCurrGen[r1]->xThetas[j] - populationCurrGen[r2]->xThetas[j]));
                    if (u->xThetas[j] <= -M_PI) u->xThetas[j] = 2 * M_PI + u->xThetas[j];
                    if (u->xThetas[j] > M_PI) u->xThetas[j] = 2 * (-M_PI) + u->xThetas[j];
                }
                else u->xThetas.push_back(populationCurrGen[i]->xThetas[j]);
            }
            // then go through for all betas
            for (int j = 0; j < aminoacids.size() - 3; j++) { // -------------------------------------------------------------- možno, da D ni pravilen --------------------------------------------------------------
                if (randValBetween0And1() < Cr || j == jRand) {
                    u->xBetas.push_back(populationCurrGen[rBest]->xBetas[j] + F * (populationCurrGen[r1]->xBetas[j] - populationCurrGen[r2]->xBetas[j]));
                    if (u->xBetas[j] <= -M_PI) u->xBetas[j] = 2 * M_PI + u->xBetas[j];
                    if (u->xBetas[j] > M_PI) u->xBetas[j] = 2 * (-M_PI) + u->xBetas[j];
                }
                else u->xBetas.push_back(populationCurrGen[i]->xBetas[j]);
            }

            u->energy = energyCalculation(aminoacids, *u);
            nfesCounter++;

            // Selection
            Solution* newU = new Solution();
            if (u->energy < populationCurrGen[i]->energy) {
                newU->Cr = populationCurrGen[i]->Cr;
                newU->F = populationCurrGen[i]->F;
                newU->xThetas = u->xThetas;
                newU->xBetas = u->xBetas;
            }
            else {
                if (i == 0) {
                    newU->Cr = populationCurrGen[i]->Cr;
                    newU->F = populationCurrGen[i]->F;
                }
                else {
                    newU->Cr = populationPrevGen[i]->Cr;
                    newU->F = populationPrevGen[i]->F;
                }
                newU->xThetas = populationCurrGen[i]->xThetas;
                newU->xBetas = populationCurrGen[i]->xBetas;
            }
            populationNextGen.push_back(newU);

            std::cout << ".";
            if (target - 0.000001f < u->energy < target + 0.000001f) finished = false;
            delete u;
        }

        // Check if the target energy (+- epsilon) has been reached, if so then end the while loop

        populationPrevGen = populationCurrGen;
        populationCurrGen = populationNextGen;
        populationNextGen.clear();
    }
    auto end = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Time taken: " << duration << " ms" << '\n';
    std::cout << "Best energy: " << bestEnergy << '\n';

    for (int i = 0; i < populationCurrGen.size(); i++) {
        delete populationCurrGen[i];
    }
    for (int i = 0; i < populationNextGen.size(); i++) {
        delete populationCurrGen[i];
    }

    return 0;
}