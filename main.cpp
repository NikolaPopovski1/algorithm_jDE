#include <iostream>
#include <string>
#include <vector>
#include <utility>
#include <chrono>
#include <cmath>
#include <tuple>
#include <cfloat>

#ifndef M_PI
#define M_PI 3.1415926f
#endif

// Constants
constexpr const float XjL = -(float)M_PI;
constexpr const float XjU = (float)M_PI;
constexpr const float epsilon = (float)1e-6;

struct Solution {
    std::vector<float> xThetas;
    std::vector<float> xBetas;
    float energy;
    float F;
    float Cr;
};

struct Coordinates {
    float x, y, z;

    Coordinates(float xValue = 0.0f, float yValue = 0.0f, float zValue = 0.0f) : x(xValue), y(yValue), z(zValue) {}
};

inline Coordinates coordinate(const std::vector<Coordinates>& allCoords, const Solution& element, int i) {
    if (i == 0) return Coordinates(0.0f, 0.0f, 0.0f);
    if (i == 1) return Coordinates(0.0f, 1.0f, 0.0f);
    if (i == 2) return Coordinates(cosf(element.xThetas[0]), 1.0f + sinf(element.xThetas[0]), 0.0f);
    if (3 <= i && i <= element.xThetas.size() + 2) {
        int idx = i - 3;
        return Coordinates(
                allCoords[i - 1].x + cosf(element.xThetas[idx]) * cosf(element.xBetas[idx]),
                allCoords[i - 1].y + sinf(element.xThetas[idx]) * cosf(element.xBetas[idx]),
                allCoords[i - 1].z + sinf(element.xBetas[idx])
        );
    }
    throw std::runtime_error("Invalid index in coordinate calculation");
}

inline float distance(const Coordinates& ri, const Coordinates& rj) {
    float dx = rj.x - ri.x;
    float dy = rj.y - ri.y;
    float dz = rj.z - ri.z;
    return sqrtf(dx * dx + dy * dy + dz * dz);
}

constexpr float c(char aminoacidI, char aminoacidJ) {
    if (aminoacidI == 'A' && aminoacidJ == 'A') return 1.0f;
    if (aminoacidI == 'B' && aminoacidJ == 'B') return 0.5f;
    return -0.5f;
}

float energyCalculation(const std::string& aminoacids, const Solution* element) {
    float resultOne = 0.0f;
    for (float theta : element->xThetas) {
        resultOne += (1.0f - cosf(theta));
    }
    resultOne /= 4.0f;

    std::vector<Coordinates> allCoords(element->xThetas.size() + 2);
    for (int i = 0; i < allCoords.size(); ++i) {
        allCoords[i] = coordinate(allCoords, *element, i);
    }

    float resultTwo = 0.0f;
    for (int i = 0; i < element->xThetas.size(); i++) {
        for (int j = i + 2; j < element->xThetas.size() + 2; j++) {
            float inv = distance(allCoords[i], allCoords[j]);
            float inv2 = inv * inv;
            float inv3 = inv * inv2;
            float inv6 = 1.0f / (inv3 * inv3);
            float inv12 = inv6 * inv6;
            resultTwo += (inv12 - c(aminoacids[i], aminoacids[j]) * inv6);
        }
    }

    return resultOne + resultTwo * 4.0f;
}

int main(int argc, char* argv[]) {
    std::string aminoacids;
    unsigned int seed = 0, nfesLmt = 0, nfesCounter = 0, runtimeLmt = 0, np = 0;
    float target = 0.0f;

    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if ((arg[0] == 'A' || arg[0] == 'B')) aminoacids = argv[i];
        else if (arg == "-seed") seed = std::stoi(argv[++i]);
        else if (arg == "-target") target = std::stof(argv[++i]);
        else if (arg == "-nfesLmt") nfesLmt = std::stoi(argv[++i]);
        else if (arg == "-runtimeLmt") runtimeLmt = std::stoi(argv[++i]);
        else if (arg == "-Np" || arg == "-np") np = std::stoi(argv[++i]);
    }

    int bestEnergy = 0;
    float F = 0.5f;
    float Cr = 0.9f;
    srand(seed);

    auto start = std::chrono::high_resolution_clock::now();
    auto now = std::chrono::high_resolution_clock::now();
    auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(now - start).count();

    // Initialisation
    Solution* populationCurrGen[np];
    for (unsigned int i = 0; i < np; i++) {
        Solution* sol = new Solution();
        sol->F = 0.5f;
        sol->Cr = 0.9f;
        sol->xThetas.reserve(aminoacids.size());
        sol->xBetas.reserve(aminoacids.size() - 3);

        for (int j = 0; j < aminoacids.size() - 2; j++) {
            sol->xThetas.push_back(XjL + 2 * XjU * ((float)rand() / RAND_MAX));
        }
        for (int j = 0; j < aminoacids.size() - 3; j++) {
            sol->xBetas.push_back(XjL + 2 * XjU * ((float)rand() / RAND_MAX));
        }

        sol->energy = energyCalculation(aminoacids, sol);
        populationCurrGen[i] = sol;
        nfesCounter++;
    }

    Solution u, uNew;
    int r1, r2, r3, rBest, jRand;
    float tmp;
    while (elapsed_ms <= runtimeLmt && nfesCounter <= nfesLmt) {
        now = std::chrono::high_resolution_clock::now();
        elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(now - start).count();

        for (unsigned int i = 0; i < np; i++) {
            u = Solution();

            // Mutation
            float randomValue = (float)rand() / (RAND_MAX + 1) + (rand() % 1);
            if (randomValue < 0.1) F = 0.1f + 0.9f * randomValue;
            else F = populationCurrGen[i]->F;

            // Crossing
            randomValue = (float)rand() / (RAND_MAX + 1) + (rand() % 1);
            if (randomValue < 0.1) Cr = randomValue;
            else Cr = populationCurrGen[i]->Cr;


            // Part of mutation
            // indexes are done properly
            do {
                r1 = rand() % np;
            } while (r1 == i);
            do {
                r2 = rand() % np;
            } while (r2 == i || r2 == r1);
            do {
                r3 = rand() % np;
            } while (r3 == i || r3 == r1 || r3 == r2);

            // Selection
            // & part of crossing

            jRand = rand() % aminoacids.size(); // -------------------------------------------------------------- možno, da D ni pravilen --------------------------------------------------------------
            // first go through for all thetas
            for (int j = 0; j < aminoacids.size() - 2; j++) { // -------------------------------------------------------------- možno, da D ni pravilen, mozno tud da implementacija ni pravilna in da bi moral implementirat  --------------------------------------------------------------
                if (((float)rand() / (RAND_MAX + 1) + (rand() % 1)) < Cr || j == jRand) {
                    u.xThetas.push_back(populationCurrGen[r3]->xThetas[j] + F * (populationCurrGen[r1]->xThetas[j] - populationCurrGen[r2]->xThetas[j])); // Ali je tu xb,j v psevdokodi mišljen x best?
                    if (u.xThetas[j] <= XjL) u.xThetas[j] = 2 * XjU + u.xThetas[j];
                    if (u.xThetas[j] > XjU) u.xThetas[j] = 2 * XjL + u.xThetas[j];
                }
                else u.xThetas.push_back(populationCurrGen[i]->xThetas[j]);
            }
            // then go through for all betas
            for (int j = 0; j < aminoacids.size() - 3; j++) { // -------------------------------------------------------------- možno, da D ni pravilen --------------------------------------------------------------
                if (((float)rand() / (RAND_MAX + 1) + (rand() % 1)) < Cr || j == jRand) {
                    u.xBetas.push_back(populationCurrGen[r3]->xBetas[j] + F * (populationCurrGen[r1]->xBetas[j] - populationCurrGen[r2]->xBetas[j]));
                    if (u.xBetas[j] <= XjL) u.xBetas[j] = 2 * XjU + u.xBetas[j];
                    if (u.xBetas[j] > XjU) u.xBetas[j] = 2 * XjL + u.xBetas[j];
                }
                else u.xBetas.push_back(populationCurrGen[i]->xBetas[j]);
            }
            u.energy = energyCalculation(aminoacids, &u);
            nfesCounter++;

            if (u.energy <= populationCurrGen[i]->energy) {
                uNew = Solution();
                for (int j = 0; j < aminoacids.size() - 2; j++) { // -------------------------------------------------------------- možno, da D ni pravilen, mozno tud da implementacija ni pravilna in da bi moral implementirat  --------------------------------------------------------------
                    uNew.xThetas.push_back(populationCurrGen[r3]->xThetas[j] + 0.5f * (u.xThetas[j] - populationCurrGen[i]->xThetas[j]));
                    if (uNew.xThetas[j] <= XjL) uNew.xThetas[j] = 2 * XjU + uNew.xThetas[j];
                    if (uNew.xThetas[j] > XjU) uNew.xThetas[j] = 2 * XjL + uNew.xThetas[j];
                }
                // then go through for all betas
                for (int j = 0; j < aminoacids.size() - 3; j++) { // -------------------------------------------------------------- možno, da D ni pravilen --------------------------------------------------------------
                    uNew.xBetas.push_back(populationCurrGen[r3]->xBetas[j] + 0.5f * (u.xBetas[j] - populationCurrGen[i]->xBetas[j]));
                    if (uNew.xBetas[j] <= XjL) uNew.xBetas[j] = 2 * XjU + uNew.xBetas[j];
                    if (uNew.xBetas[j] > XjU) uNew.xBetas[j] = 2 * XjL + uNew.xBetas[j];
                }

                uNew.energy = energyCalculation(aminoacids, &uNew);
                nfesCounter++;

                if (uNew.energy <= u.energy) {
                    populationCurrGen[i]->xThetas = uNew.xThetas;
                    populationCurrGen[i]->xBetas = uNew.xBetas;
                    populationCurrGen[i]->energy = uNew.energy;
                    populationCurrGen[i]->Cr = Cr;
                    populationCurrGen[i]->F = 0.5f;
                } else {
                    populationCurrGen[i]->xThetas = u.xThetas;
                    populationCurrGen[i]->xBetas = u.xBetas;
                    populationCurrGen[i]->energy = u.energy;
                    populationCurrGen[i]->Cr = Cr;
                    populationCurrGen[i]->F = F;
                }
            }
        }
    }
    auto end = std::chrono::high_resolution_clock::now();

    rBest = 0;
    tmp = populationCurrGen[0]->energy;
    for (int j = 1; j < np; j++) {
        if (tmp > populationCurrGen[j]->energy) {
            rBest = j;
            tmp = populationCurrGen[j]->energy;
        }
    }
    bestEnergy = rBest;

    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << populationCurrGen[bestEnergy]->energy << '\n';
    std::cout << nfesCounter << '\n';
    std::cout << duration << '\n';
    //std::cout << "Speed (nfes / duration = speed): " << nfesCounter << " / " << duration << " = " << nfesCounter/duration << '\n';

    for (int i = 0; i < np; i++) {
        delete populationCurrGen[i];
    }

    return 0;
}