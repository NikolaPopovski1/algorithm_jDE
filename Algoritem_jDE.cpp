#include <iostream>
#include <string>
#include <cstdlib>

int main(int argc, char* argv[]) {
    unsigned int seed = 0;
    double target = 0.0;
    unsigned int nfesLmt = 0;
    unsigned int runtimeLmt = 0;
    unsigned int Np = 0;

    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];

        if ((arg == "A" || arg == "B") && i + 1 < argc) {
            seed = std::stoi(argv[++i]);
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
            Np = std::stoi(argv[++i]);
        }
    }

    // Initialisation


    // Print parsed values to check
    std::cout << "Seed: " << seed << std::endl;
    std::cout << "Target: " << target << std::endl;
    std::cout << "NfesLmt: " << nfesLmt << std::endl;
    std::cout << "RuntimeLmt: " << runtimeLmt << std::endl;
    std::cout << "Np: " << Np << std::endl;

    // Example usage of seed (set the random seed)
    srand(seed);

    std::cout << "Random number (multiplied by 33): " << rand() * 33 << std::endl;

    return 0;
}
