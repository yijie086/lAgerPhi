// clean_lund.cc
// Simple LUND cleaner: remove events that contain NaN/Inf in px, py, pz, or E.
//
// Usage:
//   g++ -std=c++17 -O2 clean_lund.cc -o clean_lund
//   ./clean_lund input.lund output_clean.lund
//
// This code assumes the LUND format you showed:
// Header (10 columns):
//   1: npart (int)
//   2: int
//   3: int
//   4: double
//   5: double
//   6: beam_pid (int)
//   7: beamE (double)
//   8: int
//   9: int
//   10: weight (double)
//
// Particle line (14 columns):
//   1: idx (int)
//   2: double
//   3: ist (int)
//   4: pid (int)
//   5: parent (int)
//   6: daughter (int)
//   7: px (double)
//   8: py (double)
//   9: pz (double)
//   10: E (double)
//   11: m (double)
//   12: vx (double)
//   13: vy (double)
//   14: vz (double)

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>

bool isBadDouble(double x) {
    // Return true if x is NaN or infinite
    return std::isnan(x) || !std::isfinite(x);
}

int main(int argc, char *argv[])
{
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0]
                  << " input.lund output_clean.lund" << std::endl;
        return 1;
    }

    const char *inFile  = argv[1];
    const char *outFile = argv[2];

    std::ifstream fin(inFile);
    if (!fin.is_open()) {
        std::cerr << "Error: cannot open input file: " << inFile << std::endl;
        return 1;
    }

    std::ofstream fout(outFile);
    if (!fout.is_open()) {
        std::cerr << "Error: cannot open output file: " << outFile << std::endl;
        return 1;
    }

    std::string line;
    long long nEvents      = 0;  // total header lines read
    long long nKeptEvents  = 0;  // events written to output
    long long nBadEvents   = 0;  // events skipped due to NaN/Inf

    while (true) {
        // --- Read header line ---
        if (!std::getline(fin, line)) {
            // EOF
            break;
        }

        // Skip empty or comment lines
        if (line.empty() || line[0] == '#') {
            fout << line << "\n"; // optional: keep comments in output
            continue;
        }

        std::istringstream hs(line);

        int    npart;
        int    h_i2, h_i3, h_beamPid, h_i8, h_i9;
        double h_d4, h_d5, beamE, h_weight;

        // Parse header line according to the format
        hs >> npart >> h_i2 >> h_i3
           >> h_d4 >> h_d5
           >> h_beamPid >> beamE
           >> h_i8 >> h_i9 >> h_weight;

        if (!hs) {
            std::cerr << "Warning: cannot parse header line: " << line << std::endl;
            // We don't know how many lines belong to this event, so just continue.
            continue;
        }

        ++nEvents;

        // If beam energy itself is NaN/Inf, mark event as bad.
        bool badEvent = false;
        if (isBadDouble(beamE)) {
            badEvent = true;
        }

        // We will store particle lines temporarily in a vector
        std::vector<std::string> particleLines;
        particleLines.reserve(npart);

        // --- Read particle lines for this event ---
        for (int i = 0; i < npart; ++i) {
            if (!std::getline(fin, line)) {
                std::cerr << "Error: unexpected EOF while reading event "
                          << nEvents << " particle " << i << std::endl;
                badEvent = true;
                break;
            }

            // Allow blank/comment lines inside event? Usually not, but we handle it.
            if (line.empty() || line[0] == '#') {
                // If a comment/empty line appears where a particle is expected,
                // we consider this event bad because the format is broken.
                badEvent = true;
                // Still push the line to keep the stream in sync
                particleLines.push_back(line);
                continue;
            }

            particleLines.push_back(line);

            // Parse particle line to check px, py, pz, E
            std::istringstream ps(line);

            int    idx, ist, pid, parent, daughter;
            double dummy2;
            double px, py, pz, E, m, vx, vy, vz;

            ps >> idx >> dummy2 >> ist >> pid >> parent >> daughter
               >> px >> py >> pz >> E >> m >> vx >> vy >> vz;

            if (!ps) {
                std::cerr << "Warning: cannot parse particle line at event "
                          << nEvents << " : " << line << std::endl;
                badEvent = true;
                continue;
            }

            if (isBadDouble(px) || isBadDouble(py) ||
                isBadDouble(pz) || isBadDouble(E)) {
                std::cerr << "Warning: NaN/Inf found in event "
                          << nEvents << " (px/py/pz/E). Event will be skipped."
                          << std::endl;
                badEvent = true;
                // continue reading remaining lines to keep stream in sync
            }
        }

        if (badEvent) {
            ++nBadEvents;
            // Do not write this event to output
            continue;
        }

        // --- Write clean event: header + particle lines ---
        fout << npart << " "
             << h_i2 << " " << h_i3 << " "
             << h_d4 << " " << h_d5 << " "
             << h_beamPid << " " << beamE << " "
             << h_i8 << " " << h_i9 << " " << h_weight << "\n";

        for (const auto &pl : particleLines) {
            fout << pl << "\n";
        }

        ++nKeptEvents;
    }

    std::cout << "Done." << std::endl;
    std::cout << "  Total events (headers) read : " << nEvents     << std::endl;
    std::cout << "  Events kept (written)      : " << nKeptEvents << std::endl;
    std::cout << "  Events skipped (NaN/Inf)   : " << nBadEvents  << std::endl;

    return 0;
}
