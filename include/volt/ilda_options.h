#pragma once

#include <volt/math/vector3.h>
#include <volt/math/matrix3.h>
#include <string>

namespace Volt {

struct IldaOptions {
    bool selection_only = false;
    bool print_results = false;
    bool single_circuit = false;
    bool extract_lines = false;
    bool estimateF = true;

    int atomA = 0;
    int atomB = 0;
    int circuitAtom1 = 0;
    int circuitAtom2 = 0;

    double aA = 1.0;
    double cA = 0.0;
    double aB = 1.0;
    double cB = 0.0;

    int typeA = -1;                  // -1 => derive from grains table
    int typeB = -1;                  // -1 => derive from grains table

    double cis_tol = 0.0;
    double Rsphere = 10.0;
    double htol = 0.5;
    double btol = 0.01;
    double angtol = 5.0;             // degrees (converted to radians at use)
    double distF = 10.0;

    Vector3 n{0.0, 0.0, 0.0};

    Vector3 xA{0.0, 0.0, 0.0};
    Vector3 yA{0.0, 0.0, 0.0};
    Vector3 xB{0.0, 0.0, 0.0};
    Vector3 yB{0.0, 0.0, 0.0};

    Matrix3 EcohA{Matrix3::Zero()};
    Matrix3 EcohB{Matrix3::Zero()};

    std::string grainAtomsPath;
    std::string grainsPath;

    std::string crystalStructure = "FCC";
    double rmsdCutoff = 0.1;
};

}
