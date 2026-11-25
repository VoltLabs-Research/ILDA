#pragma once

#include <volt/ilda_types.h>
#include <volt/ilda_circuit.h>
#include <volt/core/simulation_cell.h>
#include <volt/math/vector3.h>

#include <utility>
#include <vector>

namespace Volt {

class IldaLog;

extern const Vector3 ILDA_COLORLIST[20];
inline constexpr int ILDA_NCOLORS = 20;

// The histogram and segment helpers of the defect search.
namespace IldaDefects {

double vecAngle(const Vector3& v1, const Vector3& v2);

std::pair<bool, double> parallelTest(const Vector3& v1, const Vector3& v2, double angtol);

std::vector<double> arange(double start, double stop, double step);

std::vector<double> binMids(const std::vector<double>& edges);

int findBin(double v, const std::vector<double>& edges);

int argminAbsDiff(const std::vector<double>& mids, double val);

double interp01(double x, const std::vector<double>& fp);

int getId2(double bnorm, double hA, double hB, double btol,
           const std::vector<double>& bbinmids,
           const std::vector<double>& hbinmidsA,
           const std::vector<double>& hbinmidsB,
           const std::vector<int>& bhhistIds,
           int nHA, int nHB);

enum class Edge1Kind { NoAtom, UsePcount, Real };

void setMidedge(IldaMesh::Edge* e, int value);

std::pair<int, int> segmentAtom(Edge1Kind kind, IldaMesh::Edge* edge1,
                                IldaMesh::Edge* edge2, int pcount);

}

void findDefects(
    IldaMesh& mesh,
    const CircuitContext& ctx,
    double btol,
    double htol,
    double angtol,
    double bmagmax,
    std::vector<IldaSegment>& outSegments,
    std::vector<DisconnectionMode>& outModes,
    double& outTotalLength,
    IldaLog& log
);

}
