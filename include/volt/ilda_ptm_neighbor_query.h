#pragma once

#include <volt/analysis/ptm.h>
#include <volt/math/vector3.h>
#include <volt/math/quaternion.h>
#include <volt/structures/crystal_structure_types.h>

#include <cstddef>
#include <cstdint>
#include <vector>

namespace Volt {

class IldaPtmNeighborQuery {
public:
    struct Neighbor {
        int index = -1;
        Vector3 delta{0.0, 0.0, 0.0};
        Vector3 idealVector{0.0, 0.0, 0.0};
    };

    struct AtomState {
        StructureType structureType = OTHER;
        Quaternion orientation{0.0, 0.0, 0.0, 1.0};
        std::uint64_t correspondencesCode = 0;
        double interatomicDistance = 0.0;
        std::uint8_t valid = 0;
    };

    explicit IldaPtmNeighborQuery(PTM& ptm);

    void computeAtomStates();

    void setAtomStates(std::vector<AtomState> states);

    const std::vector<AtomState>& atomStates() const { return _states; }

    void findNeighbors(std::size_t particleIndex, const Quaternion* qref);

    int count() const { return static_cast<int>(_neighbors.size()); }
    StructureType structureType() const { return _structureType; }
    Quaternion orientation() const { return _orientation; }

    const Neighbor& operator[](int slot) const { return _neighbors[static_cast<std::size_t>(slot)]; }

private:
    PTM& _ptm;
    std::vector<AtomState> _states;
    std::vector<Neighbor> _neighbors;
    StructureType _structureType = OTHER;
    Quaternion _orientation{0.0, 0.0, 0.0, 1.0};
};

double scalePTMVec(StructureType type, double a);

double scalePTMInteratomicDistance(StructureType type, double a);

}
