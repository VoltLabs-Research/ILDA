#include <volt/ilda_ptm_neighbor_query.h>

#include <ptm_constants.h>
#include <ptm_functions.h>

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>

#include <array>
#include <cmath>
#include <cstring>
#include <spdlog/spdlog.h>

namespace Volt {

IldaPtmNeighborQuery::IldaPtmNeighborQuery(PTM& ptm) : _ptm(ptm) {}

void IldaPtmNeighborQuery::setAtomStates(std::vector<AtomState> states) {
    _states = std::move(states);
}

void IldaPtmNeighborQuery::computeAtomStates() {
    const std::size_t n = _ptm.particleCount();
    _states.assign(n, AtomState{});
    if(n == 0) {
        return;
    }

    std::vector<std::uint64_t> cached(n, 0ull);

    tbb::parallel_for(tbb::blocked_range<std::size_t>(0, n), [&](const tbb::blocked_range<std::size_t>& range) {
        PTM::Kernel kernel(_ptm);
        for(std::size_t i = range.begin(); i < range.end(); ++i) {
            kernel.cacheNeighbors(i, &cached[i]);
            const StructureType type = kernel.identifyStructure(i, cached);
            AtomState& state = _states[i];
            state.structureType = type;
            if(type == OTHER) {
                continue;
            }
            state.orientation = kernel.orientation().normalized();
            state.correspondencesCode = kernel.correspondencesCode();
            state.interatomicDistance = kernel.interatomicDistance();
            state.valid = 1;
        }
    });
}

void IldaPtmNeighborQuery::findNeighbors(std::size_t particleIndex, const Quaternion* qref) {
    _neighbors.clear();

    if(particleIndex >= _states.size()) {
        _structureType = OTHER;
        _orientation = Quaternion(0.0, 0.0, 0.0, 1.0);
        return;
    }

    const AtomState& state = _states[particleIndex];
    _structureType = state.structureType;
    _orientation = state.orientation;

    if(!state.valid || state.structureType == OTHER) {
        _structureType = OTHER;
        return;
    }

    const int ptmType = PTM::toPtmStructureType(state.structureType);
    if(ptmType == PTM_MATCH_NONE) {
        _structureType = OTHER;
        return;
    }
    const int numNbrs = ptm_num_nbrs[ptmType];

    std::array<std::int8_t, PTM_MAX_INPUT_POINTS> mapping{};
    int templateIndex = 0;
    ptm_decode_correspondences(
        ptmType,
        state.correspondencesCode,
        mapping.data(),
        &templateIndex
    );

    // PARITY: PTM C API quaternion order is {w,x,y,z}; Volt Quaternion is {x,y,z,w}.
    double q[4] = {
        state.orientation.w(),
        state.orientation.x(),
        state.orientation.y(),
        state.orientation.z()
    };

    if(qref != nullptr) {
        double qtarget[4] = { qref->w(), qref->x(), qref->y(), qref->z() };
        const int newTemplateIndex = ptm_remap_template(
            ptmType,
            templateIndex,
            qtarget,
            q,
            mapping.data()
        );
        if(newTemplateIndex < 0) {
            spdlog::warn("IldaPtmNeighborQuery: ptm_remap_template failed for atom {}", particleIndex);
            _structureType = OTHER;
            return;
        }
        templateIndex = newTemplateIndex;
        _orientation = Quaternion(q[1], q[2], q[3], q[0]);
    }

    const double (*templatePoints)[3] = PTM::getTemplate(state.structureType, templateIndex);
    if(templatePoints == nullptr) {
        _structureType = OTHER;
        return;
    }

    _neighbors.reserve(static_cast<std::size_t>(numNbrs));
    for(int k = 0; k < numNbrs; ++k) {
        // PARITY: correspondences/template are 1-based ([0]=central), so offset +1.
        const int nbrListPos = mapping[static_cast<std::size_t>(k + 1)] - 1;
        const int neighborAtomIndex = _ptm.cachedNeighborIndex(particleIndex, nbrListPos);

        Neighbor neighbor;
        neighbor.index = neighborAtomIndex;
        neighbor.delta = (neighborAtomIndex >= 0)
            ? _ptm.cachedNeighborDelta(particleIndex, nbrListPos)
            : Vector3::Zero();
        neighbor.idealVector = Vector3(
            templatePoints[k + 1][0],
            templatePoints[k + 1][1],
            templatePoints[k + 1][2]
        );
        _neighbors.push_back(neighbor);
    }
}

double scalePTMVec(StructureType type, double a) {
    switch(type) {
        case FCC: return a / std::sqrt(2.0);
        case HCP: return a;
        case BCC: return a / (14.0 * std::sqrt(3.0) / 3.0 - 7.0);
        case CUBIC_DIAMOND: return a * std::sqrt(3.0 / 16.0) / 0.678083388794149;
        default:
            spdlog::warn("scalePTMVec: untested structure type {}", static_cast<int>(type));
            return a;
    }
}

double scalePTMInteratomicDistance(StructureType type, double a) {
    switch(type) {
        case FCC: return a * std::sqrt(2.0);
        case HCP: return a * std::sqrt(2.0);
        case BCC: return a * std::sqrt(4.0 / 3.0);
        case CUBIC_DIAMOND: return a * std::sqrt(16.0 / 3.0);
        default: return a;
    }
}

}
