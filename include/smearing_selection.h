#ifndef _SelectionType_h_
#define _SelectionType_h_
#include <string>
namespace erAna {
  enum ESelection {
    // default selection cuts, data distribution
    eDefault = 0,
    // default selection cuts, pure proton sims
    eDefaultProton,
    // default selection cuts, pure iron sims
    eDefaultIron,
    // default selection cuts, 50:50 sims
    eDefaultMixed,
    // remove FidFOVICRC13prel
    eNoFOV1,
    // remove FidFOVICRC13prel, E_sys up
    eNoFOV1Up,
    // remove FidFOVICRC13prel, E_sys low
    eNoFOV1Low,
    // remove FidFOVICRC13prel,
    // replace xMaxObsInExpectedFOV with xMaxInFOV/minViewAngle
    eNoFOV2,
    // default and minCosZFDE( 0.795, -0.092)
    eVertical,
    // default and !minCosZFDE( 0.795, -0.092)
    eInclined,
    eNSelections
  };

  inline std::string SelectionTypeName(const ESelection selType)
  {
    const std::string names[eNSelections] = {"default", "defaultP", "defaultFe", "defaultMix", "noFOV1", "noFOV1Up", "noFOV1low", "noFOV2", "vertical", "inclined"};
    return names[selType];
  }
}

#endif
