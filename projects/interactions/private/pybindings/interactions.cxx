
#include <vector>

#include "../../public/SIREN/interactions/Interaction.h"
#include "../../public/SIREN/interactions/CrossSection.h"
#include "../../public/SIREN/interactions/DarkNewsCrossSection.h"
#include "../../public/SIREN/interactions/DarkNewsDecay.h"
#include "../../public/SIREN/interactions/Decay.h"
#include "../../public/SIREN/interactions/DISFromSpline.h"
#include "../../public/SIREN/interactions/DummyCrossSection.h"
#include "../../public/SIREN/interactions/ElasticScattering.h"
#include "../../public/SIREN/interactions/ElectroweakDecay.h"
#include "../../public/SIREN/interactions/HNLDipoleDecay.h"
#include "../../public/SIREN/interactions/HNLDipoleFromTable.h"
#include "../../public/SIREN/interactions/HNLDipoleDISFromSpline.h"
#include "../../public/SIREN/interactions/HNLDISFromSpline.h"
#include "../../public/SIREN/interactions/HNLDecay.h"
#include "../../public/SIREN/interactions/InteractionCollection.h"
#include "../../public/SIREN/interactions/QuarkDISFromSpline.h"
#include "../../public/SIREN/interactions/CharmMesonDecay.h"
#include "../../public/SIREN/interactions/CharmMesonDecay3Body.h"
#include "../../public/SIREN/interactions/DMesonELoss.h"
#include "../../public/SIREN/interactions/PythiaDISCrossSection.h"


#include "./Interaction.h"
#include "./CrossSection.h"
#include "./DarkNewsCrossSection.h"
#include "./DarkNewsDecay.h"
#include "./DISFromSpline.h"
#include "./QuarkDISFromSpline.h"
#include "./Decay.h"
#include "./DummyCrossSection.h"
//#include "./ElasticScattering.h"
#include "./ElectroweakDecay.h"
#include "./HNLDipoleDecay.h"
#include "./HNLDipoleFromTable.h"
#include "./HNLDipoleDISFromSpline.h"
#include "./HNLDISFromSpline.h"
#include "./HNLDecay.h"
#include "./InteractionCollection.h"
#include "./CharmMesonDecay.h"
#include "./CharmMesonDecay3Body.h"
#include "./DMesonELoss.h"
#include "./PythiaDISCrossSection.h"
#include "./MarleyCrossSection.h"

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

PYBIND11_DECLARE_HOLDER_TYPE(T__,std::shared_ptr<T__>)

using namespace pybind11;

PYBIND11_MODULE(interactions,m) {
    using namespace siren::interactions;

    register_Interaction(m);
    register_CrossSection(m);
    register_Decay(m);
    register_DarkNewsCrossSection(m);
    register_DarkNewsDecay(m);
    register_DISFromSpline(m);
    register_DummyCrossSection(m);
    //register_ElasticScattering();
    register_ElectroweakDecay(m);
    register_HNLDipoleDecay(m);
    register_HNLDipoleFromTable(m);
    register_HNLDipoleDISFromSpline(m);
    register_HNLDISFromSpline(m);
    register_HNLDecay(m);
    register_InteractionCollection(m);
    register_QuarkDISFromSpline(m);
    register_CharmMesonDecay(m);
    register_CharmMesonDecay3Body(m);
    register_DMesonELoss(m);
    register_PythiaDISCrossSection(m);
    register_MarleyCrossSection(m);
}
