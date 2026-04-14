
#include <vector>

#include "../../public/SIREN/interactions/Interaction.h"
#include "../../public/SIREN/interactions/CrossSection.h"
#include "../../public/SIREN/interactions/NeutrissimoDecay.h"
#include "../../public/SIREN/interactions/InteractionCollection.h"
#include "../../public/SIREN/interactions/DISFromSpline.h"
#include "../../public/SIREN/interactions/CharmDISFromSpline.h"
#include "../../public/SIREN/interactions/QuarkDISFromSpline.h"
#include "../../public/SIREN/interactions/HNLFromSpline.h"
#include "../../public/SIREN/interactions/DipoleFromTable.h"
#include "../../public/SIREN/interactions/DarkNewsCrossSection.h"
#include "../../public/SIREN/interactions/DarkNewsDecay.h"
#include "../../public/SIREN/interactions/Hadronization.h"
#include "../../public/SIREN/interactions/CharmHadronization.h"
#include "../../public/SIREN/interactions/CharmMesonDecay.h"
#include "../../public/SIREN/interactions/CharmMesonDecay3Body.h"
#include "../../public/SIREN/interactions/DMesonELoss.h"
#include "../../public/SIREN/interactions/PythiaDISCrossSection.h"

#include "./Interaction.h"
#include "./CrossSection.h"
#include "./DipoleFromTable.h"
#include "./DarkNewsCrossSection.h"
#include "./DarkNewsDecay.h"
#include "./DISFromSpline.h"
#include "./CharmDISFromSpline.h"
#include "./QuarkDISFromSpline.h"
#include "./HNLFromSpline.h"
#include "./Decay.h"
#include "./NeutrissimoDecay.h"
#include "./InteractionCollection.h"
#include "./DummyCrossSection.h"
#include "./Hadronization.h"
#include "./CharmHadronization.h"
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
    register_Hadronization(m);

    register_CharmHadronization(m);
    register_CharmMesonDecay(m);
    register_CharmMesonDecay3Body(m);
    register_DMesonELoss(m);

    register_DipoleFromTable(m);
    register_DarkNewsCrossSection(m);
    register_DarkNewsDecay(m);
    register_DISFromSpline(m);
    register_CharmDISFromSpline(m);
    register_QuarkDISFromSpline(m);
    register_HNLFromSpline(m);
    register_NeutrissimoDecay(m);
    register_InteractionCollection(m);
    register_DummyCrossSection(m);
    register_PythiaDISCrossSection(m);
    register_MarleyCrossSection(m);
}
