#ifndef LI_LeptonInjector_H
#define LI_LeptonInjector_H

#include <queue>
#include <memory> // adds shared pointer
#include <iostream>

#include <photospline/bspline.h>
#include <photospline/splinetable.h>
#include <photospline/cinter/splinetable.h>

#include <cereal/cereal.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/array.hpp>
#include <cereal/types/set.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>
#include "serialization/array.h"

#include "LeptonInjector/Random.h"
#include "LeptonInjector/Particle.h"
#include "LeptonInjector/Constants.h"
#include "LeptonInjector/DataWriter.h"
#include "LeptonInjector/EventProps.h"
#include "LeptonInjector/Coordinates.h"
#include "LeptonInjector/Distributions.h"
#include "LeptonInjector/LeptonInjector.h"
#include "LeptonInjector/BasicInjectionConfiguration.h"

#include "phys-services/CrossSection.h"

#include "earthmodel-service/Path.h"
#include "earthmodel-service/Vector3D.h"
#include "earthmodel-service/EarthModel.h"

namespace LeptonInjector {

class LeptonWeighter {
private:
    std::vector<std::shared_ptr<InjectorBase>> injectors;
    std::shared_ptr<earthmodel::EarthModel> earth_model;
    std::shared_ptr<CrossSectionCollection> cross_sections;
    double InteractionProbability(std::shared_ptr<InjectorBase const> injector, InteractionRecord const & record) const;
    double InteractionProbability(std::pair<earthmodel::Vector3D, earthmodel::Vector3D> bounds, InteractionRecord const & record) const;
public:
    LeptonWeighter(std::vector<std::shared_ptr<InjectorBase>> injectors, std::shared_ptr<earthmodel::EarthModel> earth_model, std::shared_ptr<CrossSectionCollection> cross_sections);
    double EventWeight(InteractionRecord const & record) const;
}

} //namespace LeptonInjector

CEREAL_CLASS_VERSION(LeptonInjector::LeptonWeighter, 0);

#endif // LI_LeptonInjector_H

