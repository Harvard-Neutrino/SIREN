#include "SIREN/injection/TreeWeighter.h"

#include <iterator>                                              // for ite...
#include <array>                                                  // for array
#include <cassert>                                                // for assert
#include <cmath>                                                  // for exp
#include <initializer_list>                                       // for ini...
#include <iostream>                                               // for ope...
#include <set>                                                    // for set
#include <stdexcept>
#include <tuple>
#include <cassert>
#include <fstream>
#include <algorithm>
                                             // for out...
#include "SIREN/interactions/CrossSection.h"            // for Cro...
#include "SIREN/interactions/InteractionCollection.h"  // for Cro...
#include "SIREN/dataclasses/InteractionRecord.h"         // for Int...
#include "SIREN/dataclasses/InteractionSignature.h"      // for Int...
#include "SIREN/detector/DetectorModel.h"                   // for Ear...
#include "SIREN/detector/Coordinates.h"
#include "SIREN/distributions/Distributions.h"           // for Inj...
#include "SIREN/geometry/Geometry.h"                     // for Geo...
#include "SIREN/injection/Injector.h"                // for Inj...
#include "SIREN/injection/Process.h"                     // for Phy...
#include "SIREN/injection/WeightingUtils.h"              // for Cro...
#include "SIREN/math/Vector3D.h"                         // for Vec...


#include "SIREN/injection/Injector.h"

#include "SIREN/distributions/primary/vertex/VertexPositionDistribution.h"

#include "SIREN/interactions/CrossSection.h"
#include "SIREN/interactions/InteractionCollection.h"

#include "SIREN/dataclasses/InteractionSignature.h"

#include <rk/rk.hh>

namespace siren {
namespace injection {

using detector::DetectorPosition;
using detector::DetectorDirection;

//---------------
// class TreeWeighter
//---------------

void TreeWeighter::Initialize() {
    int i = 0;
    primary_process_weighters.reserve(injectors.size());
    secondary_process_weighter_maps.reserve(injectors.size());
    for(auto const & injector : injectors) {
        assert(primary_physical_process->MatchesHead(injector->GetPrimaryProcess()));
        primary_process_weighters.push_back(std::make_shared<PrimaryProcessWeighter>(PrimaryProcessWeighter(primary_physical_process, injector->GetPrimaryProcess(), detector_model)));
        std::map<siren::dataclasses::ParticleType, std::shared_ptr<SecondaryProcessWeighter>>
            injector_sec_process_weighter_map;
        std::map<siren::dataclasses::ParticleType, std::shared_ptr<siren::injection::SecondaryInjectionProcess>>
            injector_sec_process_map = injector->GetSecondaryProcessMap();
        for(auto const & sec_phys_process : secondary_physical_processes) {
            try{
                std::shared_ptr<siren::injection::SecondaryInjectionProcess> sec_inj_process = injector_sec_process_map.at(sec_phys_process->GetPrimaryType());
                assert(sec_phys_process->MatchesHead(sec_inj_process)); // make sure cross section collection matches
                injector_sec_process_weighter_map[sec_phys_process->GetPrimaryType()] =
                    std::make_shared<SecondaryProcessWeighter>(
                            SecondaryProcessWeighter(
                                sec_phys_process,
                                sec_inj_process,detector_model
                            )
                    );
            } catch(const std::out_of_range& oor) {
                std::cout << "Out of Range error: " << oor.what() << '\n';
                std::cout << "Initialization Incomplete: Particle " <<  sec_phys_process->GetPrimaryType() << " does not exist in injector\n";
                return;
            }
        }
        if(injector_sec_process_weighter_map.size() != injector_sec_process_map.size()) {
            std::cout << "Initialization Incomplete: No one-to-one mapping between injection and physical distributions for injector " << i << "\n";
            return;
        }
        secondary_process_weighter_maps.push_back(injector_sec_process_weighter_map);
    }
}

double TreeWeighter::EventWeight(siren::dataclasses::InteractionTree const & tree) const {
    // The weight is given by
    //
    // [sum_{injectors i}
    //  x prod_{tree datum d}
    //  x (prod_{generation dist j} p_gen^{idj})
    //  / (prod_{physical dist j} p_phys^{idj}) ] ^-1
    //
    double inv_weight = 0;
    for(unsigned int idx = 0; idx < injectors.size(); ++idx) {
        double physical_probability = 1.0;
        double generation_probability = injectors[idx]->EventsToInject();//GenerationProbability(tree);
        for(auto const & datum : tree.tree) {
            std::tuple<siren::math::Vector3D, siren::math::Vector3D> bounds;
            if(datum->depth() == 0) {
                bounds = injectors[idx]->PrimaryInjectionBounds(datum->record);
                physical_probability *= primary_process_weighters[idx]->PhysicalProbability(bounds, datum->record);
                generation_probability *= primary_process_weighters[idx]->GenerationProbability(*datum);
            }
            else {
                try {
                    bounds = injectors[idx]->SecondaryInjectionBounds(datum->record);
                    double phys_prob = secondary_process_weighter_maps[idx].at(datum->record.signature.primary_type)->PhysicalProbability(bounds, datum->record);
                    double gen_prob = secondary_process_weighter_maps[idx].at(datum->record.signature.primary_type)->GenerationProbability(*datum);
                    physical_probability *= phys_prob;
                    generation_probability *= gen_prob;
                } catch(const std::out_of_range& oor) {
                    std::cout << "Out of Range error: " << oor.what() << '\n';
                    return 0;
                }
            }
        }
        inv_weight += generation_probability / physical_probability;
    }
    return 1./inv_weight;
}

void TreeWeighter::SaveWeighter(std::string const & filename) const {
    std::ofstream os(filename+".siren_weighter", std::ios::binary);
    ::cereal::BinaryOutputArchive archive(os);
    this->save(archive,0);
}

void TreeWeighter::LoadWeighter(std::string const & filename) {
    std::cout << "Weighter loading not yet supported... sorry!\n";
    exit(0);
    std::ifstream is(filename+".siren_weighter", std::ios::binary);
    ::cereal::BinaryInputArchive archive(is);
    //this->load(archive,0);
}

TreeWeighter::TreeWeighter(std::vector<std::shared_ptr<Injector>> injectors, std::shared_ptr<siren::detector::DetectorModel> detector_model, std::shared_ptr<siren::injection::PhysicalProcess> primary_physical_process, std::vector<std::shared_ptr<siren::injection::PhysicalProcess>> secondary_physical_processes)
    : injectors(injectors)
      , detector_model(detector_model)
      , primary_physical_process(primary_physical_process)
      , secondary_physical_processes(secondary_physical_processes)
{
    Initialize();
}

TreeWeighter::TreeWeighter(std::vector<std::shared_ptr<Injector>> injectors, std::shared_ptr<siren::detector::DetectorModel> detector_model, std::shared_ptr<siren::injection::PhysicalProcess> primary_physical_process)
    : injectors(injectors)
      , detector_model(detector_model)
      , primary_physical_process(primary_physical_process)
      , secondary_physical_processes(std::vector<std::shared_ptr<siren::injection::PhysicalProcess>>())
{
    Initialize();
}

TreeWeighter::TreeWeighter(std::vector<std::shared_ptr<Injector>> _injectors, std::string filename) {
    LoadWeighter(filename);
    if(_injectors.size() > 0) {
        // overwrite the serialized injectors if the user have provided any
        injectors = _injectors;
    }
    Initialize();
}

} // namespace injection
} // namespace siren
