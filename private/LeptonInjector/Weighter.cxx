#include <cassert>
#include <fstream>
#include <algorithm>

#include "LeptonInjector/LeptonInjector.h"
#include "LeptonInjector/EventProps.h"

#include "phys-services/CrossSection.h"

#include "earthmodel-service/Path.h"

#include <rk/rk.hh>

namespace LeptonInjector{

//---------------
// class LeptonWeighter
//---------------

namespace {
    template <class InIt>
    typename std::iterator_traits<InIt>::value_type accumulate(InIt begin, InIt end) {
        typedef typename std::iterator_traits<InIt>::value_type real;
        real sum = real(0);
        real running_error = real(0);
        real temp;
        real difference;

        for (; begin != end; ++begin) {
            difference = *begin;
            difference -= running_error;
            temp = sum;
            temp += difference;
            running_error = temp;
            running_error -= sum;
            running_error -= difference;
            sum = std::move(temp);
        }
        return sum;
    };

    template<typename T>
    T accumulate(std::initializer_list<T> list) {
        return accumulate(list.begin(), list.end());
    }
}

double LeptonWeighter::InteractionProbability(std::shared_ptr<InjectorBase const> injector, InteractionRecord const & record) const {
    std::pair<earthmodel::Vector3D, earthmodel::Vector3D> bounds = injector->InjectionBounds(record);
    return InteractionProbability(bounds, record);
}

double LeptonWeighter::InteractionProbability(std::pair<earthmodel::Vector3D, earthmodel::Vector3D> bounds, InteractionRecord const & record) const {
    std::map<Particle::ParticleType, std::vector<std::shared_ptr<CrossSection>>> const & cross_sections_by_target = cross_sections->GetCrossSectionsByTarget();
    std::vector<double> interaction_probabilities;
    GetColumnDepthInCGS(Vector3D const & p0, Vector3D const & p1, std::set<LeptonInjector::Particle::ParticleType> targets)
}

LeptonWeighter::LeptonWeighter(std::vector<std::shared_ptr<InjectorBase>> injectors, std::shared_ptr<earthmodel::EarthModel> earth_model, std::shared_ptr<CrossSectionCollection> cross_sections) {
    
}

double LeptonWeighter::EventWeight(InteractionRecord const & record) const {

}

} // namespace LeptonInjector

