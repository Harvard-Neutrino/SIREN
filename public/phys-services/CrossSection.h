#ifndef LI_CrossSection_H
#define LI_CrossSection_H

#include <map>
#include <set>
#include <array>
#include <string>

#include <photospline/splinetable.h>

#include "LeptonInjector/Particle.h"
#include "LeptonInjector/Random.h"

namespace LeptonInjector {

struct InteractionSignature {
    LeptonInjector::Particle::ParticleType primary_type;
    LeptonInjector::Particle::ParticleType target_type;
    std::vector<LeptonInjector::Particle::ParticleType> secondary_types;
};

struct InteractionRecord {
    InteractionSignature signature;
    std::array<double, 4> primary_momentum = {0, 0, 0, 0};
    std::array<double, 4> target_momentum = {0, 0, 0, 0};
    std::vector<std::array<double, 4>> secondary_momenta;
};

class CrossSection {
private:
    //
public:
    CrossSection() {};
    virtual double TotalCrossSection(InteractionRecord const &) const = 0;
    virtual double TotalCrossSection(LeptonInjector::Particle::ParticleType primary, double energy) const = 0;
    virtual double DifferentialCrossSection(InteractionRecord const &) const = 0;
    virtual void SampleFinalState(InteractionRecord &, std::shared_ptr<LeptonInjector::LI_random>) const = 0;

    virtual std::vector<Particle::ParticleType> GetPossibleTargets() const = 0;
    virtual std::vector<Particle::ParticleType> GetPossiblePrimaries() const = 0;
    virtual std::vector<InteractionSignature> GetPossibleSignatures() const = 0;
};

class DISFromSpline : public CrossSection {
private:
    photospline::splinetable<> differential_cross_section_;
    photospline::splinetable<> total_cross_section_;

    std::map<std::pair<LeptonInjector::Particle::ParticleType, LeptonInjector::Particle::ParticleType>, InteractionSignature> signatures_;
    std::set<LeptonInjector::Particle::ParticleType> primary_types_;
    std::set<LeptonInjector::Particle::ParticleType> target_types_;

    double minimum_Q2_;
    double target_mass_;
    int interaction_type_;

public:
    DISFromSpline(std::string differential_filename, std::string total_filename, int interaction, double target_mass, double minumum_Q2, std::set<LeptonInjector::Particle::ParticleType> primary_types, std::set<LeptonInjector::Particle::ParticleType> target_types);
    DISFromSpline(std::string differential_filename, std::string total_filename, std::set<LeptonInjector::Particle::ParticleType> primary_types, std::set<LeptonInjector::Particle::ParticleType> target_types);
    DISFromSpline(std::string differential_filename, std::string total_filename, int interaction, double target_mass, double minumum_Q2, std::vector<LeptonInjector::Particle::ParticleType> primary_types, std::vector<LeptonInjector::Particle::ParticleType> target_types);
    DISFromSpline(std::string differential_filename, std::string total_filename, std::vector<LeptonInjector::Particle::ParticleType> primary_types, std::vector<LeptonInjector::Particle::ParticleType> target_types);

    double TotalCrossSection(InteractionRecord const &) const;
    double TotalCrossSection(LeptonInjector::Particle::ParticleType primary, double energy) const;
    double DifferentialCrossSection(InteractionRecord const &) const;
    double DifferentialCrossSection(double energy, double x, double y, double secondary_lepton_mass) const;
    void SampleFinalState(InteractionRecord &, std::shared_ptr<LeptonInjector::LI_random> random) const;

    std::vector<Particle::ParticleType> GetPossibleTargets() const;
    std::vector<Particle::ParticleType> GetPossiblePrimaries() const;
    std::vector<InteractionSignature> GetPossibleSignatures() const;

    void LoadFromFile(std::string differential_filename, std::string total_filename);

private:
    void ReadParamsFromSplineTable();
    void InitializeSignatures();
};

template<typename T>
struct TableData2D {
    std::vector<T> x;
    std::vector<T> f;
};

template<typename T>
struct TableData2D {
    std::vector<T> x;
    std::vector<T> y;
    std::vector<T> f;
};

template<typename T>
struct IndexFinderIrregular {
    std::set data;
    IndexFinderIrregular(std::set<T> x) data(x){};

    unsigned int operator()(T const & x) {
        // Lower bound returns pointer to element that is greater than or equal to x
        // i.e. x \in (a,b] --> pointer to b, x \in (b,c] --> pointer to c
        // begin is the first element
        // distance(begin, pointer to y) --> y
        // therefore this function returns the index of the upper bin edge
        return std::distance(data.begin(), data.lower_bound(x));
    }
};

template<typename T>
struct IndexFinderRegular {
    T low;
    T high;
    T range;
    unsigned int n_points;
    T delta;

    IndexFinderIrregular(std::set<T> x) {
        n_points = x.size();
        low = *x.begin();
        high = *(x.last() - 1);
        range = high - low;
        delta = range / (n_points - 1);
    };

    unsigned int operator()(T const & x) {
        int i = (int)alt_floor<DataType>()((x - low) / range * (n_points - 1));
        return i;
    }
};



template<typename T>
struct Interpolator1D {
private:
    T min_x;
    T range;
    unsigned int n_points;
    T delta;

    bool is_log = false;
    bool is_regular = true;

    IndexFinderRegular regular_index;
    IndexFinderIrregular irregular_index;

	std::map<unsigned int, T> function;

public:

    EquispacedLogInterpolator():
    }

    static
    T MaxDist(std::set<T> x, T avg_diff) {
        std::vector<T> dist(x.size() - 1);
        for(unsigned int i=1; i<x.size(); ++i) {
            dist[i-1] = std::abs(x[i] - x[i-1] - avg_diff);
        }
        return *std::max_element(dist.first(), dist.last());
    };

	void Add1DTable(Table1DData & table) {
        std::set<T> x(table.x.begin(). table.x.end());
        std::map<T, unsigned int> xmap;
        for(unsigned int n = 0; auto i : x) {
            xmap[i] = n;
            ++n;
        }

        regular_index = IndexFinderRegular(x);
        T relative_max_dist = MaxDist(x, regular_index.delta) / regular_index.delta;
        T log_relative_max_dist;

        std::vector<T> log_x(x.begin(), x.end());
        std::transform(log_x.begin(), log_x.end(), log_x.begin() [](T t)->T{return log(t);});
        std::set<T> log_x_set(log_x.begin(), log_x.end());

        if(relative_max_dist < 1e-4) {
            is_regular = true;
            is_log = false;
        }

        if(not is_regular) {
            regular_index = IndexFinderRegular(log_x_set);
            log_relative_max_dist = MaxDist(log_x_set, regular_index.delta) / regular_index.delta;
            if(log_relative_max_dist < 1e-4) {
                is_regular = true;
                if_log = true;
            }
        }

        if(not is_regular) {
            is_log = log_relative_max_dist < relative_max_dist;
            if(is_log) {
                irregular_index = IndexFinderIrregular(log_x_set);
            } else {
                irregular_index = IndexFinderIrregular(x);
            }
        }

        std::vector<T> function_values(table.f.begin(), table.f.end());

        if(is_regular) {
            min_x = irregular_index.low;
            delta = irregular_index.delta;
        } else {
            min_x = regular_index.low;
            delta = regular_index.delta;
        }

        if(is_log) {
            std::transform(function_values.begin(), function_values.end(), function_values.begin() [](T t)->T{return log(t);});
        }

        for(unsigned int i=0; i<table.x.size(); ++i) {
            differential_cross_section_[
                        xmap[table.x[i]]
            ] = table.f[i];
        }
    }

	void Add2DTable(Table2DData & table) {
        std::set<T> x(table.x.begin(). table.x.end());
        std::set<T> y(table.y.begin(). table.y.end());
        std::map<T, unsigned int> xmap;
        std::map<T, unsigned int> ymap;
        for(unsigned int n = 0; auto i : x) {
            xmap[i] = n;
            ++n;
        }
        for(unsigned int n = 0; auto i : y) {
            ymap[i] = n;
            ++n;
        }
        for(unsigned int i=0; i<table.x.size(); ++i) {
            differential_cross_section_[
                std::pair<unsigned int, unsigned int>(
                        xmap[table.x[i]],
                        ymap[table.y[i]]
                        )
            ] = table.f[i];
        }
    }

    template<typename DataType>
    DataType operator()(DataType const & x) {
        DataType val = x;

        if(is_log) {
            val = log(val);
        }

        unsigned int index;

        if(is_regular) {
            index = regular_index(x);
        } else {
            index = irregular_index(x);
        }

        T result;

        if(index < 0) {
            result = DataType(function[0]);
        }
        else if(index >= values.size() - 1) {
            result = DataType(function[function.size() - 1]);
        }
        else {
            T ya = function[index];
            T yb = function[i+1];
            T xa = min_x + index*delta;
            result = (DataType(ya + (yb - ya) * (x - xa) / delta));
        }

        if(is_log) {
            result = exp(result);
        }

        return result;
    }
};


class DipoleFromTable : public CrossSection {
public:
private:
    std::map<Particle::ParticleType, DifferentialTableData> differential
    double hnl_mass;
public:

    DipoleFromTable(double target_mass) {};
    double TotalCrossSection(InteractionRecord const &) const;
    double TotalCrossSection(LeptonInjector::Particle::ParticleType primary, double energy) const;
    double DifferentialCrossSection(InteractionRecord const &) const;
    void SampleFinalState(InteractionRecord &, std::shared_ptr<LeptonInjector::LI_random>) const;

    std::vector<Particle::ParticleType> GetPossibleTargets() const;
    std::vector<Particle::ParticleType> GetPossiblePrimaries() const;
    std::vector<InteractionSignature> GetPossibleSignatures() const;
private:
    AddDifferentialCrossSectionFile(std::string filename, Particle::ParticleType target);
    AddTotalCrossSectionFile(std::string filename, Particle::ParticleType target);
};

} // namespace LeptonInjector

#endif // LI_CrossSection_H

