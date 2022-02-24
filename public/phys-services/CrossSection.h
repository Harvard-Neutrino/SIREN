#pragma once
#ifndef LI_CrossSection_H
#define LI_CrossSection_H

#include <map>
#include <set>
#include <array>
#include <string>

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

#include "LeptonInjector/Particle.h"
#include "LeptonInjector/Random.h"

namespace LeptonInjector {

struct InteractionSignature {
    LeptonInjector::Particle::ParticleType primary_type;
    LeptonInjector::Particle::ParticleType target_type;
    std::vector<LeptonInjector::Particle::ParticleType> secondary_types;
    bool operator==(InteractionSignature const & other) const;
    friend std::ostream& operator<<(std::ostream& os, InteractionSignature const& signature);
    template<class Archive>
    void serialize(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(cereal::make_nvp("PrimaryType", primary_type));
            archive(cereal::make_nvp("TargetType", target_type));
            archive(cereal::make_nvp("SecondaryTypes", secondary_types));
        } else {
            throw std::runtime_error("InteractionSignature only supports version <= 0!");
        }
    }
};

struct InteractionRecord {
    InteractionSignature signature;
    double primary_mass = 0;
    std::array<double, 4> primary_momentum = {0, 0, 0, 0};
    double primary_helicity = 0;
    double target_mass = 0;
    std::array<double, 4> target_momentum = {0, 0, 0, 0};
    double target_helicity = 0;
    std::array<double, 3> interaction_vertex = {0, 0, 0};
    std::vector<double> secondary_masses;
    std::vector<std::array<double, 4>> secondary_momenta;
    std::vector<double> secondary_helicity;
    std::vector<double> interaction_parameters;
    friend std::ostream& operator<<(std::ostream& os, InteractionRecord const& record);
    template<class Archive>
    void serialize(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(::cereal::make_nvp("InteractionSignature", signature));
            archive(::cereal::make_nvp("PrimaryMass", primary_mass));
            archive(::cereal::make_nvp("PrimaryMomentum", primary_momentum));
            archive(::cereal::make_nvp("PrimaryHelicity", primary_helicity));
            archive(::cereal::make_nvp("TargetMass", target_mass));
            archive(::cereal::make_nvp("TargetMomentum", target_momentum));
            archive(::cereal::make_nvp("TargetHelicity", target_helicity));
            archive(::cereal::make_nvp("InteractionVertex", interaction_vertex));
            archive(::cereal::make_nvp("SecondaryMasses", secondary_masses));
            archive(::cereal::make_nvp("SecondaryMomenta", secondary_momenta));
            archive(::cereal::make_nvp("SecondaryHelicity", secondary_helicity));
            archive(::cereal::make_nvp("InteractionParameters", interaction_parameters));
        } else {
            throw std::runtime_error("InteractionRecord only supports version <= 0!");
        }
    };
};

struct DecaySignature {
    LeptonInjector::Particle::ParticleType primary_type;
    std::vector<LeptonInjector::Particle::ParticleType> secondary_types;
    bool operator==(DecaySignature const & other) const;
    friend std::ostream& operator<<(std::ostream& os, DecaySignature const& signature);
    template<class Archive>
    void serialize(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(cereal::make_nvp("PrimaryType", primary_type));
            archive(cereal::make_nvp("SecondaryTypes", secondary_types));
        } else {
            throw std::runtime_error("DecaySignature only supports version <= 0!");
        }
    }
};

struct DecayRecord {
    DecaySignature signature;
    double primary_mass;
    std::array<double, 4> primary_momentum;
    double primary_helicity;
    std::array<double, 3> decay_vertex = {0, 0, 0};
    std::vector<double> secondary_masses;
    std::vector<std::array<double, 4>> secondary_momenta;
    std::vector<double> secondary_helicity;
    std::vector<double> decay_parameters;
    friend std::ostream& operator<<(std::ostream& os, DecayRecord const& record);
    template<class Archive>
    void serialize(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(::cereal::make_nvp("DecaySignature", signature));
            archive(::cereal::make_nvp("PrimaryMass", primary_mass));
            archive(::cereal::make_nvp("PrimaryMomentum", primary_momentum));
            archive(::cereal::make_nvp("PrimaryHelicity", primary_helicity));
            archive(::cereal::make_nvp("DecayVertex", decay_vertex));
            archive(::cereal::make_nvp("SecondaryMasses", secondary_masses));
            archive(::cereal::make_nvp("SecondaryMomenta", secondary_momenta));
            archive(::cereal::make_nvp("SecondaryHelicity", secondary_helicity));
            archive(::cereal::make_nvp("DecayParameters", decay_parameters));
        } else {
            throw std::runtime_error("DecayRecord only supports version <= 0!");
        }
    };
};

class CrossSection {
friend cereal::access;
private:
public:
    CrossSection() {};
    virtual double TotalCrossSection(InteractionRecord const &) const = 0;
    virtual double TotalCrossSection(LeptonInjector::Particle::ParticleType primary, double energy, Particle::ParticleType target) const = 0;
    virtual double DifferentialCrossSection(InteractionRecord const &) const = 0;
    virtual double InteractionThreshold(InteractionRecord const &) const = 0;
    virtual void SampleFinalState(InteractionRecord &, std::shared_ptr<LeptonInjector::LI_random>) const = 0;

    virtual std::vector<Particle::ParticleType> GetPossibleTargets() const = 0;
    virtual std::vector<Particle::ParticleType> GetPossibleTargetsFromPrimary(Particle::ParticleType primary_type) const = 0;
    virtual std::vector<Particle::ParticleType> GetPossiblePrimaries() const = 0;
    virtual std::vector<InteractionSignature> GetPossibleSignatures() const = 0;

    virtual std::vector<InteractionSignature> GetPossibleSignaturesFromParents(Particle::ParticleType primary_type, Particle::ParticleType target_type) const = 0;
    virtual double FinalStateProbability(InteractionRecord const & record) const = 0;
    virtual std::vector<std::string> DensityVariables() const = 0;
    template<class Archive>
    void save(Archive & archive, std::uint32_t const version) const {};
    template<class Archive>
    void load(Archive & archive, std::uint32_t const version) {};
};

class CrossSectionCollection {
private:
    Particle::ParticleType primary_type;
    std::vector<std::shared_ptr<CrossSection>> cross_sections;
    std::map<Particle::ParticleType, std::vector<std::shared_ptr<CrossSection>>> cross_sections_by_target;
    std::set<Particle::ParticleType> target_types;
    static const std::vector<std::shared_ptr<CrossSection>> empty;
    void InitializeTargetTypes();
public:
    CrossSectionCollection();
    CrossSectionCollection(Particle::ParticleType primary_type, std::vector<std::shared_ptr<CrossSection>> cross_sections);
    std::vector<std::shared_ptr<CrossSection>> const & GetCrossSections() const {return cross_sections;};
    std::vector<std::shared_ptr<CrossSection>> const & GetCrossSectionsForTarget(Particle::ParticleType p) const;
    std::map<Particle::ParticleType, std::vector<std::shared_ptr<CrossSection>>> const & GetCrossSectionsByTarget() const {
        return cross_sections_by_target;
    };
    std::set<Particle::ParticleType> const & TargetTypes() const {
        return target_types;
    };
    virtual bool MatchesPrimary(InteractionRecord const & record) const;
public:
    template<class Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(cereal::make_nvp("PrimaryType", primary_type));
            archive(cereal::make_nvp("CrossSections", cross_sections));
        } else {
            throw std::runtime_error("CrossSectionCollection only supports version <= 0!");
        }
    }

    template<class Archive>
    void load(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(cereal::make_nvp("PrimaryType", primary_type));
            archive(cereal::make_nvp("CrossSections", cross_sections));
        } else {
            throw std::runtime_error("CrossSectionCollection only supports version <= 0!");
        }
    }
};

class DISFromSpline : public CrossSection {
friend cereal::access;
private:
    photospline::splinetable<> differential_cross_section_;
    photospline::splinetable<> total_cross_section_;

    std::vector<InteractionSignature> signatures_;
    std::set<LeptonInjector::Particle::ParticleType> primary_types_;
    std::set<LeptonInjector::Particle::ParticleType> target_types_;
    std::map<LeptonInjector::Particle::ParticleType, std::vector<LeptonInjector::Particle::ParticleType>> targets_by_primary_types_;
    std::map<std::pair<LeptonInjector::Particle::ParticleType, LeptonInjector::Particle::ParticleType>, std::vector<InteractionSignature>> signatures_by_parent_types_;

    int interaction_type_;
    double target_mass_;
    double minimum_Q2_;

public:
    DISFromSpline() {};
    DISFromSpline(std::vector<char> differential_data, std::vector<char> total_data, int interaction, double target_mass, double minumum_Q2, std::set<LeptonInjector::Particle::ParticleType> primary_types, std::set<LeptonInjector::Particle::ParticleType> target_types);
    DISFromSpline(std::vector<char> differential_data, std::vector<char> total_data, int interaction, double target_mass, double minumum_Q2, std::vector<LeptonInjector::Particle::ParticleType> primary_types, std::vector<LeptonInjector::Particle::ParticleType> target_types);
    DISFromSpline(std::string differential_filename, std::string total_filename, int interaction, double target_mass, double minumum_Q2, std::set<LeptonInjector::Particle::ParticleType> primary_types, std::set<LeptonInjector::Particle::ParticleType> target_types);
    DISFromSpline(std::string differential_filename, std::string total_filename, std::set<LeptonInjector::Particle::ParticleType> primary_types, std::set<LeptonInjector::Particle::ParticleType> target_types);
    DISFromSpline(std::string differential_filename, std::string total_filename, int interaction, double target_mass, double minumum_Q2, std::vector<LeptonInjector::Particle::ParticleType> primary_types, std::vector<LeptonInjector::Particle::ParticleType> target_types);
    DISFromSpline(std::string differential_filename, std::string total_filename, std::vector<LeptonInjector::Particle::ParticleType> primary_types, std::vector<LeptonInjector::Particle::ParticleType> target_types);

    double TotalCrossSection(InteractionRecord const &) const;
    double TotalCrossSection(LeptonInjector::Particle::ParticleType primary, double energy) const;
    double TotalCrossSection(LeptonInjector::Particle::ParticleType primary, double energy, Particle::ParticleType target) const;
    double DifferentialCrossSection(InteractionRecord const &) const;
    double DifferentialCrossSection(double energy, double x, double y, double secondary_lepton_mass) const;
    double InteractionThreshold(InteractionRecord const &) const;
    void SampleFinalState(InteractionRecord &, std::shared_ptr<LeptonInjector::LI_random> random) const;

    std::vector<Particle::ParticleType> GetPossibleTargets() const;
    std::vector<Particle::ParticleType> GetPossibleTargetsFromPrimary(Particle::ParticleType primary_type) const;
    std::vector<Particle::ParticleType> GetPossiblePrimaries() const;
    std::vector<InteractionSignature> GetPossibleSignatures() const;
    std::vector<InteractionSignature> GetPossibleSignaturesFromParents(Particle::ParticleType primary_type, Particle::ParticleType target_type) const;

    virtual double FinalStateProbability(InteractionRecord const & record) const;

    void LoadFromFile(std::string differential_filename, std::string total_filename);
    void LoadFromMemory(std::vector<char> & differential_data, std::vector<char> & total_data);

    double GetMinimumQ2() const {return minimum_Q2_;};
    double GetTargetMass() const {return target_mass_;};
    int GetInteractionType() const {return interaction_type_;};

public:
    virtual std::vector<std::string> DensityVariables() const override;
    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            splinetable_buffer buf;
            buf.size = 0;
            auto result_obj = differential_cross_section_.write_fits_mem();
            buf.data = result_obj.first;
            buf.size = result_obj.second;

            std::vector<char> diff_blob;
            diff_blob.resize(buf.size);
            std::copy((char*)buf.data, (char*)buf.data + buf.size, &diff_blob[0]);

            archive(::cereal::make_nvp("DifferentialCrossSectionSpline", diff_blob));

            buf.size = 0;
            result_obj = total_cross_section_.write_fits_mem();
            buf.data = result_obj.first;
            buf.size = result_obj.second;

            std::vector<char> total_blob;
            total_blob.resize(buf.size);
            std::copy((char*)buf.data, (char*)buf.data + buf.size, &total_blob[0]);

            archive(::cereal::make_nvp("TotalCrossSectionSpline", total_blob));
            archive(::cereal::make_nvp("PrimaryTypes", primary_types_));
            archive(::cereal::make_nvp("TargetTypes", target_types_));
            archive(::cereal::make_nvp("InteractionType", interaction_type_));
            archive(::cereal::make_nvp("TargetMass", target_mass_));
            archive(::cereal::make_nvp("MinimumQ2", minimum_Q2_));
            archive(cereal::virtual_base_class<CrossSection>(this));
        } else {
            throw std::runtime_error("DISFromSpline only supports version <= 0!");
        }
    }
    template<typename Archive>
    void load(Archive & archive, std::uint32_t version) {
        if(version == 0) {
            std::vector<char> differential_data;
            std::vector<char> total_data;
            archive(::cereal::make_nvp("DifferentialCrossSectionSpline", differential_data));
            archive(::cereal::make_nvp("TotalCrossSectionSpline", total_data));
            archive(::cereal::make_nvp("PrimaryTypes", primary_types_));
            archive(::cereal::make_nvp("TargetTypes", target_types_));
            archive(::cereal::make_nvp("InteractionType", interaction_type_));
            archive(::cereal::make_nvp("TargetMass", target_mass_));
            archive(::cereal::make_nvp("MinimumQ2", minimum_Q2_));
            LoadFromMemory(differential_data, total_data);
            InitializeSignatures();
            archive(cereal::virtual_base_class<CrossSection>(this));
        } else {
            throw std::runtime_error("DISFromSpline only supports version <= 0!");
        }
    }
private:
    void ReadParamsFromSplineTable();
    void InitializeSignatures();
};


template<typename T>
struct TableData1D {
    std::vector<T> x;
    std::vector<T> f;
    template<class Archive>
    void serialize(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(cereal::make_nvp("TableX", x));
            archive(cereal::make_nvp("TableF", f));
        } else {
            throw std::runtime_error("TableData1D only supports version <= 0!");
        }
    }
};

template<typename T>
struct TableData2D {
    std::vector<T> x;
    std::vector<T> y;
    std::vector<T> f;
    template<class Archive>
    void serialize(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(cereal::make_nvp("TableX", x));
            archive(cereal::make_nvp("TableY", y));
            archive(cereal::make_nvp("TableF", f));
        } else {
            throw std::runtime_error("TableData2D only supports version <= 0!");
        }
    }
};

template<typename T>
struct IndexFinderIrregular {
    std::vector<T> data;
    std::vector<T> diff;
    T low;
    T high;
    T range;
    unsigned int n_points;

    IndexFinderIrregular() {};
    IndexFinderIrregular(std::set<T> x): data(x.begin(), x.end()) {
        std::sort(data.begin(), data.end());
        low = data.front();
        high = data.back();
        range = high - low;
        diff.resize(data.size() - 1);
        for(unsigned int i=1; i<data.size(); ++i) {
            diff[i-1] = data[i] - data[i-1];
        }
        n_points = data.size();
    };

    template<class Archive>
    void serialize(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(cereal::make_nvp("Data", data));
            archive(cereal::make_nvp("Diff", diff));
            archive(cereal::make_nvp("Low", low));
            archive(cereal::make_nvp("High", high));
            archive(cereal::make_nvp("Range", range));
            archive(cereal::make_nvp("NPoints", n_points));
        } else {
            throw std::runtime_error("IndexFinderIrregular only supports version <= 0!");
        }
    }

    std::tuple<unsigned int, T, T, T> operator()(T const & x) const {
        // Lower bound returns pointer to element that is greater than or equal to x
        // i.e. x \in (a,b] --> pointer to b, x \in (b,c] --> pointer to c
        // begin is the first element
        // distance(begin, pointer to y) --> y
        // therefore this function returns the index of the lower bin edge
        unsigned int index = std::distance(data.begin(), std::lower_bound(data.begin(), data.end(), x)) - 1;
        if(index < 0)
            index = 0;
        else if(index >= n_points - 1)
            index = n_points - 2;
        return std::tuple<unsigned int, T, T, T>(index,
                x,
                data[index],
                diff[index]);
    }
};

template<typename DataType>
struct alt_floor {
    DataType operator()(DataType const & x) {
        return std::floor(x);
    }
};

template<typename T>
struct IndexFinderRegular {
    T low;
    T high;
    T range;
    unsigned int n_points;
    T delta;

    IndexFinderRegular() {};
    IndexFinderRegular(std::set<T> x) {
        std::vector<T> points(x.begin(), x.end());
        std::sort(points.begin(), points.end());
        n_points = points.size();
        low = points.front();
        high = points.back();
        range = high - low;
        delta = range / (n_points - 1);
    };

    template<class Archive>
    void serialize(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(cereal::make_nvp("Low", low));
            archive(cereal::make_nvp("High", high));
            archive(cereal::make_nvp("Range", range));
            archive(cereal::make_nvp("NPoints", n_points));
            archive(cereal::make_nvp("Delta", delta));
        } else {
            throw std::runtime_error("IndexFinderRegular only supports version <= 0!");
        }
    }

    std::tuple<unsigned int, T, T, T> operator()(T const & x) const {
        int i = (int)alt_floor<T>()((x - low) / range * (n_points - 1));
        if(i < 0)
            i = 0;
        else if(i >= n_points - 1)
            i = n_points - 2;
        T lower_edge = low + i * delta;
        return std::tuple<unsigned int, T, T, T>(i, x, lower_edge, delta);
    }
};


template<typename T>
struct Indexer1D {
private:
    T min_x;
    T max_x;
    T range;

    std::vector<T> points;

    bool is_log = true;
    bool is_regular = false;

    IndexFinderRegular<T> regular_index;
    IndexFinderIrregular<T> irregular_index;

public:

    Indexer1D() {};
    Indexer1D(TableData1D<T> & table) {
        AddTable(table);
    };

    template<class Archive>
    void serialize(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(cereal::make_nvp("MinX", min_x));
            archive(cereal::make_nvp("MaxX", max_x));
            archive(cereal::make_nvp("Range", range));
            archive(cereal::make_nvp("Points", points));
            archive(cereal::make_nvp("IsLog", is_log));
            archive(cereal::make_nvp("IsRegular", is_regular));
            archive(cereal::make_nvp("RegularIndexer", regular_index));
            archive(cereal::make_nvp("IrregularIndexer", irregular_index));
        } else {
            throw std::runtime_error("Indexer1D only supports version <= 0!");
        }
    }

    static
    T MaxDist(std::vector<T> x, T avg_diff) {
        std::vector<T> dist(x.size() - 1);
        for(unsigned int i=1; i<x.size(); ++i) {
            dist[i-1] = std::abs(std::abs(x[i] - x[i-1]) - avg_diff);
            if(std::isinf(dist[i-1])) {
                return std::numeric_limits<T>::infinity();
            }
        }
        return *std::max_element(dist.begin(), dist.end());
    };

    void AddTable(TableData1D<T> & table) {
        is_regular = false;
        std::set<T> x_set(table.x.begin(), table.x.end());
        std::vector<T> x(x_set.begin(), x_set.end());
        std::sort(x.begin(), x.end());
        unsigned int n_points = x.size();
        assert(n_points >= 2);

        std::vector<T> log_x(x.begin(), x.end());
        std::transform(log_x.begin(), log_x.end(), log_x.begin(), [](T t)->T{return log(t);});
        std::set<T> log_x_set(log_x.begin(), log_x.end());

        regular_index = IndexFinderRegular<T>(log_x_set);
        T relative_max_dist;

        T log_relative_max_dist = MaxDist(log_x, regular_index.delta) / regular_index.delta;

        if(log_relative_max_dist < 1e-4 and not std::isinf(regular_index.delta)) {
            is_regular = true;
            is_log = true;
        }

        if(not is_regular) {
            regular_index = IndexFinderRegular<T>(x_set);
            relative_max_dist = MaxDist(x, regular_index.delta) / regular_index.delta;
            if(relative_max_dist < 1e-4 and not std::isinf(regular_index.delta)) {
                is_regular = true;
                is_log = false;
            }
        }

        if(not is_regular) {
            is_log = log_relative_max_dist < relative_max_dist;
            if(is_log) {
                irregular_index = IndexFinderIrregular<T>(log_x_set);
            } else {
                irregular_index = IndexFinderIrregular<T>(x_set);
            }
        }

        if(is_log) {
            points = std::vector<T>(log_x_set.begin(), log_x_set.end());
        } else {
            points = std::vector<T>(x.begin(), x.end());
        }

        if(is_regular) {
            min_x = regular_index.low;
            max_x = regular_index.high;
            range = regular_index.range;
            irregular_index.data.clear();
        } else {
            min_x = irregular_index.low;
            max_x = irregular_index.high;
            range = irregular_index.range;
        }
        if(is_log) {
            min_x = exp(min_x);
            max_x = exp(max_x);
            range = max_x - min_x;
        }
    }

    std::tuple<unsigned int, T, T, T> operator()(T const & x) const {
        T val = x;

        if(is_log) {
            val = log(val);
        }

        if(is_regular) {
            return regular_index(val);
        } else {
            return irregular_index(val);
        }
    }

    T Min() const {
        return min_x;
    }
    T Range() const {
        return range;
    }
    T Max() const {
        return max_x;
    }
    bool IsRegular() const {
        return is_regular;
    }
    bool IsLog() const {
        return is_log;
    }
    unsigned int NPoints() const {
        return points.size();
    }
};

template<typename T>
struct Interpolator1D {
private:
    TableData1D<T> original_table;
    Indexer1D<T> indexer;
    std::map<unsigned int, T> function;
    std::vector<bool> zero_mask;
    bool is_log = false;
public:

    Interpolator1D() {};
    Interpolator1D(TableData1D<T> & table) {
        AddTable(table);
    };

    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(cereal::make_nvp("Table", original_table));
        } else {
            throw std::runtime_error("Interpolator1D only supports version <= 0!");
        }
    };
    template<typename Archive>
    void load(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            TableData1D<T> table;
            archive(cereal::make_nvp("Table", table));
            AddTable(table);
        } else {
            throw std::runtime_error("Interpolator1D only supports version <= 0!");
        }
    };

    void AddTable(TableData1D<T> & table) {
        original_table = table;
        std::set<T> x(table.x.begin(), table.x.end());
        std::map<T, unsigned int> xmap;
        unsigned int n = 0;
        for(auto i : x) {
            xmap[i] = n;
            ++n;
        }

        assert(x.size() >= 2);
        assert(table.f.size() >= 2);
        assert(x.size() == table.f.size());
        indexer = Indexer1D<T>(table);

        is_log = indexer.IsLog();

        std::vector<T> function_values(table.f.begin(), table.f.end());
        if(is_log) {
            zero_mask.reserve(function_values.size());
            std::transform(function_values.begin(), function_values.end(), zero_mask.begin(), [](T t)->T{return t<=0;});
            std::transform(function_values.begin(), function_values.end(), function_values.begin(), [](T t)->T{if(t > 0){return log(t);} else {return t;}});
        }

        for(unsigned int i=0; i<table.x.size(); ++i) {
            function[
                        xmap[table.x[i]]
            ] = function_values[i];
        }
    }

    T operator()(T const & x) const {

        std::tuple<unsigned int, T, T, T> index_result = indexer(x);
        unsigned int index = std::get<0>(index_result);
        T val = std::get<1>(index_result);
        T xa = std::get<2>(index_result);
        T delta = std::get<3>(index_result);

        if(index < 0) {
            index = 0;
        }
        else if(index >= indexer.NPoints() - 1) {
            index = indexer.NPoints() - 2;
        }

        T fa = function.at(index);
        T fb = function.at(index+1);
        T result;
        if(is_log) {
            if(zero_mask[index]) {
                if(zero_mask[index+1]) {
                    result = fa + (fb - fa) * exp(val - xa - delta);
                } else {
                    result = fa + (exp(fb) - fa) * exp(val - xa - delta);
                }
            } else {
                if(zero_mask[index+1]) {
                    result = exp(fa) + (fb - exp(fa)) * exp(val - xa - delta);
                } else {
                    result = exp(fa + (fb - fa) * (val - xa) / delta);
                }
            }
        } else {
            result = fa + (fb - fa) * (val - xa) / delta;
        }

        if(result < 0) {
            result = 0;
        }

        return result;
    }

    T MinX() const {
        return indexer.Min();
    }

    T RangeX() const {
        return indexer.Range();
    }

    T MaxX() const {
        return indexer.Max();
    }

    bool IsLog() const {
        return is_log;
    }
};


template<typename T>
struct Interpolator2D {
private:
    TableData2D<T> original_table;
    Indexer1D<T> indexer_x;
    Indexer1D<T> indexer_y;
    std::map<std::pair<unsigned int, unsigned int>, bool> zero_mask;
    std::map<std::pair<unsigned int, unsigned int>, T> function;
    bool is_log = false;
public:

    Interpolator2D() {};
    Interpolator2D(TableData2D<T> & table) {
        SetTable(table);
    };

    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(cereal::make_nvp("Table", original_table));
        } else {
            throw std::runtime_error("Interpolator2D only supports version <= 0!");
        }
    };
    template<typename Archive>
    void load(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            TableData2D<T> table;
            archive(cereal::make_nvp("Table", table));
            SetTable(table);
        } else {
            throw std::runtime_error("Interpolator2D only supports version <= 0!");
        }
    };

    void SetTable(TableData2D<T> & table) {
        original_table = table;
        std::set<T> x(table.x.begin(), table.x.end());
        std::set<T> y(table.y.begin(), table.y.end());
        std::map<T, unsigned int> xmap;
        std::map<T, unsigned int> ymap;
        unsigned int n = 0;
        for(auto i : x) {
            xmap[i] = n;
            ++n;
        }
        n = 0;
        for(auto i : y) {
            ymap[i] = n;
            ++n;
        }

        TableData1D<T> x_data;
        TableData1D<T> y_data;

        assert(table.x.size() >= 2);
        assert(table.y.size() >= 2);
        assert(table.f.size() >= 2);

        x_data.x = table.x;
        x_data.f = table.f;
        y_data.x = table.y;
        y_data.f = table.f;

        assert(x_data.x.size() >= 2);
        assert(x_data.f.size() >= 2);
        assert(y_data.x.size() >= 2);
        assert(y_data.f.size() >= 2);

        indexer_x = Indexer1D<T>(x_data);
        indexer_y = Indexer1D<T>(y_data);

        is_log = indexer_x.IsLog() or indexer_y.IsLog();

        std::vector<T> function_values(table.f.begin(), table.f.end());
        std::vector<bool> z_mask;
        if(is_log) {
            z_mask.reserve(function_values.size());
            std::transform(function_values.begin(), function_values.end(), z_mask.begin(), [](T t)->T{return t<=0;});
            std::transform(function_values.begin(), function_values.end(), function_values.begin(), [](T t)->T{if(t > 0){return log(t);} else {return t;}});
            //std::transform(function_values.begin(), function_values.end(), function_values.begin(), [](T t)->T{return log(t);});
        }

        for(unsigned int i=0; i<table.x.size(); ++i) {
            function[
                std::pair<unsigned int, unsigned int>(
                        xmap[table.x[i]],
                        ymap[table.y[i]]
                        )
            ] = function_values[i];
        }
        if(is_log) {
            for(unsigned int i=0; i<table.x.size(); ++i) {
                zero_mask[
                    std::pair<unsigned int, unsigned int>(
                            xmap[table.x[i]],
                            ymap[table.y[i]]
                            )
                ] = z_mask[i];
            }
        }
    }

    T operator()(T const & x, T const & y) const {
        std::tuple<unsigned int, T, T, T> index_result_x = indexer_x(x);
        unsigned int index_x = std::get<0>(index_result_x);
        T val_x = std::get<1>(index_result_x);
        T xa = std::get<2>(index_result_x);
        T delta_x = std::get<3>(index_result_x);

        std::tuple<unsigned int, T, T, T> index_result_y = indexer_y(y);
        unsigned int index_y = std::get<0>(index_result_y);
        T val_y = std::get<1>(index_result_y);
        T ya = std::get<2>(index_result_y);
        T delta_y = std::get<3>(index_result_y);

        if(index_x < 0) {
            index_x = 0;
        }
        else if(index_x >= indexer_x.NPoints() - 1) {
            index_x = indexer_x.NPoints() - 2;
        }

        if(index_y < 0) {
            index_y = 0;
        }
        else if(index_y >= indexer_y.NPoints() - 1) {
            index_y = indexer_y.NPoints() - 2;
        }

        T da_x = val_x - xa;
        T db_x = delta_x - da_x;
        T da_y = val_y - ya;
        T db_y = delta_y - da_y;

        T faa = function.at(std::pair<unsigned int, unsigned int>(index_x, index_y));
        T fab = function.at(std::pair<unsigned int, unsigned int>(index_x, index_y+1));
        T fba = function.at(std::pair<unsigned int, unsigned int>(index_x+1, index_y));
        T fbb = function.at(std::pair<unsigned int, unsigned int>(index_x+1, index_y+1));

        T result;

        if(is_log) {
            bool zaa = zero_mask.at(std::pair<unsigned int, unsigned int>(index_x, index_y));
            bool zab = zero_mask.at(std::pair<unsigned int, unsigned int>(index_x, index_y+1));
            bool zba = zero_mask.at(std::pair<unsigned int, unsigned int>(index_x+1, index_y));
            bool zbb = zero_mask.at(std::pair<unsigned int, unsigned int>(index_x+1, index_y+1));

            bool has_zero = zaa or zab or zba or zbb;
            if(has_zero) {
                if(indexer_x.IsLog()) {
                    da_x = exp(da_x);
                    db_x = exp(db_x);
                    delta_x = exp(delta_x);
                }
                if(indexer_y.IsLog()) {
                    da_y = exp(da_y);
                    db_y = exp(db_y);
                    delta_y = exp(delta_y);
                }
                if(!zaa)
                    faa = exp(faa);
                if(!zab)
                    fab = exp(fab);
                if(!zba)
                    fba = exp(fba);
                if(!zbb)
                    fbb = exp(fbb);
                result = (
                        db_x * db_y * faa +
                        db_x * da_y * fab +
                        da_x * db_y * fba +
                        da_x * da_y * fbb
                        ) / (delta_x * delta_y);
            } else {
                result = exp((
                        db_x * db_y * faa +
                        db_x * da_y * fab +
                        da_x * db_y * fba +
                        da_x * da_y * fbb
                        ) / (delta_x * delta_y));
            }
        } else {
            result = (
                    db_x * db_y * faa +
                    db_x * da_y * fab +
                    da_x * db_y * fba +
                    da_x * da_y * fbb
                    ) / (delta_x * delta_y);
        }

        if(result < 0) {
            result = 0;
        }

        return result;
    }

    T MinX() const {
        return indexer_x.Min();
    }

    T RangeX() const {
        return indexer_x.Range();
    }

    T MaxX() const {
        return indexer_x.Max();
    }

    T MinY() const {
        return indexer_y.Min();
    }

    T RangeY() const {
        return indexer_y.Range();
    }

    T MaxY() const {
        return indexer_y.Max();
    }
};


class DipoleFromTable : public CrossSection {
friend cereal::access;
protected:
DipoleFromTable() {};
public:
    enum HelicityChannel {Conserving, Flipping};
private:
    bool z_samp = true;
    std::map<Particle::ParticleType, Interpolator2D<double>> differential;
    std::map<Particle::ParticleType, Interpolator1D<double>> total;
    const std::set<Particle::ParticleType> primary_types = {Particle::ParticleType::NuE, Particle::ParticleType::NuMu, Particle::ParticleType::NuTau, Particle::ParticleType::NuEBar, Particle::ParticleType::NuMuBar, Particle::ParticleType::NuTauBar};
    double hnl_mass;
    HelicityChannel channel;
public:
    double GetHNLMass() const {return hnl_mass;};
    static double DipoleyMin(double Enu, double mHNL, double target_mass);
    static double DipoleyMax(double Enu, double mHNL, double target_mass);
    DipoleFromTable(double hnl_mass, HelicityChannel channel) : hnl_mass(hnl_mass), channel(channel) {};
    DipoleFromTable(double hnl_mass, HelicityChannel channel, std::set<Particle::ParticleType> const & primary_types) : hnl_mass(hnl_mass), channel(channel), primary_types(primary_types) {};
    double TotalCrossSection(InteractionRecord const &) const;
    double TotalCrossSection(LeptonInjector::Particle::ParticleType primary, double energy, Particle::ParticleType target) const;
    double DifferentialCrossSection(InteractionRecord const &) const;
    double DifferentialCrossSection(Particle::ParticleType primary_type, double primary_energy, Particle::ParticleType target_type, double y) const;
    double InteractionThreshold(InteractionRecord const &) const;
    void SampleFinalState(InteractionRecord &, std::shared_ptr<LeptonInjector::LI_random>) const;

    std::vector<Particle::ParticleType> GetPossibleTargets() const;
    std::vector<Particle::ParticleType> GetPossibleTargetsFromPrimary(Particle::ParticleType primary_type) const;
    std::vector<Particle::ParticleType> GetPossiblePrimaries() const;
    std::vector<InteractionSignature> GetPossibleSignatures() const;
    std::vector<InteractionSignature> GetPossibleSignaturesFromParents(Particle::ParticleType primary_type, Particle::ParticleType target_type) const;

    virtual double FinalStateProbability(InteractionRecord const & record) const;

    void AddDifferentialCrossSectionFile(std::string filename, Particle::ParticleType target);
    void AddTotalCrossSectionFile(std::string filename, Particle::ParticleType target);
    void AddDifferentialCrossSection(Particle::ParticleType target, Interpolator2D<double>);
    void AddTotalCrossSection(Particle::ParticleType target, Interpolator1D<double>);
public:
    virtual std::vector<std::string> DensityVariables() const override;
    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(::cereal::make_nvp("DifferentialCrossSection", differential));
            archive(::cereal::make_nvp("TotalCrossSection", total));
            archive(::cereal::make_nvp("PrimaryTypes", primary_types));
            archive(::cereal::make_nvp("HNLMass", hnl_mass));
            archive(::cereal::make_nvp("HelicityChannel", static_cast<int>(channel)));
            archive(cereal::virtual_base_class<CrossSection>(this));
        } else {
            throw std::runtime_error("DipoleFromTable only supports version <= 0!");
        }
    }
    template<typename Archive>
    void load(Archive & archive, std::uint32_t version) {
        if(version == 0) {
            archive(::cereal::make_nvp("DifferentialCrossSection", differential));
            archive(::cereal::make_nvp("TotalCrossSection", total));
            std::set<LeptonInjector::Particle::ParticleType> prim;
            archive(::cereal::make_nvp("PrimaryTypes", prim));
            archive(::cereal::make_nvp("HNLMass", hnl_mass));
            archive(::cereal::make_nvp("HelicityChannel", channel));
            archive(cereal::virtual_base_class<CrossSection>(this));
        } else {
            throw std::runtime_error("DipoleFromTable only supports version <= 0!");
        }
    }
};

} // namespace LeptonInjector

CEREAL_CLASS_VERSION(LeptonInjector::CrossSection, 0);

CEREAL_CLASS_VERSION(LeptonInjector::DISFromSpline, 0);
CEREAL_REGISTER_TYPE(LeptonInjector::DISFromSpline);
CEREAL_REGISTER_POLYMORPHIC_RELATION(LeptonInjector::CrossSection, LeptonInjector::DISFromSpline);

CEREAL_CLASS_VERSION(LeptonInjector::TableData1D<double>, 0);
CEREAL_CLASS_VERSION(LeptonInjector::TableData2D<double>, 0);
CEREAL_CLASS_VERSION(LeptonInjector::IndexFinderRegular<double>, 0);
CEREAL_CLASS_VERSION(LeptonInjector::IndexFinderIrregular<double>, 0);
CEREAL_CLASS_VERSION(LeptonInjector::Indexer1D<double>, 0);
CEREAL_CLASS_VERSION(LeptonInjector::Interpolator1D<double>, 0);
CEREAL_CLASS_VERSION(LeptonInjector::Interpolator2D<double>, 0);

CEREAL_CLASS_VERSION(LeptonInjector::DipoleFromTable, 0);
CEREAL_REGISTER_TYPE(LeptonInjector::DipoleFromTable);
CEREAL_REGISTER_POLYMORPHIC_RELATION(LeptonInjector::CrossSection, LeptonInjector::DipoleFromTable);

CEREAL_CLASS_VERSION(LeptonInjector::InteractionSignature, 0);
CEREAL_CLASS_VERSION(LeptonInjector::InteractionRecord, 0);
CEREAL_CLASS_VERSION(LeptonInjector::CrossSectionCollection, 0);


#endif // LI_CrossSection_H

