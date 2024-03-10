#pragma once
#ifndef SIREN_Interpolator_H
#define SIREN_Interpolator_H

#include <map>
#include <set>
#include <cmath>
#include <tuple>
#include <limits>
#include <vector>
#include <cassert>
#include <cstdint>
#include <stdlib.h>
#include <algorithm>
#include <stdexcept>

#include <cereal/cereal.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/set.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>


namespace siren {
namespace utilities {

template<typename T>
struct TableData1D {
    std::vector<T> x;
    std::vector<T> f;

    template<typename U>
    bool operator==(TableData1D<U> const & other) const {
        return std::tie(x, f) == std::tie(other.x, other.f);
    }

    template<typename U>
    bool operator<(TableData1D<U> const & other) const {
        return std::tie(x, f) < std::tie(other.x, other.f);
    }

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

    template<typename U>
    bool operator==(TableData2D<U> const & other) const {
        return std::tie(x, y, f) == std::tie(other.x, other.y, other.f);
    }

    template<typename U>
    bool operator<(TableData2D<U> const & other) const {
        return std::tie(x, y, f) < std::tie(other.x, other.y, other.f);
    }

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
        else if(i >= int(n_points) - 1)
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

    template<typename U>
    bool operator==(Interpolator1D<U> const & other) const {
        return original_table == other.original_table;
    }

    template<typename U>
    bool operator<(Interpolator1D<U> const & other) const {
        return original_table < other.original_table;
    }

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

    template<typename U>
    bool operator==(Interpolator2D<U> const & other) const {
        return original_table == other.original_table;
    }

    template<typename U>
    bool operator<(Interpolator2D<U> const & other) const {
        return original_table < other.original_table;
    }

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

} // namespace utilities
} // namespace siren

CEREAL_CLASS_VERSION(siren::utilities::TableData1D<double>, 0);
CEREAL_CLASS_VERSION(siren::utilities::TableData2D<double>, 0);
CEREAL_CLASS_VERSION(siren::utilities::IndexFinderRegular<double>, 0);
CEREAL_CLASS_VERSION(siren::utilities::IndexFinderIrregular<double>, 0);
CEREAL_CLASS_VERSION(siren::utilities::Indexer1D<double>, 0);
CEREAL_CLASS_VERSION(siren::utilities::Interpolator1D<double>, 0);
CEREAL_CLASS_VERSION(siren::utilities::Interpolator2D<double>, 0);

#endif // SIREN_Interpolator_H
