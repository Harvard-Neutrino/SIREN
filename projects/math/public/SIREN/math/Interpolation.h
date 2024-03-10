#ifndef SIREN_Interpolation_H
#define SIREN_Interpolation_H

#include <map>
#include <set>
#include <array>
#include <cmath>
#include <tuple>
#include <limits>
#include <memory>
#include <vector>
#include <cstdint>
#include <cstdlib>
#include <utility>
#include <iostream>
#include <iterator>
#include <typeinfo>
#include <algorithm>
#include <typeindex>
#include <stdexcept>
#include <functional>

#include <cereal/cereal.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/set.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>

#include "delabella.h"

namespace siren {
namespace math {

template<typename T>
class Transform {
public:
    virtual T Function(T x) const = 0;
    virtual T Inverse(T x) const = 0;
    virtual ~Transform() {}
    virtual bool equal(Transform<T> const & other) const = 0;
    virtual bool less(Transform<T> const & other) const = 0;
    virtual bool operator==(Transform<T> const & other) const {
		if(this == &other)
			return true;
		else
			return this->equal(other);
    }
	virtual bool operator<(Transform<T> const & other) const {
		if(typeid(this) == typeid(&other))
			return this->less(other);
		else
			return std::type_index(typeid(this)) < std::type_index(typeid(&other));
	}
    template<class Archive>
    void serialize(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
        } else {
            throw std::runtime_error("Transform only supports version <= 0!");
        }
    }
};

template<typename T>
class IdentityTransform : public Transform<T> {
public:
    IdentityTransform() {}
    virtual T Function(T x) const override {
        return x;
    }
    virtual T Inverse(T x) const override {
        return x;
    }
    virtual bool equal(Transform<T> const & other) const override {
        const IdentityTransform* x = dynamic_cast<const IdentityTransform*>(&other);
        if(!x)
            return false;
        return true;
    }
    virtual bool less(Transform<T> const & other) const override {
        return false;
    }
    template<class Archive>
    void serialize(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(cereal::virtual_base_class<Transform<T>>(this));
        } else {
            throw std::runtime_error("IdentityTransform only supports version <= 0!");
        }
    }
};

template<typename T>
class LogTransform : public Transform<T> {
public:
    LogTransform() {}
    virtual T Function(T x) const override {
        return log(x);
    }
    virtual T Inverse(T x) const override {
        return exp(x);
    }
    virtual bool equal(Transform<T> const & other) const override {
        const LogTransform* x = dynamic_cast<const LogTransform*>(&other);
        if(!x)
            return false;
        return true;
    }
    virtual bool less(Transform<T> const & other) const override {
        return false;
    }
    template<class Archive>
    void serialize(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(cereal::virtual_base_class<Transform<T>>(this));
        } else {
            throw std::runtime_error("LogTransform only supports version <= 0!");
        }
    }
};

template<typename T>
class SymLogTransform : public Transform<T> {
    T min_x;
    T log_min_x;
public:
    SymLogTransform(T min_x) : min_x(std::abs(min_x)), log_min_x(log(std::abs(min_x))) {
        if(min_x == 0) {
            throw std::runtime_error("SymLogTransform cannot be initialized with a minimum value of x=0");
        }
    }
    SymLogTransform(std::vector<T> const & x) {
        bool first = true;
        min_x = 0;
        for(size_t i=0; i<x.size(); ++i) {
            if(x[i] != 0) {
                if(first) {
                    first = false;
                    min_x = std::abs(x[i]);
                }
                min_x = std::min(min_x, std::abs(x[i]));
            }
        }
        log_min_x = log(min_x);
        if(min_x == 0) {
            throw std::runtime_error("SymLogTransform cannot be initialized with a minimum value of x=0");
        }
    }
    virtual T Function(T x) const override {
        T abs_x = std::abs(x);
        if(abs_x < min_x) {
            return x;
        } else {
            return std::copysign(log(abs_x) - log_min_x + min_x, x);
        }
    }
    virtual T Inverse(T x) const override {
        T abs_x = std::abs(x);
        if(abs_x < min_x) {
            return x;
        } else {
            return std::copysign(exp(abs_x - min_x + log_min_x), x);
        }
    }
    virtual bool equal(Transform<T> const & other) const override {
        const SymLogTransform* x = dynamic_cast<const SymLogTransform*>(&other);
        if(!x)
            return false;
        return true;
    }
    virtual bool less(Transform<T> const & other) const override {
        const SymLogTransform* x = dynamic_cast<const SymLogTransform*>(&other);
        return min_x < x->min_x;
    }
    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(::cereal::make_nvp("MinX", min_x));
            archive(cereal::virtual_base_class<Transform<T>>(this));
        } else {
            throw std::runtime_error("SymLogTransform only supports version <= 0!");
        }
    }
    template<typename Archive>
    static void load_and_construct(Archive & archive, cereal::construct<SymLogTransform> & construct, std::uint32_t const version) {
        if(version == 0) {
			T min_x;
            archive(::cereal::make_nvp("MinX", min_x));
            construct(min_x);
            archive(cereal::virtual_base_class<Transform<T>>(construct.ptr()));
        } else {
            throw std::runtime_error("SymLogTransform only supports version <= 0!");
        }
    }
};

template<typename T>
class RangeTransform : public Transform<T> {
    T min_x;
    T range;
public:
    RangeTransform(T min_x, T max_x) : min_x(min_x), range(max_x - min_x) {
        if(range == 0) {
            throw std::runtime_error("RangeTransform cannot be initialized with a range of zero");
        }
    }
    virtual T Function(T x) const override {
        return (x - min_x) / range;
    }
    virtual T Inverse(T x) const override {
        return x * range + min_x;
    }
    virtual bool equal(Transform<T> const & other) const override {
        const RangeTransform* x = dynamic_cast<const RangeTransform*>(&other);
        if(!x)
            return false;
        return true;
    }
    virtual bool less(Transform<T> const & other) const override {
        const RangeTransform* x = dynamic_cast<const RangeTransform*>(&other);
        return std::tie(min_x, range) < std::tie(x->min_x, x->range);
    }
    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(::cereal::make_nvp("MinX", min_x));
            archive(::cereal::make_nvp("MaxX", min_x + range));
            archive(cereal::virtual_base_class<Transform<T>>(this));
        } else {
            throw std::runtime_error("RangeTransform only supports version <= 0!");
        }
    }
    template<typename Archive>
    static void load_and_construct(Archive & archive, cereal::construct<RangeTransform> & construct, std::uint32_t const version) {
        if(version == 0) {
			T min_x;
			T max_x;
            archive(::cereal::make_nvp("MinX", min_x));
            archive(::cereal::make_nvp("MaxX", max_x));
            construct(min_x, max_x);
            archive(cereal::virtual_base_class<Transform<T>>(construct.ptr()));
        } else {
            throw std::runtime_error("RangeTransform only supports version <= 0!");
        }
    }
};

template<typename T>
struct LinearInterpolationOperator {
public:
    LinearInterpolationOperator() {};
    virtual T operator()(T const & x0, T const & y0, T const & x1, T const & y1, T const & x) const {
        T delta_x = x1 - x0;
        T delta_y = y1 - y0;
        return (x - x0) * delta_y / delta_x + y0;
    }
    virtual bool equal(LinearInterpolationOperator<T> const & other) const {
        return true;
    }
    virtual bool less(LinearInterpolationOperator<T> const & other) const {
        return false;
    }
    virtual bool operator==(LinearInterpolationOperator<T> const & other) const {
		if(this == &other)
			return true;
		else
			return this->equal(other);
    }
	virtual bool operator<(LinearInterpolationOperator<T> const & other) const {
		if(typeid(this) == typeid(&other))
			return this->less(other);
		else
			return std::type_index(typeid(this)) < std::type_index(typeid(&other));
	}
    template<class Archive>
    void serialize(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
        } else {
            throw std::runtime_error("LinearInterpolationOperator only supports version <= 0!");
        }
    }
    virtual ~LinearInterpolationOperator() {}
};

template<typename T>
struct DropLinearInterpolationOperator : public LinearInterpolationOperator<T> {
public:
    DropLinearInterpolationOperator() {}
    virtual T operator()(T const & x0, T const & y0, T const & x1, T const & y1, T const & x) const override {
        if(y0 == 0 or y1 == 0)
            return 0;
        return LinearInterpolationOperator<T>::operator()(x0, y0, x1, y1, x);
    }
    virtual bool equal(LinearInterpolationOperator<T> const & other) const override {
        const DropLinearInterpolationOperator* x = dynamic_cast<const DropLinearInterpolationOperator*>(&other);
        if(!x)
            return false;
        return true;
    }
    virtual bool less(LinearInterpolationOperator<T> const & other) const override {
        return false;
    }
    template<class Archive>
    void serialize(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(cereal::virtual_base_class<LinearInterpolationOperator<T>>(this));
        } else {
            throw std::runtime_error("DropLinearInterpolationOperator only supports version <= 0!");
        }
    }
};

template<typename T>
std::tuple<std::shared_ptr<Transform<T>>, std::shared_ptr<Transform<T>>> SelectInterpolationSpace1D(
        std::vector<T> const & x,
        std::vector<T> const & y,
        std::shared_ptr<LinearInterpolationOperator<T>> interp) {
    T min_x = 1;
    T min_y = 1;
    bool first_x = true;
    bool first_y = false;
    for(size_t i=0; i<x.size(); ++i) {
        if(x[i] != 0) {
            if(first_x) {
                min_x = std::abs(x[i]);
                first_x = false;
            }
            min_x = std::min(min_x, std::abs(x[i]));
        }
        if(y[i] != 0) {
            if(first_y) {
                min_y = std::abs(y[i]);
                first_y = false;
            }
            min_y = std::min(min_y, std::abs(y[i]));
        }
    }
    std::vector<T> symlog_x(x.size());
    std::vector<T> symlog_y(y.size());
    SymLogTransform<T> symlog_t_x(min_x);
    SymLogTransform<T> symlog_t_y(min_y);
    for(size_t i=0; i<x.size(); ++i) {
        symlog_x[i] = symlog_t_x.Function(x[i]);
        symlog_y[i] = symlog_t_y.Function(y[i]);
    }

    T tot_00 = 0;
    T tot_01 = 0;
    T tot_10 = 0;
    T tot_11 = 0;

    for(size_t i=1; i<x.size()-1; ++i) {
        T y_i_estimate = interp->operator()(x[i-1], y[i-1], x[i+1], y[i+1], x[i]);
        T den_00 = std::max(std::abs(y[i]), std::abs(y_i_estimate)); if(den_00 < 1e-8) den_00 = 1;
        T delta_00 = (y_i_estimate - y[i]) / den_00;
        if((not std::isnan(delta_00)) and (not std::isinf(delta_00)))
            tot_00 += delta_00 * delta_00;

        y_i_estimate = interp->operator()(x[i-1], symlog_y[i-1], x[i+1], symlog_y[i+1], x[i]);
        T den_01 = std::max(std::abs(y[i]), std::abs(symlog_t_y.Inverse(y_i_estimate))); if(den_01 < 1e-8) den_01 = 1;
        T delta_01 = (symlog_t_y.Inverse(y_i_estimate) - y[i]) / den_01;
        if((not std::isnan(delta_01)) and (not std::isinf(delta_01)))
            tot_01 += delta_01 * delta_01;

        y_i_estimate = interp->operator()(symlog_x[i-1], y[i-1], symlog_x[i+1], y[i+1], symlog_x[i]);
        T den_10 = std::max(std::abs(y[i]), std::abs(y_i_estimate)); if(den_10 < 1e-8) den_10 = 1;
        T delta_10 = (y_i_estimate - y[i]) / den_10;
        if((not std::isnan(delta_10)) and (not std::isinf(delta_10)))
            tot_10 += delta_10 * delta_10;

        y_i_estimate = interp->operator()(symlog_x[i-1], symlog_y[i-1], symlog_x[i+1], symlog_y[i+1], symlog_x[i]);
        T den_11 = std::max(std::abs(y[i]), std::abs(symlog_t_y.Inverse(y_i_estimate))); if(den_11 < 1e-8) den_11 = 1;
        T delta_11 = (symlog_t_y.Inverse(y_i_estimate) - y[i]) / den_11;
        if((not std::isnan(delta_11)) and (not std::isinf(delta_11)))
            tot_11 += delta_11 * delta_11;
    }

    std::vector<T> tots = {tot_00, tot_01, tot_10, tot_11};
    int min_combination = std::distance(tots.begin(), std::min_element(tots.begin(), tots.end()));

    switch(min_combination) {
        case 0:
            return {std::make_shared<IdentityTransform<T>>(), std::make_shared<IdentityTransform<T>>()};
            break;
        case 1:
            return {std::make_shared<IdentityTransform<T>>(), std::make_shared<SymLogTransform<T>>(min_y)};
            break;
        case 2:
            return {std::make_shared<SymLogTransform<T>>(min_x), std::make_shared<IdentityTransform<T>>()};
            break;
        case 3:
            return {std::make_shared<SymLogTransform<T>>(min_x), std::make_shared<SymLogTransform<T>>(min_y)};
            break;
        default:
            return {std::make_shared<IdentityTransform<T>>(), std::make_shared<IdentityTransform<T>>()};
            break;
    }
};

template<typename T>
class Indexer1D {
public:
    virtual std::tuple<int, int> operator()(T const & x) const = 0;
    virtual bool equal(Indexer1D<T> const & other) const = 0;
    virtual bool less(Indexer1D<T> const & other) const = 0;
    virtual bool operator==(Indexer1D<T> const & other) const {
		if(this == &other)
			return true;
		else
			return this->equal(other);
    }
	virtual bool operator<(Indexer1D<T> const & other) const {
		if(typeid(this) == typeid(&other))
			return this->less(other);
		else
			return std::type_index(typeid(this)) < std::type_index(typeid(&other));
	}
    template<class Archive>
    void serialize(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
        } else {
            throw std::runtime_error("Indexer1D only supports version <= 0!");
        }
    }
    virtual ~Indexer1D() {}
};

template<typename T>
class RegularIndexer1D : public Indexer1D<T> {
    T low;
    T high;
    T range;
    bool reversed;
    unsigned int n_points;
    T delta;
public:

    RegularIndexer1D() {};
    RegularIndexer1D(std::vector<T> data) {
        if(data.size() < 2)
            throw std::runtime_error("RegularIndexer1D needs at least 2 points");
        n_points = data.size();
        std::vector<T> unsorted = data;
        std::sort(data.begin(), data.end());
        std::vector<T> reversed_data(data.size());
        std::reverse_copy(data.begin(), data.end(), reversed_data.begin());
        if(data == unsorted)  {
            // Assume sorted to start with
            reversed = false;
        } else if(data == reversed_data) {
            // Assume reversed to start with
            reversed = true;
        } else {
            // Unsorted
            reversed = false;
        }
        low = data.front();
        high = data.back();
        range = high - low;
        if(std::abs(range) <= 0)
            throw std::runtime_error("RegularIndexer1D cannot function with zero range");
        delta = range / (n_points - 1);
    };

    virtual bool equal(Indexer1D<T> const & other) const override {
        const RegularIndexer1D* x = dynamic_cast<const RegularIndexer1D*>(&other);
        if(!x)
            return false;
        return std::tie(low, high, range, reversed, n_points, delta) == std::tie(x->low, x->high, x->range, x->reversed, x->n_points, x->delta);
    }
    virtual bool less(Indexer1D<T> const & other) const override {
        const RegularIndexer1D* x = dynamic_cast<const RegularIndexer1D*>(&other);
        return std::tie(low, high, range, reversed, n_points, delta) < std::tie(x->low, x->high, x->range, x->reversed, x->n_points, x->delta);
    }
    template<class Archive>
    void serialize(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(cereal::make_nvp("Low", low));
            archive(cereal::make_nvp("High", high));
            archive(cereal::make_nvp("Range", range));
            archive(cereal::make_nvp("Reversed", reversed));
            archive(cereal::make_nvp("NPoints", n_points));
            archive(cereal::make_nvp("Delta", delta));
            archive(cereal::virtual_base_class<Indexer1D<T>>(this));
        } else {
            throw std::runtime_error("RegularIndexer1D only supports version <= 0!");
        }
    }

    virtual std::tuple<int, int> operator()(T const & x) const override {
        int index = (int)std::floor((x - low) / range * (n_points - 1));
        if(reversed)
            index = (n_points - 1) - index;
        if(index < 0)
            index = 0;
        else if(index >= int(n_points) - 1)
            index = n_points - 2;
        return {index, index+1};
    }
};

template<typename T>
class IrregularIndexer1D : public Indexer1D<T> {
    std::vector<T> data;
    T low;
    T high;
    bool reversed;
    unsigned int n_points;
public:
    IrregularIndexer1D() {};
    IrregularIndexer1D(std::vector<T> x): data(x.begin(), x.end()) {
        if(x.size() < 2)
            throw std::runtime_error("IrregularIndexer1D needs at least 2 points");
        std::vector<T> unsorted = data;
        std::sort(data.begin(), data.end());
        std::vector<T> reversed_data(data.size());
        std::reverse_copy(data.begin(), data.end(), reversed_data.begin());
        if(data == unsorted)  {
            // Assume sorted to start with
            reversed = false;
        } else if(data == reversed_data) {
            // Assume reversed to start with
            reversed = true;
        } else {
            // Unsorted
            reversed = false;
        }
        low = data.front();
        high = data.back();
        T range = high - low;
        n_points = data.size();
        if(std::abs(range) <= 0)
            throw std::runtime_error("IrregularIndexer1D cannot function with zero range");
    };

    virtual bool equal(Indexer1D<T> const & other) const override {
        const IrregularIndexer1D* x = dynamic_cast<const IrregularIndexer1D*>(&other);
        if(!x)
            return false;
        return std::tie(data, low, high, reversed, n_points) == std::tie(x->data, x->low, x->high, x->reversed, x->n_points);
    }
    virtual bool less(Indexer1D<T> const & other) const override {
        const IrregularIndexer1D* x = dynamic_cast<const IrregularIndexer1D*>(&other);
        return std::tie(data, low, high, reversed, n_points) < std::tie(x->data, x->low, x->high, x->reversed, x->n_points);
    }
    template<class Archive>
    void serialize(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(cereal::make_nvp("Data", data));
            archive(cereal::make_nvp("Low", low));
            archive(cereal::make_nvp("High", high));
            archive(cereal::make_nvp("Reversed", reversed));
            archive(cereal::make_nvp("NPoints", n_points));
            archive(cereal::virtual_base_class<Indexer1D<T>>(this));
        } else {
            throw std::runtime_error("IrregularIndexer1D only supports version <= 0!");
        }
    }

    std::tuple<int, int> operator()(T const & x) const override {
        if(x <= low) {
            return {0, 1};
        } else if (x >= high) {
            return {n_points - 2, n_points - 1};
        }
        // Lower bound returns pointer to element that is greater than or equal to x
        // i.e. x \in (a,b] --> pointer to b, x \in (b,c] --> pointer to c
        // begin is the first element
        // distance(begin, pointer to y) --> y
        // therefore this function returns the index of the lower bin edge
        using ConstIterator = typename std::vector<T>::const_iterator;
        ConstIterator it = std::lower_bound(data.begin(), data.end(), x);
        unsigned int index = std::distance(data.begin(), it) - 1;
        if(reversed)
            index = (n_points - 1) - index;
        if(index < 0)
            index = 0;
        else if(index >= n_points - 1)
            index = n_points - 2;
        return {index, index + 1};
    }
};

template<typename T>
class TransformIndexer1D : public Indexer1D<T> {
    std::shared_ptr<Indexer1D<T>> indexer;
    std::shared_ptr<Transform<T>> transform;
public:
    TransformIndexer1D() {}
    TransformIndexer1D(std::shared_ptr<Indexer1D<T>> indexer, std::shared_ptr<Transform<T>> transform) {
        this->indexer = indexer;
        this->transform = transform;
    }

    virtual bool equal(Indexer1D<T> const & other) const override {
        const TransformIndexer1D* x = dynamic_cast<const TransformIndexer1D*>(&other);
        if(!x)
            return false;
        return std::tie(*indexer, *transform) == std::tie(*(x->indexer), *(x->transform));
    }
    virtual bool less(Indexer1D<T> const & other) const override {
        const TransformIndexer1D* x = dynamic_cast<const TransformIndexer1D*>(&other);
        return std::tie(*indexer, *transform) < std::tie(*(x->indexer), *(x->transform));
    }
    template<class Archive>
    void serialize(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(cereal::make_nvp("Indexer", indexer));
            archive(cereal::make_nvp("Transform", transform));
            archive(cereal::virtual_base_class<Indexer1D<T>>(this));
        } else {
            throw std::runtime_error("TransformIndexer1D only supports version <= 0!");
        }
    }

    std::tuple<int, int> operator()(T const & x) const override {
        T t = transform->Function(x);
        return indexer->operator()(t);
    }
};

template<typename T>
std::shared_ptr<Indexer1D<T>> SelectIndexer1D(
        std::vector<T> data,
        std::shared_ptr<Transform<T>> x_transform) {
    std::vector<T> x = data;
    std::sort(x.begin(), x.end());

    // Compute and compare the squared error between different indexing methods and select the most accurate

    T metric_x_irregular = 0;
    T metric_t_irregular = 0;
    T metric_symlog_irregular = 0;
    for(size_t i=0; i<x.size() - 1; ++i) {
        metric_x_irregular += std::pow(std::numeric_limits<T>::epsilon() * std::abs(x[i]), 2);
    }
    metric_x_irregular -= metric_x_irregular * std::numeric_limits<T>::epsilon() * x.size();
    if(x_transform == nullptr) {
        metric_t_irregular = std::numeric_limits<T>::max();
    } else {
        for(size_t i=0; i<x.size() - 1; ++i) {
            T t = x_transform->Function(x[i]);
            T x_p = x_transform->Inverse(t);
            metric_t_irregular += std::pow(std::max(std::abs(x[i] - x_p), std::numeric_limits<T>::epsilon() * std::abs(x[i])), 2);
        }
        metric_t_irregular -= metric_t_irregular * std::numeric_limits<T>::epsilon() * x.size();
    }
    if(dynamic_cast<SymLogTransform<T>*>(x_transform.get()) != nullptr) {
        metric_symlog_irregular = std::numeric_limits<T>::max();
    } else {
        SymLogTransform<T> symlog(x);
        for(size_t i=0; i<x.size() - 1; ++i) {
            T t = symlog.Function(x[i]);
            T x_p = symlog.Inverse(t);
            metric_symlog_irregular += std::pow(std::max(std::abs(x[i] - x_p), std::numeric_limits<T>::epsilon() * std::abs(x[i])), 2);
        }
        metric_symlog_irregular -= metric_symlog_irregular * std::numeric_limits<T>::epsilon() * x.size();
    }

    T thresh = 1.0 - std::numeric_limits<T>::epsilon() * 100 * x.size();

    T metric_x_regular = 0;
    T metric_t_regular = 0;
    T metric_symlog_regular = 0;
    T min_x = *std::min_element(x.begin(), x.end());
    T max_x = *std::max_element(x.begin(), x.end());
    T delta_x = (max_x - min_x) / (x.size() - 1);
    for(size_t i=0; i<x.size() - 1; ++i) {
        T x_reg = delta_x * i + min_x;
        metric_x_regular += std::pow(std::max(std::numeric_limits<T>::epsilon() * std::abs(x[i])*thresh, std::min(std::abs(x[i] - x_reg), std::numeric_limits<T>::epsilon() * std::abs(x[i]))), 2);
    }
    metric_x_regular /= thresh;

    if(x_transform == nullptr) {
        metric_t_regular = std::numeric_limits<T>::max();
    } else {
        T min_t = x_transform->Function(min_x);
        T max_t = x_transform->Function(max_x);
        T delta_t = (max_t - min_t) / (x.size() - 1);
        for(size_t i=0; i<x.size() - 1; ++i) {
            T t_reg = delta_t * i + min_t;
            T x_p = x_transform->Inverse(t_reg);
            metric_t_regular += std::pow(std::max(std::numeric_limits<T>::epsilon() * std::abs(x[i])*thresh, std::min(std::abs(x[i] - x_p), std::numeric_limits<T>::epsilon() * std::abs(x[i]))), 2);
        }
        metric_t_regular /= thresh;
    }

    if(dynamic_cast<SymLogTransform<T>*>(x_transform.get()) != nullptr) {
        metric_symlog_regular = std::numeric_limits<T>::max();
    } else {
        SymLogTransform<T> symlog(x);
        T min_sym = symlog.Function(min_x);
        T max_sym = symlog.Function(max_x);
        T delta_sym = (max_sym - min_sym) / (x.size() - 1);
        for(size_t i=0; i<x.size() - 1; ++i) {
            T t_reg = delta_sym * i + min_sym;
            T x_p = symlog.Inverse(t_reg);
            metric_symlog_regular += std::pow(std::max(std::numeric_limits<T>::epsilon() * std::abs(x[i])*thresh, std::min(std::abs(x[i] - x_p), std::numeric_limits<T>::epsilon() * std::abs(x[i]))), 2);
        }
        metric_symlog_regular /= thresh;
    }
    std::vector<T> metrics = {metric_x_irregular, metric_x_regular, metric_t_irregular, metric_t_regular, metric_symlog_irregular, metric_symlog_regular};
    size_t min_index = std::distance(metrics.begin(), std::min_element(metrics.begin(), metrics.end()));
    if(min_index <= 1) { // No transformation
        if(min_index == 0) { // Irregular
            return std::shared_ptr<Indexer1D<T>>(new IrregularIndexer1D<T>(data));
        } else { // Regular
            return std::shared_ptr<Indexer1D<T>>(new RegularIndexer1D<T>(data));
        }
    } else {
        std::shared_ptr<Transform<T>> transform;
        std::shared_ptr<Indexer1D<T>> indexer;
        if(min_index <= 3) { // Supplied transform
            transform = x_transform;
        } else { // Symlog transform
            transform = std::shared_ptr<Transform<T>>(new SymLogTransform<T>(data));
        }
        std::vector<T> t(data.size());
        for(size_t i=0; i<data.size(); ++i) {
            t[i] = transform->Function(data[i]);
        }
        if(min_index % 2 == 0) { // Irregular
            indexer = std::shared_ptr<Indexer1D<T>>(new IrregularIndexer1D<T>(t));
        } else {
            indexer = std::shared_ptr<Indexer1D<T>>(new RegularIndexer1D<T>(t));
        }
        return std::shared_ptr<Indexer1D<T>>(new TransformIndexer1D<T>(indexer, transform));
    }
};

template<typename T>
class Interpolator1D {
    virtual T operator()(T x) const = 0;
};

template<typename T>
class LinearInterpolator1D : public Interpolator1D<T> {
    // Data
    std::vector<T> t_x;
    std::vector<T> t_y;

    // Indexer
    std::shared_ptr<Indexer1D<T>> indexer;

    // Interpolation info
    std::shared_ptr<LinearInterpolationOperator<T>> interp;
    std::shared_ptr<Transform<T>> x_transform;
    std::shared_ptr<Transform<T>> y_transform;

    LinearInterpolator1D(std::vector<T> x, std::vector<T> y, std::shared_ptr<Indexer1D<T>> indexer, std::shared_ptr<LinearInterpolationOperator<T>> interp, std::shared_ptr<Transform<T>> x_transform, std::shared_ptr<Transform<T>> y_transform)
        : t_x(), t_y(), indexer(indexer), interp(interp), x_transform(x_transform), y_transform(y_transform) {
        t_x.reserve(x.size());
        t_y.reserve(y.size());
        for(size_t i=0; i<x.size(); ++i) {
            t_x.push_back(x_transform->operator()(x[i]));
            t_y.push_back(y_transform->operator()(y[i]));
        }
    }

    virtual T operator()(T x) const override {
        std::tuple<int, int> indices = indexer->operator()(x);
        T const & t_x0 = t_x[std::get<0>(indices)];
        T const & t_y0 = t_y[std::get<0>(indices)];
        T const & t_x1 = t_x[std::get<1>(indices)];
        T const & t_y1 = t_y[std::get<1>(indices)];

        T t_x = x_transform->Function(x);

        T t_res = interp->operator()(t_x0, t_y0, t_x1, t_y1, t_x);
        T res = y_transform->Inverse(t_res);
        return res;
    }
};

template<typename T>
struct BiLinearInterpolationOperator {
public:
    BiLinearInterpolationOperator() {}
    virtual T operator()(T const & x0, T const & y0, T const & x1, T const & y1, T const & z00, T const & z01, T const & z10, T const & z11, T const & x, T const & y) const {

        T delta_x = x1 - x0;
        T delta_y = y1 - y0;
        T dx0 = x - x0;
        T dx1 = x1 - x;
        T dy0 = y - y0;
        T dy1 = y1 - y;

        T res =
            dx1 * dy1 * z00 +
            dx1 * dy0 * z01 +
            dx0 * dy1 * z10 +
            dx0 * dy0 * z11;

        return res / (delta_x * delta_y);
    }
};

template<typename T>
struct DropBiLinearInterpolationOperator : public BiLinearInterpolationOperator<T> {
    DropBiLinearInterpolationOperator() {}
    virtual T operator()(T const & x0, T const & x1, T const & y0, T const & y1, T const & z00, T const & z01, T const & z10, T const & z11, T const & x, T const & y) const override {
        if((z00 == 0 and (z01 == 0 or z10 == 0)) or (z11 == 0 and (z01 == 0 or z10 == 0))) {
            return 0.0;
        }
        return BiLinearInterpolationOperator<T>::operator()(x0, x1, y0, y1, z00, z01, z10, z11, x, y);
    }
};

template<typename T>
class DelaunayIndexer2D;

template<typename T>
struct SimplexLinearInterpolationOperator {
    static inline T sq(T const & x) {
        return x*x;
    }
    using Tri = typename DelaunayIndexer2D<T>::Tri;
    using Point = typename DelaunayIndexer2D<T>::Point;
public:
    virtual T operator()(T const & x, T const & y, Tri const * tri, T z1, T z2, T z3) const {
        Point const & v1 = tri->p0;
        Point const & v2 = tri->p1;
        Point const & v3 = tri->p2;
        T denominator = (v2.y - v3.y) * (v1.x - v3.x) + (v3.x - v2.x) * (v1.y - v3.y);
        T w1 = ((v2.y - v3.y) * (x - v3.x) + (v3.x - v2.x) * (y - v3.y)) / denominator;
        T w2 = ((v3.y - v1.y) * (x - v3.x) + (v1.x - v3.x) * (y - v3.y)) / denominator;
        T w3 = 1.0 - w1 - w2;
        return w1 * z1 + w2 * z2 + w3 * z3;
    }
    virtual ~SimplexLinearInterpolationOperator() {}
};

template<typename T>
struct DropSimplexLinearInterpolationOperator : public SimplexLinearInterpolationOperator<T> {
    using Tri = typename DelaunayIndexer2D<T>::Tri;
public:
    virtual T operator()(T const & x, T const & y, Tri const * tri, T z1, T z2, T z3) const override {
        if(z1 == 0 or z2 == 0 or z3 == 0)
            return 0.0;
        return SimplexLinearInterpolationOperator<T>::operator()(x, y, tri, z1, z1, z3);
    }
};

struct Index2D {
    int x0;
    int y0;
    int x1;
    int y1;
    int z00;
    int z01;
    int z10;
    int z11;
    Index2D(int x0, int y0, int x1, int y1, int z00, int z01, int z10, int z11)
        : x0(x0), y0(y0), x1(x1), y1(y1), z00(z00), z01(z01), z10(z10), z11(z11) {}
};

template<typename T>
class Indexer2D {
public:
    virtual Index2D operator()(T const & x, T const & y) const = 0;
    virtual ~Indexer2D() = default;
};

template<typename T, typename IndexerT>
class GridIndexer2D : public Indexer2D<T> {
    std::map<std::tuple<int, int>, int> xy_to_z_map;
    IndexerT x_indexer;
    IndexerT y_indexer;
public:
    GridIndexer2D() {};
    GridIndexer2D(std::vector<T> x, std::vector<T> y) {
        std::set<T> unique_x(x.begin(), x.end());
        std::set<T> unique_y(y.begin(), y.end());
        std::vector<T> sorted_x(unique_x.begin(), unique_x.end());
        std::vector<T> sorted_y(unique_y.begin(), unique_y.end());
        std::sort(sorted_x.begin(), sorted_x.end());
        std::sort(sorted_y.begin(), sorted_y.end());
        x_indexer = IndexerT(sorted_x);
        y_indexer = IndexerT(sorted_y);
        std::map<T, size_t> x_map;
        std::map<T, size_t> y_map;
        for(size_t i=0; i<sorted_x.size(); ++i)
            x_map.insert({sorted_x[i], i});
        for(size_t i=0; i<sorted_y.size(); ++i)
            y_map.insert({sorted_y[i], i});
        for(size_t i=0; i<x.size(); ++i)
            xy_to_z_map.insert({{x_map[x[i]], y_map[y[i]]}, i});
    };

    virtual Index2D operator()(T const & x, T const & y) const override {
        std::tuple<int, int> x_indices = x_indexer(x);
        std::tuple<int, int> y_indices = y_indexer(y);
        int x0 = std::get<0>(x_indices);
        int y0 = std::get<0>(y_indices);
        int x1 = std::get<1>(x_indices);
        int y1 = std::get<1>(y_indices);
        std::tuple<int, int> xy00(x0, y0);
        std::tuple<int, int> xy01(x0, y1);
        std::tuple<int, int> xy10(x1, y0);
        std::tuple<int, int> xy11(x1, y1);
        int z00 = xy_to_z_map.at(xy00);
        int z01 = xy_to_z_map.at(xy01);
        int z10 = xy_to_z_map.at(xy10);
        int z11 = xy_to_z_map.at(xy11);
        return Index2D(x0, y0, x1, y1, z00, z01, z10, z11);
    }
};

template<typename T>
using RegularGridIndexer2D = GridIndexer2D<T, RegularIndexer1D<T>>;

template<typename T>
using IrregularGridIndexer2D = GridIndexer2D<T, IrregularIndexer1D<T>>;

namespace quadtree {
// QuadTree implementation from https://github.com/pvigier/Quadtree

template<typename T>
class Vector2 {
public:
    T x;
    T y;

    constexpr Vector2<T>(T X = 0, T Y = 0) noexcept : x(X), y(Y) {

    }

    constexpr Vector2<T>& operator+=(const Vector2<T>& other) noexcept {
        x += other.x;
        y += other.y;
        return *this;
    }

    constexpr Vector2<T>& operator/=(T t) noexcept {
        x /= t;
        y /= t;
        return *this;
    }
};

template<typename T>
constexpr Vector2<T> operator+(Vector2<T> lhs, const Vector2<T>& rhs) noexcept {
    lhs += rhs;
    return lhs;
}

template<typename T>
constexpr Vector2<T> operator/(Vector2<T> vec, T t) noexcept {
    vec /= t;
    return vec;
}

template<typename T>
class Box {
public:
    T left;
    T top;
    T width; // Must be positive
    T height; // Must be positive

    constexpr Box(T Left = 0, T Top = 0, T Width = 0, T Height = 0) noexcept :
        left(Left), top(Top), width(Width), height(Height) {

    }

    constexpr Box(const Vector2<T>& position, const Vector2<T>& size) noexcept :
        left(position.x), top(position.y), width(size.x), height(size.y) {

    }

    constexpr T getRight() const noexcept {
        return left + width;
    }

    constexpr T getBottom() const noexcept {
        return top + height;
    }

    constexpr Vector2<T> getTopLeft() const noexcept {
        return Vector2<T>(left, top);
    }

    constexpr Vector2<T> getCenter() const noexcept {
        return Vector2<T>(left + width / 2, top + height / 2);
    }

    constexpr Vector2<T> getSize() const noexcept {
        return Vector2<T>(width, height);
    }

    constexpr bool contains(const Box<T>& box) const noexcept {
        return left <= box.left && box.getRight() <= getRight() &&
            top <= box.top && box.getBottom() <= getBottom();
    }

    constexpr bool intersects(const Box<T>& box) const noexcept {
        return !(left >= box.getRight() || getRight() <= box.left ||
            top >= box.getBottom() || getBottom() <= box.top);
    }
};

template<typename T, typename GetBox, typename Equal = std::equal_to<T>, typename Float = float>
class Quadtree {
public:
    Quadtree() {}
    Quadtree(const Box<Float>& box, const GetBox& getBox = GetBox(),
        const Equal& equal = Equal()) :
        mBox(box), mRoot(std::make_unique<Node>()), mGetBox(getBox), mEqual(equal) {

    }

    void add(const T& value) {
        add(mRoot.get(), 0, mBox, value);
    }

    void remove(const T& value) {
        remove(mRoot.get(), mBox, value);
    }

    std::vector<T> query(const Box<Float>& box) const {
        auto values = std::vector<T>();
        query(mRoot.get(), mBox, box, values);
        return values;
    }

    std::vector<std::pair<T, T>> findAllIntersections() const {
        auto intersections = std::vector<std::pair<T, T>>();
        findAllIntersections(mRoot.get(), intersections);
        return intersections;
    }

    Box<Float> getBox() const {
        return mBox;
    }

private:
    static constexpr auto Threshold = std::size_t(16);
    static constexpr auto MaxDepth = std::size_t(8);

    struct Node {
        std::array<std::unique_ptr<Node>, 4> children;
        std::vector<T> values;
    };

    Box<Float> mBox;
    std::unique_ptr<Node> mRoot;
    GetBox mGetBox;
    Equal mEqual;

    bool isLeaf(const Node* node) const {
        return !static_cast<bool>(node->children[0]);
    }

    Box<Float> computeBox(const Box<Float>& box, int i) const {
        auto origin = box.getTopLeft();
        auto childSize = box.getSize() / static_cast<Float>(2);
        switch (i) {
            // North West
            case 0:
                return Box<Float>(origin, childSize);
            // Norst East
            case 1:
                return Box<Float>(Vector2<Float>(origin.x + childSize.x, origin.y), childSize);
            // South West
            case 2:
                return Box<Float>(Vector2<Float>(origin.x, origin.y + childSize.y), childSize);
            // South East
            case 3:
                return Box<Float>(origin + childSize, childSize);
            default:
                return Box<Float>();
        }
    }

    int getQuadrant(const Box<Float>& nodeBox, const Box<Float>& valueBox) const {
        auto center = nodeBox.getCenter();
        // West
        if (valueBox.getRight() < center.x) {
            // North West
            if (valueBox.getBottom() < center.y)
                return 0;
            // South West
            else if (valueBox.top >= center.y)
                return 2;
            // Not contained in any quadrant
            else
                return -1;
        }
        // East
        else if (valueBox.left >= center.x) {
            // North East
            if (valueBox.getBottom() < center.y)
                return 1;
            // South East
            else if (valueBox.top >= center.y)
                return 3;
            // Not contained in any quadrant
            else
                return -1;
        }
        // Not contained in any quadrant
        else
            return -1;
    }

    void add(Node* node, std::size_t depth, const Box<Float>& box, const T& value) {
        assert(node != nullptr);
        if(not box.contains(mGetBox(value))) {
            std::cout << box.left << " " << box.getRight() << " " << box.top << " " << box.getBottom() << std::endl;
            std::cout << mGetBox(value).left << " " << mGetBox(value).getRight() << " " << mGetBox(value).top << " " << mGetBox(value).getBottom() << std::endl;
        }
        assert(box.contains(mGetBox(value)));
        if (isLeaf(node)) {
            // Insert the value in this node if possible
            if (depth >= MaxDepth || node->values.size() < Threshold)
                node->values.push_back(value);
            // Otherwise, we split and we try again
            else {
                split(node, box);
                add(node, depth, box, value);
            }
        }
        else {
            auto i = getQuadrant(box, mGetBox(value));
            // Add the value in a child if the value is entirely contained in it
            if (i != -1)
                add(node->children[static_cast<std::size_t>(i)].get(), depth + 1, computeBox(box, i), value);
            // Otherwise, we add the value in the current node
            else
                node->values.push_back(value);
        }
    }

    void split(Node* node, const Box<Float>& box) {
        assert(node != nullptr);
        assert(isLeaf(node) && "Only leaves can be split");
        // Create children
        for (auto& child : node->children)
            child = std::make_unique<Node>();
        // Assign values to children
        auto newValues = std::vector<T>(); // New values for this node
        for (const auto& value : node->values) {
            auto i = getQuadrant(box, mGetBox(value));
            if (i != -1)
                node->children[static_cast<std::size_t>(i)]->values.push_back(value);
            else
                newValues.push_back(value);
        }
        node->values = std::move(newValues);
    }

    bool remove(Node* node, const Box<Float>& box, const T& value) {
        assert(node != nullptr);
        assert(box.contains(mGetBox(value)));
        if (isLeaf(node)) {
            // Remove the value from node
            removeValue(node, value);
            return true;
        }
        else {
            // Remove the value in a child if the value is entirely contained in it
            auto i = getQuadrant(box, mGetBox(value));
            if (i != -1) {
                if (remove(node->children[static_cast<std::size_t>(i)].get(), computeBox(box, i), value))
                    return tryMerge(node);
            }
            // Otherwise, we remove the value from the current node
            else
                removeValue(node, value);
            return false;
        }
    }

    void removeValue(Node* node, const T& value) {
        // Find the value in node->values
        auto it = std::find_if(std::begin(node->values), std::end(node->values),
            [this, &value](const auto& rhs){ return mEqual(value, rhs); });
        assert(it != std::end(node->values) && "Trying to remove a value that is not present in the node");
        // Swap with the last element and pop back
        *it = std::move(node->values.back());
        node->values.pop_back();
    }

    bool tryMerge(Node* node) {
        assert(node != nullptr);
        assert(!isLeaf(node) && "Only interior nodes can be merged");
        auto nbValues = node->values.size();
        for (const auto& child : node->children) {
            if (!isLeaf(child.get()))
                return false;
            nbValues += child->values.size();
        }
        if (nbValues <= Threshold) {
            node->values.reserve(nbValues);
            // Merge the values of all the children
            for (const auto& child : node->children) {
                for (const auto& value : child->values)
                    node->values.push_back(value);
            }
            // Remove the children
            for (auto& child : node->children)
                child.reset();
            return true;
        } else
            return false;
    }

    void query(Node* node, const Box<Float>& box, const Box<Float>& queryBox, std::vector<T>& values) const {
        assert(node != nullptr);
        assert(queryBox.intersects(box));
        for (const auto& value : node->values) {
            if (queryBox.intersects(mGetBox(value)))
                values.push_back(value);
        }
        if (!isLeaf(node)) {
            for (auto i = std::size_t(0); i < node->children.size(); ++i) {
                auto childBox = computeBox(box, static_cast<int>(i));
                if (queryBox.intersects(childBox))
                    query(node->children[i].get(), childBox, queryBox, values);
            }
        }
    }

    void findAllIntersections(Node* node, std::vector<std::pair<T, T>>& intersections) const {
        // Find intersections between values stored in this node
        // Make sure to not report the same intersection twice
        for (auto i = std::size_t(0); i < node->values.size(); ++i) {
            for (auto j = std::size_t(0); j < i; ++j) {
                if (mGetBox(node->values[i]).intersects(mGetBox(node->values[j])))
                    intersections.emplace_back(node->values[i], node->values[j]);
            }
        }
        if (!isLeaf(node)) {
            // Values in this node can intersect values in descendants
            for (const auto& child : node->children) {
                for (const auto& value : node->values)
                    findIntersectionsInDescendants(child.get(), value, intersections);
            }
            // Find intersections in children
            for (const auto& child : node->children)
                findAllIntersections(child.get(), intersections);
        }
    }

    void findIntersectionsInDescendants(Node* node, const T& value, std::vector<std::pair<T, T>>& intersections) const {
        // Test against the values stored in this node
        for (const auto& other : node->values) {
            if (mGetBox(value).intersects(mGetBox(other)))
                intersections.emplace_back(value, other);
        }
        // Test against values stored into descendants of this node
        if (!isLeaf(node)) {
            for (const auto& child : node->children)
                findIntersectionsInDescendants(child.get(), value, intersections);
        }
    }
};

}

template<typename T>
struct GetBox {
    using Tri = typename DelaunayIndexer2D<T>::Tri;
    inline quadtree::Box<T> operator()(Tri const & tri) const {
        T x_min = std::min(tri.p0.x, std::min(tri.p1.x, tri.p2.x));
        T x_max = std::max(tri.p0.x, std::max(tri.p1.x, tri.p2.x));
        T y_min = std::min(tri.p0.y, std::min(tri.p1.y, tri.p2.y));
        T y_max = std::max(tri.p0.y, std::max(tri.p1.y, tri.p2.y));
        return quadtree::Box<T>(x_min, y_min, x_max - x_min, y_max - y_min);
    }
    inline quadtree::Box<T> operator()(Tri const * tri) const {
        if(tri == nullptr)
            return quadtree::Box<T>(std::numeric_limits<T>::infinity(), std::numeric_limits<T>::infinity(), 0, 0);
        else
            return this->operator()(*tri);
    }
};

template<typename T>
class DelaunayIndexer2D {
    using Simplex = typename IDelaBella2<T>::Simplex;
    using Vertex = typename IDelaBella2<T>::Vertex;
public:
    struct Point {
        T x;
        T y;
        int i=-1;;
        bool operator==(Point const & o) const {
            return x == o.x and y == o.y and i == o.i;
        }
    };

    struct Edge {
        int a;
        int b;
    };

    struct Tri {
        Point p0;
        Point p1;
        Point p2;
        bool operator==(Tri const & o) const {
            return p0 == o.p0 and p1 == o.p1 and p2 == o.p2;
        }
    };
private:
    quadtree::Quadtree<Tri const *, GetBox<T>, std::equal_to<Tri const *>, T> q_tree;
    std::vector<Tri> all_tris;

    static inline bool PointInTriangle (Point const & p, Tri const & tri) {
        Point const & v0 = tri.p0;
        Point const & v1 = tri.p1;
        Point const & v2 = tri.p2;
        T s = (p.y * (v0.x - v2.x) + v0.y * v2.x - v0.x * v2.y + p.x * (-v0.y + v2.y))
            / (-v1.y * v2.x + v0.y * (-v1.x + v2.x) + v0.x * (v1.y - v2.y) + v1.x * v2.y);
        T t = (p.y * (v0.x - v1.x) + v0.y * v1.x - v0.x * v1.y + p.x * (-v0.y + v1.y))
            / (v0.y * (v1.x - v2.x) + v1.y * v2.x - v1.x * v2.y + v0.x * (-v1.y + v2.y));
        return s > 0 and t > 0 and (1 - s - t) > 0;
    }

    static inline bool PointInTriangle (Point const & p, Tri const * tri) {
        if(tri == nullptr)
            return false;
        else
            return PointInTriangle(p, *tri);
    }

public:
    DelaunayIndexer2D() {}
    DelaunayIndexer2D(DelaunayIndexer2D const &) = default;
    DelaunayIndexer2D(DelaunayIndexer2D &&) = default;
    DelaunayIndexer2D(DelaunayIndexer2D &) = default;
    DelaunayIndexer2D & operator=(DelaunayIndexer2D const &) = default;
    DelaunayIndexer2D & operator=(DelaunayIndexer2D &&) = default;
    DelaunayIndexer2D & operator=(DelaunayIndexer2D &) = default;
    DelaunayIndexer2D(std::vector<T> const & x, std::vector<T> const & y) {
        all_tris.reserve(x.size() / 3 + 1);
        T min_x = *std::min_element(x.begin(), x.end());
        T max_x = *std::max_element(x.begin(), x.end());
        T min_y = *std::min_element(y.begin(), y.end());
        T max_y = *std::max_element(y.begin(), y.end());
        T width_x = max_x - min_x;
        T width_y = max_y - min_y;
        Tri tri;
        quadtree::Box<T> bounding_box(min_x - width_x * 1e-4, min_y - width_y * 1e-4, width_x * (1 + 2e-4), width_y * (1 + 2e-4));
        q_tree = quadtree::Quadtree<Tri const *, GetBox<T>, std::equal_to<Tri const *>, T>(bounding_box);
        int verts;
        size_t npoly;
        Simplex const * dela;

        IDelaBella2<T> * idb(IDelaBella2<T>::Create());
        Point * cloud = new Point[x.size()];
        std::map<std::tuple<double, double>, size_t> index_map;
        for(size_t i=0; i<x.size(); ++i) {
            cloud[i].x = x[i];
            cloud[i].y = y[i];
            cloud[i].i = i;
            index_map.insert({{x[i], y[i]}, i});
        }

        verts = idb->Triangulate(x.size(), &(cloud->x), &(cloud->y), sizeof(Point));
        if(verts <= 0) {
            throw std::runtime_error("Could not triangulate input grid");
        }

        npoly = idb->GetNumPolygons();
        dela = idb->GetFirstDelaunaySimplex();

        for(size_t i=0; i<npoly; ++i) {
            tri.p0.x = dela->v[0]->x; tri.p0.y = dela->v[0]->y; tri.p0.i = dela->v[0]->i;
            tri.p1.x = dela->v[1]->x; tri.p1.y = dela->v[1]->y; tri.p1.i = dela->v[1]->i;
            tri.p2.x = dela->v[2]->x; tri.p2.y = dela->v[2]->y; tri.p2.i = dela->v[2]->i;
            all_tris.emplace_back(tri);
            dela = dela->next;
        }

        for(size_t i=0; i<npoly; ++i) {
            q_tree.add(&all_tris[i]);
        }
        // Frees memory allocated in idb
        idb->Destroy(); // Also frees memory allocated above in: Point * cloud = new Point[x.size()];
    }

    virtual Tri const * operator()(T const & x, T const & y) const {
        Point p; p.x = x; p.y = y;
        quadtree::Box<T> box(x, y, 0, 0);
        std::vector<Tri const *> tris = q_tree.query(box);
        for(size_t i=0; i<tris.size(); ++i) {
            if(PointInTriangle(p, tris[i])) {
                return tris[i];
            }
        }
        for(size_t i=0; i<all_tris.size(); ++i) {
            if(PointInTriangle(p, all_tris[i])) {
                return &all_tris[i];
            }
        }
        return nullptr;
    }
};

template<typename T>
class Interpolator2D {
public:
    virtual T operator()(T x, T y) const = 0;
};

template<typename T>
class GridLinearInterpolator2D : public Interpolator2D<T> {
    // Data
    std::vector<T> t_x;
    std::vector<T> t_y;
    std::vector<T> t_z;

    // Indexer
    IrregularGridIndexer2D<T> indexer;

    // Interpolation info
    std::shared_ptr<BiLinearInterpolationOperator<T>> interp;
    std::shared_ptr<Transform<T>> x_transform;
    std::shared_ptr<Transform<T>> y_transform;
    std::shared_ptr<Transform<T>> z_transform;

public:
    GridLinearInterpolator2D(std::vector<T> x, std::vector<T> y, std::vector<T> z, std::shared_ptr<BiLinearInterpolationOperator<T>> interp, std::shared_ptr<Transform<T>> x_transform, std::shared_ptr<Transform<T>> y_transform, std::shared_ptr<Transform<T>> z_transform)
        : t_x(), t_y(), t_z(), indexer(x, y), interp(interp), x_transform(x_transform), y_transform(y_transform), z_transform(z_transform) {
        t_x.reserve(x.size());
        t_y.reserve(y.size());
        t_z.reserve(y.size());
        for(size_t i=0; i<x.size(); ++i) {
            t_x.push_back(x_transform->operator()(x[i]));
            t_y.push_back(y_transform->operator()(y[i]));
            t_z.push_back(z_transform->operator()(z[i]));
        }
    }

    virtual T operator()(T x, T y) const override {
        Index2D indices = indexer(x);
        T const & t_x0 = t_x[indices.x0];
        T const & t_y0 = t_y[indices.y0];
        T const & t_x1 = t_x[indices.x1];
        T const & t_y1 = t_y[indices.x1];
        T const & t_z00 = t_z[indices.z00];
        T const & t_z01 = t_z[indices.z01];
        T const & t_z10 = t_z[indices.z10];
        T const & t_z11 = t_z[indices.z11];

        T t_x = x_transform->Function(x);
        T t_y = y_transform->Function(y);

        T t_res = interp->operator()(t_x0, t_y0, t_x1, t_y1, t_z00, t_z01, t_z10, t_z11, t_x, t_y);
        T res = z_transform->Inverse(t_res);
        return res;
    }
};

template<typename T>
class LinearDelaunayInterpolator2D : public Interpolator2D<T> {
    using Simplex = typename IDelaBella2<T>::Simplex;
    using Vertex = typename IDelaBella2<T>::Vertex;
    using Tri = typename DelaunayIndexer2D<T>::Tri;
    // Data
    std::vector<T> t_x;
    std::vector<T> t_y;
    std::vector<T> t_z;

    // Indexer
    DelaunayIndexer2D<T> indexer;

    // Interpolation info
    std::shared_ptr<SimplexLinearInterpolationOperator<T>> interp;
    std::shared_ptr<Transform<T>> x_transform;
    std::shared_ptr<Transform<T>> y_transform;
    std::shared_ptr<Transform<T>> z_transform;
public:
    LinearDelaunayInterpolator2D(std::vector<T> const & x, std::vector<T> const & y, std::vector<T> const & z, std::shared_ptr<Transform<T>> x_transform, std::shared_ptr<Transform<T>> y_transform, std::shared_ptr<Transform<T>> z_transform, std::shared_ptr<SimplexLinearInterpolationOperator<T>> interp)
        : t_x(), t_y(), t_z(), indexer(x, y), interp(interp), x_transform(x_transform), y_transform(y_transform), z_transform(z_transform) {
        t_x.reserve(x.size());
        t_y.reserve(y.size());
        t_z.reserve(y.size());
        for(size_t i=0; i<x.size(); ++i) {
            t_x.push_back(x_transform->Function(x[i]));
            t_y.push_back(y_transform->Function(y[i]));
            t_z.push_back(z_transform->Function(z[i]));
        }
        indexer = DelaunayIndexer2D<T>(t_x, t_y);
    }

    virtual T operator()(T x, T y) const override {
        Tri const *tri = indexer(x, y);
        if(tri == nullptr)
            throw std::runtime_error("Simplex ptr is NULL");
        T t_z0 = t_z[tri->p0.i];
        T t_z1 = t_z[tri->p1.i];
        T t_z2 = t_z[tri->p2.i];

        T t_res = interp->operator()(x, y, tri, t_z0, t_z1, t_z2);
        return z_transform->Inverse(t_res);
    }
};

} // math
} // LI

CEREAL_CLASS_VERSION(siren::math::Transform<double>, 0);

CEREAL_CLASS_VERSION(siren::math::IdentityTransform<double>, 0);
CEREAL_REGISTER_TYPE(siren::math::IdentityTransform<double>);
CEREAL_REGISTER_POLYMORPHIC_RELATION(siren::math::Transform<double>, siren::math::IdentityTransform<double>);

CEREAL_CLASS_VERSION(siren::math::LogTransform<double>, 0);
CEREAL_REGISTER_TYPE(siren::math::LogTransform<double>);
CEREAL_REGISTER_POLYMORPHIC_RELATION(siren::math::Transform<double>, siren::math::LogTransform<double>);

CEREAL_CLASS_VERSION(siren::math::SymLogTransform<double>, 0);
CEREAL_REGISTER_TYPE(siren::math::SymLogTransform<double>);
CEREAL_REGISTER_POLYMORPHIC_RELATION(siren::math::Transform<double>, siren::math::SymLogTransform<double>);

CEREAL_CLASS_VERSION(siren::math::RangeTransform<double>, 0);
CEREAL_REGISTER_TYPE(siren::math::RangeTransform<double>);
CEREAL_REGISTER_POLYMORPHIC_RELATION(siren::math::Transform<double>, siren::math::RangeTransform<double>);

CEREAL_CLASS_VERSION(siren::math::LinearInterpolationOperator<double>, 0);

CEREAL_CLASS_VERSION(siren::math::DropLinearInterpolationOperator<double>, 0);
CEREAL_REGISTER_TYPE(siren::math::DropLinearInterpolationOperator<double>);
CEREAL_REGISTER_POLYMORPHIC_RELATION(siren::math::LinearInterpolationOperator<double>, siren::math::DropLinearInterpolationOperator<double>);

CEREAL_CLASS_VERSION(siren::math::Indexer1D<double>, 0);

CEREAL_CLASS_VERSION(siren::math::RegularIndexer1D<double>, 0);
CEREAL_REGISTER_TYPE(siren::math::RegularIndexer1D<double>);
CEREAL_REGISTER_POLYMORPHIC_RELATION(siren::math::Indexer1D<double>, siren::math::RegularIndexer1D<double>);

CEREAL_CLASS_VERSION(siren::math::IrregularIndexer1D<double>, 0);
CEREAL_REGISTER_TYPE(siren::math::IrregularIndexer1D<double>);
CEREAL_REGISTER_POLYMORPHIC_RELATION(siren::math::Indexer1D<double>, siren::math::IrregularIndexer1D<double>);

CEREAL_CLASS_VERSION(siren::math::TransformIndexer1D<double>, 0);
CEREAL_REGISTER_TYPE(siren::math::TransformIndexer1D<double>);
CEREAL_REGISTER_POLYMORPHIC_RELATION(siren::math::Indexer1D<double>, siren::math::TransformIndexer1D<double>);

#endif // SIREN_Interpolation_H
