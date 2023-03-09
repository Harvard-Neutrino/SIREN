#ifndef LI_Interpolation_H
#define LI_Interpolation_H

namespace LI {
namespace math {

template<typename T>
class Transform {
public:
    virtual T Function(T x) const = 0;
    virtual T Inverse(T x) const = 0;
};

template<typename T>
class IdentityTransform : public Transform {
public:
    IdentityTransform() {}
    virtual T Function(T x) const override {
        return x;
    }
    virtual T Inverse(T x) const override {
        return x;
    }
};

template<typename T>
class GenericTransform : public Transform {
    std::function<T(T)> function;
    std::function<T(T)> inverse;
public:
    GenericTransform(std::function<T(T)> function, std::function<T(T)> inverse)
        : function(function), inverse(inverse) {}
    virtual T Function(T x) const override {
        return function(x);
    }
    virtual T Inverse(T x) const override {
        return inverse(x);
    }
};

template<typename T>
class LogTransform : public Transform {
public:
    LogTransform() {}
    virtual T Function(T x) const override {
        return log(x);
    }
    virtual T Inverse(T x) const override {
        return exp(x);
    }
};

template<typename T>
class SymLogTransform : public Transform {
    T min_x;
    T log_min_x;
public:
    SymLogTransform(T min_x) : min_x(std::abs(min_x)), log_min_x(log(std::abs(min_x))) {}
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
            return std::copysign(exp(x - min_x + log_min_x), x);
        }
    }
};

template<typename T>
class RangeTransform : public Transform {
    T min_x;
    T range;
public:
    RangeTransform(T min_x, T max_x) : min_x(min_x), range(max_x - min_x) {}
    virtual T Function(T x) const override {
        return (x - min_x) / range;
    }
    virtual T Inverse(T x) const override {
        return x * range + min_x;
    }
};

template<typename T>
class FunctionalRangeTransform : public Transform {
    std::function<T(T)> min_function;
    std::function<T(T)> max_function;
public:
    RangeTransform(std::function<T(T)> min_function, std::function<T(T)> max_function)
        : min_function(min_function), max_function(max_function) {}
    virtual T Function(T x) const override {
        T min_x = min_function(x);
        T range = max_function(x) - min_x;
        return (x - min_x) / range;
    }
    virtual T Inverse(T x) const override {
        T min_x = min_function(x);
        T range = max_function(x) - min_x;
        return x * range + min_x;
    }
};

class LinearInterpolator {
    virtual T operator(T const & x0, T const & x1, T const & y0, T const & y1, T const & x) const {
        T delta_x = x1 - x0;
        T delta_y = y1 - y0;
        return (x - x0) * delta_y / delta_x;
    }
};

class DropLinearInterpolator : public LinearInterpolator {
    virtual T operator(T const & x0, T const & x1, T const & y0, T const & y1, T const & x) const override {
        if(y0 == 0 or y1 == 0)
            return 0;
        T delta_x = x1 - x0;
        T delta_y = y1 - y0;
        return (x - x0) * delta_y / delta_x;
    }
};

template<typename T>
std::tuple<std::shared_ptr<Transform<T>>, std::shared_ptr<Transform<T>>> DetermineInterpolationSpace(
        std::vector<T> const & x,
        std::vector<T> const & y,
        std::shared_ptr<LinearInterpolator> interp) {
    std::vector<T> symlog_x(x.size());
    std::vector<T> symlog_y(y.size());
    T min_x = 1;
    T min_y = 1;
    bool have_x = false;
    bool have_y = false;
    for(size_t i=0; i<x.size(); ++i) {
        if(x[i] != 0)
            if(have_x) {
                min_x = std::min(min_x, std::abs(x[i]));
                have_x = true;
            } else
                min_x = std::abs(x[i]);
        if(y[i] != 0)
            if(have_y) {
                min_y = std::min(min_y, std::abs(y[i]));
                have_y = true;
            } else
                min_y = std::abs(y[i]);
    }
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
        T y_i_estimate = interp->operator()(x[i-1], x[i+1], y[i-1], y[i+1], x[i]);
        T delta = y_i_estimate - y[i];
        tot_00 += delta * delta;

        y_i_estimate = interp->operator()(x[i-1], x[i+1], symlog_y[i-1], symlog_y[i+1], x[i]);
        T delta = symlog_t_y.Inverse(y_i_estimate) - y[i];
        tot_01 += delta * delta;

        T y_i_estimate = interp->operator()(symlog_x[i-1], symlog_x[i+1], y[i-1], y[i+1], symlog_x[i]);
        T delta = y_i_estimate - y[i];
        tot_10 += delta * delta;

        y_i_estimate = interp->operator()(symlog_x[i-1], symlog_x[i+1], symlog_y[i-1], symlog_y[i+1], symlog_x[i]);
        T delta = symlog_t_y.Inverse(y_i_estimate) - y[i];
        tot_11 += delta * delta;
    }

    std::vector<T> tots = {tot_00, tot_01, tot_10, tot_11};
    int min_combination = std::distance(tots.begin(), std::max_element(tots.begin(), tots.end()));

    switch(min_combination) {
        case 1:
            return {std::make_shared<IdentityTransform<T>>(), std::make_shared<SymLogTransform<T>>(min_y)};
        case 2:
            return {std::make_shared<SymLogTransform<T>>(min_x), std::make_shared<IdentityTransform<T>>()};
        case 3:
            return {std::make_shared<SymLogTransform<T>>(min_x), std::make_shared<SymLogTransform<T>>(min_y)};
        default:
            return {std::make_shared<IdentityTransform<T>>(), std::make_shared<IdentityTransform<T>>()};
    }
};

}
}

#endif // LI_Interpolation_H
