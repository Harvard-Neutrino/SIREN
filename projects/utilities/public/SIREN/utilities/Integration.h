#pragma once
#ifndef SIREN_Integration_H
#define SIREN_Integration_H

#include <cmath>
#include <vector>
#include <cassert>

namespace siren {
namespace utilities {

/**
* An object which performs 1-D intgration of a function
* using the trapezoid rule, and allows the approximation
* to be refined incrementally.
*/
template<typename FuncType>
struct TrapezoidalIntegrator{
private:
  /**
   * The function being integrated
   */
  const FuncType& f;
  /**
   * The lower bound of integration
   */
  double a;
  /**
   * The upper bound of integration
   */
  double b;
  /**
   * Counter expressing the number of times the integral
   * approximation has been refined
   */
  unsigned int currentDetail;
  /**
   * The current approximation of the integral
   */
  double value;

  /**
   * Add one level of detail to the integral approximation
   */
  void update(){
     currentDetail++;
     if(currentDetail==1)
        value=(b-a)*(f(a)+f(b))/2;
     else{
        unsigned long npoints=1ul<<(currentDetail-2);
        double dx=(b-a)/npoints;
        double x=a+dx/2;
        double sum=0.0;
        for(unsigned long i=0; i<npoints; i++, x+=dx)
           sum+=f(x);
        value=(value+(b-a)*sum/npoints)/2;
     }
  }

public:
  /**
   * @param f the function to be integrated
   * @param a the lower bound of integration
   * @param b the upper bound of integration
   */
  TrapezoidalIntegrator(const FuncType& f, double a, double b):
  f(f),a(a),b(b),currentDetail(0),value(0){
     if(a>b)
        std::swap(a,b);
  }

  /**
   * Get the integral approximation, updating it with higher
   * detail if necessary
   *
   * @param detail how finely to approximate the integral.
   *               A detail level of n requires 1+2^n function evaluations,
   *               but reuses any evaluations already performed when lower
   *               detail levels were calculated.
   */
  double integrate(unsigned int detail){
     detail++;
     while(currentDetail<detail)
        update();
     return(value);
  }

  /**
   * Get the current detail level of the ingral approximation.
   */
  unsigned int getDetail() const{ return(currentDetail); }
};


/**
* @brief Performs a fifth order Romberg integration of a function to a chosen tolerance.
*
* This routine is rather simplistic and not suitable for complicated functions,
* particularly not ones with discontinuities, but it is very fast for smooth functions.
*
* @param func the function to integrate
* @param a the lower bound of integration
* @param b the upper bound of integration
* @param tol the (absolute) tolerance on the error of the integral
*/
template<typename FuncType>
double rombergIntegrate(const FuncType& func, double a, double b, double tol=1e-6){
  const unsigned int order=5;
  const unsigned int maxIter=20;
  if(tol<0)
     throw(std::runtime_error("Integration tolerance must be positive"));

  std::vector<double> stepSizes, estimates, c(order), d(order);
  stepSizes.push_back(1);

  TrapezoidalIntegrator<FuncType> t(func,a,b);
  for(unsigned int i=0; i<maxIter; i++){
     //refine the integral estimate
     estimates.push_back(t.integrate(t.getDetail()));

     if(i>=(order-1)){ //if enough estimates have been accumulated
        //extrapolate to zero step size
        const unsigned int baseIdx=i-(order-1);
        std::copy(estimates.begin()+baseIdx,estimates.begin()+baseIdx+order, c.begin());
        std::copy(estimates.begin()+baseIdx,estimates.begin()+baseIdx+order, d.begin());

        double ext=estimates.back(), extErr;
        for(unsigned int m=1; m<order; m++){
           for(unsigned int j=0; j<order-m; j++){
              double ho=stepSizes[j+baseIdx];
              double hp=stepSizes[j+m+baseIdx];
              double w=c[j+1]-d[j];
              double den=ho-hp;
              assert(den!=0.0);
              den=w/den;
              c[j]=ho*den;
              d[j]=hp*den;
           }
           extErr=d[order-1-m];
           ext+=extErr;
        }

        //declare victory if the tolerance criterion is met
        if(std::abs(extErr)<=tol*std::abs(ext))
           return(ext);
     }
     //prepare for next step
     stepSizes.push_back(stepSizes.back()/4);
  }
  throw(std::runtime_error("Integral failed to converge"));
}

}
}

#endif // SIREN_Integration_H
