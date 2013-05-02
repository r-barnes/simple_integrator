#ifndef __integrator_
#define __integrator_

#include <functional>
#include <cmath>
#include <iostream>
#include <array>

template <class T>
class integrator {
  private:
    T stateval;
    typedef std::function<void(const T &state, T &dxdt, double t)> dx_type;
    double dtval, dtmax, dtmin;
    dx_type dx;
    double t;
    int goodsteps;
    int stepcount;
  public:
    integrator(const T &stateval, dx_type dx, double dtmax, double dtmin);
    double time() const;
    T& state();
    double dt() const;
    void dt(double h);
    void step();
    int steps() const;
};

template <class T>
integrator<T>::integrator(const T &stateval, dx_type dx, double dtmax, double dtmin) : stateval(stateval), dx(dx), dtmax(dtmax), dtmin(dtmin) {
  dtval=dtmin;
  stepcount=0;
  t=0;
}

template<class T>
double integrator<T>::time() const { return t; }

template<class T>
T& integrator<T>::state() { return stateval; }

template<class T>
double integrator<T>::dt() const { return dtval; }

template<class T>
void integrator<T>::dt(double h) { dtval=h; }

template<class T>
int integrator<T>::steps() const { return stepcount; }

template<class T>
void integrator<T>::step() {
  T e1, e2;

  ++stepcount;

  dx(stateval           , e1, t);
  dx(stateval+e1*dtval/2, e2, t+dtval/2.);
  double abs_e1=abs(stateval+e1*dtval);
  double abs_e2=abs(stateval+e1*dtval/2+e2*dtval/2);

  if( (std::abs(abs_e1-abs_e2)/(std::abs(abs_e1+abs_e2)/2))>0.05 ){
    stateval+=e1*dtval/2.;
    t+=dtval/2.;
    dtval/=2.;
    if(dtval<dtmin) dtval=dtmin;
    goodsteps=0;
  } else {
    stateval+=e1*dtval;
    t+=dtval;
    ++goodsteps;
  }

  if(dtval<dtmax && goodsteps>=8)
    dtval*=2;
  if(dtval>dtmax)
    dtval=dtmax;
}

#endif
