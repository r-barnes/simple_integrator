#ifndef __integrator_
#define __integrator_

#include <functional>
#include <cmath>
#include <array>
#include <cassert>
#include "calendar_queue.hpp"

template <class T>
class integrator {
  protected:
    T stateval;
    typedef std::function<void(const T &state, T &dxdt, double t)> dx_type;
    dx_type dx;
    double dtval, dtmax, dtmin;
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
  assert(dtmin>0);
  assert(dtmax>0);
  assert(dtmin<=dtmax);

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
void integrator<T>::dt(double h) {
  assert(h>0);
  assert(h>=dtmin);
  assert(h<=dtmax);
  dtval=h;
}

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






template <class T>
class event_integrator : public integrator<T> {
  private:
    CalendarQueue calq;
    bool at_event;
    std::string at_event_name;
  public:
    event_integrator(const T &stateval, typename integrator<T>::dx_type dx, double dtmax, double dtmin) : integrator<T>(stateval, dx, dtmax, dtmin) {
      at_event=false;
      at_event_name="";
    }
    void insert_event(double t, const std::string &event, double recur_int);
    void insert_event(double t, const char event[], double recur_int);
    bool is_event() const;
    std::string event() const;
    void step();
};

template<class T>
void event_integrator<T>::insert_event(double t, const std::string &event, double recur_int=0){
  calq.insert(t,event,recur_int);
}

template<class T>
void event_integrator<T>::insert_event(double t, const char event[], double recur_int=0){
  calq.insert(t,event,recur_int);
}

template<class T>
bool event_integrator<T>::is_event() const {
  return at_event;
}

template<class T>
std::string event_integrator<T>::event() const {
  if( is_event() )
    return calq.current_event();
  else
    return "";
}

template<class T>
void event_integrator<T>::step() {
  T e1, e2;

  if(at_event && calq.current_time()==integrator<T>::t){
    at_event_name=calq.current_event();
    calq.pop();
    return;
  }

  at_event=false;

  ++integrator<T>::stepcount;

  if(!calq.empty()){
    while(
        integrator<T>::dtval > integrator<T>::dtmin
        && integrator<T>::t + integrator<T>::dtval > calq.current_time()
    ) {
      integrator<T>::dtval/=2.;
    }

    if(integrator<T>::dtval < integrator<T>::dtmin)
      integrator<T>::dtval = integrator<T>::dtmin;

    if(integrator<T>::dtval == integrator<T>::dtmin && integrator<T>::t + integrator<T>::dtval > calq.current_time()){
      integrator<T>::dtval=calq.current_time()-integrator<T>::t;
      integrator<T>::dx(integrator<T>::stateval, e1, integrator<T>::t);
      integrator<T>::stateval+=e1*integrator<T>::dtval;
      integrator<T>::t=calq.current_time();
      integrator<T>::dtval = integrator<T>::dtmin;
      at_event=true;
      at_event_name=calq.current_event();
      calq.pop();
      return;
    }
  }    

  integrator<T>::dx(integrator<T>::stateval           , e1, integrator<T>::t);
  integrator<T>::dx(integrator<T>::stateval + e1*integrator<T>::dtval/2, e2, integrator<T>::t+integrator<T>::dtval/2.);
  double abs_e1=abs(integrator<T>::stateval + e1*integrator<T>::dtval);
  double abs_e2=abs(integrator<T>::stateval + e1*integrator<T>::dtval/2 + e2*integrator<T>::dtval/2);

  if( (std::abs(abs_e1-abs_e2)/(std::abs(abs_e1+abs_e2)/2))>0.05 ){
    integrator<T>::stateval+=e1*integrator<T>::dtval/2.;
    integrator<T>::t+=integrator<T>::dtval/2.;
    integrator<T>::dtval/=2.;
    if(integrator<T>::dtval<integrator<T>::dtmin)
      integrator<T>::dtval=integrator<T>::dtmin;
    integrator<T>::goodsteps=0;
  } else {
    integrator<T>::stateval+=e1*integrator<T>::dtval;
    integrator<T>::t+=integrator<T>::dtval;
    ++integrator<T>::goodsteps;
  }

  if(integrator<T>::dtval<integrator<T>::dtmax && integrator<T>::goodsteps>=8)
    integrator<T>::dtval*=2;
  if(integrator<T>::dtval>integrator<T>::dtmax)
    integrator<T>::dtval=integrator<T>::dtmax;
}

#endif
