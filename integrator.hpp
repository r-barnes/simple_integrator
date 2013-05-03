#ifndef __integrator_
#define __integrator_

#include <functional>
#include <cmath>
#include <array>
#include <cassert>
#include <queue>

namespace si_lib {

/**
  The calendar queue class stores a (time,event_name) pair.
  It is essentially a priority queue whose top priority is always increasing.
*/
class CalendarQueue {
	private:
    /**
      This class stores an event. An event is represented by a floating-point
      time and a string identifying what the event is.
    */
		class DiscreteEvent {
			public:
				double t;
				std::string event;
        double recur_int;
				DiscreteEvent (double t, const char event[], double recur_int)       : t(t), event(event), recur_int(recur_int) {}
				DiscreteEvent (double t, const std::string &event, double recur_int) : t(t), event(event), recur_int(recur_int) {}
        bool operator< (const DiscreteEvent& a) const {
          return t>a.t;
        }
		};

		typedef std::priority_queue<DiscreteEvent, std::vector<DiscreteEvent> > deq_type;
		deq_type deq;
	public:
    ///Schedule \a event at \a t
		void insert(double t, const char event[], double recur_int=0){
      assert(recur_int>=0);
      //assert(empty() || t>=current_time());
			deq.push(DiscreteEvent(t,event,recur_int));
		}
    ///Schedule \a event at \a t
		void insert(double t, const std::string &event, double recur_int=0){
      assert(recur_int>=0);
      //assert(empty() || t>=current_time());
			deq.push(DiscreteEvent(t,event,recur_int));
		}
    ///Time the top element of the calendar queue is scheduled for
		double current_time() const {
			return deq.top().t;
		}
    ///Event indicated by the top element of the calendar queue
		std::string current_event() const {
			return deq.top().event;
		}
    ///Remove the top element of the calendar queue, and cast it into oblivion
		void pop() {
      if(deq.empty())
        return;

      if(deq.top().recur_int>0)
        insert(current_time()+deq.top().recur_int, deq.top().event, deq.top().recur_int);
			deq.pop();
		}
    ///Take the top element of the queue and insert a copy of it \a dt
    ///into the future
		double reschedule_top(double dt) {
			deq.push(DiscreteEvent(deq.top().t+dt, deq.top().event, deq.top().recur_int));
			return current_time()+dt;
		}
    ///Clear all events from the calendar queue
		void clear() {
			deq=deq_type();
		};
    ///Is the calendar queue empty?
		bool empty() const {
			return deq.empty();
		}
};







template <class T, int N>
class arraystate : public std::array<T,N> {
  public:
    arraystate() {}
    arraystate(const std::vector<T> &init){
      for(unsigned int i=0;i<std::array<T,N>::size();++i)
        std::array<T,N>::operator[](i)=init[i];
    }

    arraystate<T,N>& operator+=(const arraystate<T,N> &other) {
      for(unsigned int i=0;i<std::array<T,N>::size();++i)
        std::array<T,N>::operator[](i)+=other[i];
      return *this;
    }
};

template <class T, int N>
arraystate<T,N> operator+( const arraystate<T,N> &a , double b ){
  arraystate<T,N> result(a);
  for(double& i: result)
    i+=b;
  return result;
}

template <class T, int N>
arraystate<T,N> operator+( double a, const arraystate<T,N> &b ){
  return b+a;
}

template <class T, int N>
arraystate<T,N> operator+( const arraystate<T,N> &a, const arraystate<T,N> &b ){
  arraystate<T,N> result(a);
  for(unsigned int i=0;i<result.size();++i)
    result[i]+=b[i];
  return result;
}

template <class T, int N>
arraystate<T,N> operator*( const arraystate<T,N> &a, double b ){
  arraystate<T,N> result(a);
  for(double& i: result)
    i*=b;
  return result;
}

template <class T, int N>
arraystate<T,N> operator*( double a, const arraystate<T,N> &b ){
  return b*a;
}

template <class T, int N>
arraystate<T,N> operator/( const arraystate<T,N> &a, double b){
  arraystate<T,N> result(a);
  for(double& i: result)
    i/=b;
  return result;
}

template <class T, int N>
double abs( const arraystate<T,N> &a ){
  double result=0;
  for(const double& i: a)
    result+=std::abs(i);
  return result;
}

template <class T, int N>
std::ostream& operator<<(std::ostream &out, const arraystate<T,N> &a){
  for(typename arraystate<T,N>::const_iterator i=a.begin();i!=a.end();++i)
    out<<*i<<" ";
  return out;
}





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

}

#endif
