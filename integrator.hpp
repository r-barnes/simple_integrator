#ifndef __Integrator_
#define __Integrator_

#include <functional>
#include <array>
#include <cassert>
#include <queue>
#include <cmath>

namespace si_lib {

double abs(const double &a){
  if(a<0)
    return -a;
  else
    return a;
}

/**
  The calendar queue class stores a (time,event_name,recurrence_interval) tuple.
  It is essentially a priority queue whose top priority is always increasing.
*/
class CalendarQueue {
	private:
    /**
      This class stores an event. An event is represented by a floating-point
      time and a string identifying what the event is. The event may recur at
      a regular interval.
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




/**
  The ArrayState class wraps std::array providing algebra operations used by
  the Integrator.
*/
template <class T, int N>
class ArrayState : public std::array<T,N> {
  public:
    ArrayState() {}
    ArrayState(const std::vector<T> &init){
      for(unsigned int i=0;i<std::array<T,N>::size();++i)
        std::array<T,N>::operator[](i)=init[i];
    }

    ArrayState<T,N>& operator+=(const ArrayState<T,N> &other) {
      for(unsigned int i=0;i<std::array<T,N>::size();++i)
        std::array<T,N>::operator[](i)+=other[i];
      return *this;
    }
};

template <class T, int N>
ArrayState<T,N> operator+( const ArrayState<T,N> &a , double b ){
  ArrayState<T,N> result(a);
  for(T& i: result)
    i+=b;
  return result;
}

template <class T, int N>
ArrayState<T,N> operator+( double a, const ArrayState<T,N> &b ){
  return b+a;
}

template <class T, int N>
ArrayState<T,N> operator+( const ArrayState<T,N> &a, const ArrayState<T,N> &b ){
  ArrayState<T,N> result(a);
  for(unsigned int i=0;i<result.size();++i)
    result[i]+=b[i];
  return result;
}

template <class T, int N>
ArrayState<T,N> operator*( const ArrayState<T,N> &a, double b ){
  ArrayState<T,N> result(a);
  for(T& i: result)
    i*=b;
  return result;
}

template <class T, int N>
ArrayState<T,N> operator*( double a, const ArrayState<T,N> &b ){
  return b*a;
}

template <class T, int N>
ArrayState<T,N> operator/( const ArrayState<T,N> &a, double b){
  ArrayState<T,N> result(a);
  for(T& i: result)
    i/=b;
  return result;
}

template <class T, int N>
double abs( const ArrayState<T,N> &a ){
  double result=0;
  for(const T& i: a)
    result+=abs(i);
  return result;
}

template <class T, int N>
std::ostream& operator<<(std::ostream &out, const ArrayState<T,N> &a){
  for(typename ArrayState<T,N>::const_iterator i=a.begin();i!=a.end();++i)
    out<<*i<<" ";
  return out;
}






/**
  The Integrator class supplies a simple Euler's method integrator with an
  adaptive step-size scheme
*/
template <class T>
class Integrator {
  protected:
    T stateval;
    typedef std::function<void(const T &state, T &dxdt, double t)> dx_type;
    dx_type dx;
    double dtval, dtmin, dtmax;
    double t;
    int goodsteps;
    int stepcount;
  public:
    Integrator(const T &stateval, dx_type dx, double dtmin, double dtmax);
    double time() const;
    T& state();
    double dt() const;
    void dt(double h);
    void step();
    int steps() const;
};

template <class T>
Integrator<T>::Integrator(const T &stateval, dx_type dx, double dtmin, double dtmax) : stateval(stateval), dx(dx), dtmin(dtmin), dtmax(dtmax) {
  assert(dtmin>0);
  assert(dtmax>0);
  assert(dtmin<=dtmax);

  dtval=dtmin;
  stepcount=0;
  t=0;
}

template<class T>
double Integrator<T>::time() const { return t; }

template<class T>
T& Integrator<T>::state() { return stateval; }

template<class T>
double Integrator<T>::dt() const { return dtval; }

template<class T>
void Integrator<T>::dt(double h) {
  assert(h>0);
  assert(h>=dtmin);
  assert(h<=dtmax);
  dtval=h;
}

template<class T>
int Integrator<T>::steps() const { return stepcount; }

template<class T>
void Integrator<T>::step() {
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





/**
  The EventIntegrator class supplies a simple Euler's method integrator with an
  adaptive step-size scheme. It also allows for discrete events. It approaches
  these events with an exponentially-decreasing stepsize and moves away from
  them in the same fashion.
*/
template <class T>
class EventIntegrator : public Integrator<T> {
  private:
    CalendarQueue calq;
    bool at_event;
    std::string at_event_name;
  public:
    EventIntegrator(const T &stateval, typename Integrator<T>::dx_type dx, double dtmin, double dtmax) : Integrator<T>(stateval, dx, dtmin, dtmax) {
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
void EventIntegrator<T>::insert_event(double t, const std::string &event, double recur_int=0){
  calq.insert(t,event,recur_int);
}

template<class T>
void EventIntegrator<T>::insert_event(double t, const char event[], double recur_int=0){
  calq.insert(t,event,recur_int);
}

template<class T>
bool EventIntegrator<T>::is_event() const {
  return at_event;
}

template<class T>
std::string EventIntegrator<T>::event() const {
  if( is_event() )
    return calq.current_event();
  else
    return "";
}

template<class T>
void EventIntegrator<T>::step() {
  T e1, e2;

  if(at_event && calq.current_time()==Integrator<T>::t){
    at_event_name=calq.current_event();
    calq.pop();
    return;
  }

  at_event=false;

  ++Integrator<T>::stepcount;

  if(!calq.empty()){
    while(
        Integrator<T>::dtval > Integrator<T>::dtmin
        && Integrator<T>::t + Integrator<T>::dtval > calq.current_time()
    ) {
      Integrator<T>::dtval/=2.;
    }

    if(Integrator<T>::dtval < Integrator<T>::dtmin)
      Integrator<T>::dtval = Integrator<T>::dtmin;

    if(Integrator<T>::dtval == Integrator<T>::dtmin && Integrator<T>::t + Integrator<T>::dtval > calq.current_time()){
      Integrator<T>::dtval=calq.current_time()-Integrator<T>::t;
      Integrator<T>::dx(Integrator<T>::stateval, e1, Integrator<T>::t);
      Integrator<T>::stateval+=e1*Integrator<T>::dtval;
      Integrator<T>::t=calq.current_time();
      Integrator<T>::dtval = Integrator<T>::dtmin;
      at_event=true;
      at_event_name=calq.current_event();
      calq.pop();
      return;
    }
  }    

  Integrator<T>::dx(Integrator<T>::stateval           , e1, Integrator<T>::t);
  Integrator<T>::dx(Integrator<T>::stateval + e1*Integrator<T>::dtval/2, e2, Integrator<T>::t+Integrator<T>::dtval/2.);
  double abs_e1=abs(Integrator<T>::stateval + e1*Integrator<T>::dtval);
  double abs_e2=abs(Integrator<T>::stateval + e1*Integrator<T>::dtval/2 + e2*Integrator<T>::dtval/2);

  if( (std::abs(abs_e1-abs_e2)/(std::abs(abs_e1+abs_e2)/2))>0.05 ){
    Integrator<T>::stateval+=e1*Integrator<T>::dtval/2.;
    Integrator<T>::t+=Integrator<T>::dtval/2.;
    Integrator<T>::dtval/=2.;
    if(Integrator<T>::dtval<Integrator<T>::dtmin)
      Integrator<T>::dtval=Integrator<T>::dtmin;
    Integrator<T>::goodsteps=0;
  } else {
    Integrator<T>::stateval+=e1*Integrator<T>::dtval;
    Integrator<T>::t+=Integrator<T>::dtval;
    ++Integrator<T>::goodsteps;
  }

  if(Integrator<T>::dtval<Integrator<T>::dtmax && Integrator<T>::goodsteps>=8)
    Integrator<T>::dtval*=2;
  if(Integrator<T>::dtval>Integrator<T>::dtmax)
    Integrator<T>::dtval=Integrator<T>::dtmax;
}

}

#endif
