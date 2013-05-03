#ifndef _calendar_queue_header
#define _calendar_queue_header
#include <vector>
#include <string>
#include <queue>
#include <cassert>

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

#endif
