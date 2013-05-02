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
				DiscreteEvent (double t, const char event[])       : t(t), event(event) {}
				DiscreteEvent (double t, const std::string &event) : t(t), event(event) {}
        bool operator< (const DiscreteEvent& a) const {
          return t>a.t;
        }
		};

		typedef std::priority_queue<DiscreteEvent, std::vector<DiscreteEvent> > deq_type;
		deq_type deq;
	public:
    ///Schedule \a event at \a t
		void insert(double t, const char event[]){
      //assert(empty() || t>=current_time());
			deq.push(DiscreteEvent(t,event));
		}
    ///Schedule \a event at \a t
		void insert(double t, const std::string &event){
      //assert(empty() || t>=current_time());
			deq.push(DiscreteEvent(t,event));
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
			if(!deq.empty())
				deq.pop();
		}
    ///Take the top element of the queue and insert a copy of it \a dt
    ///into the future
		double reschedule_top(double dt) {
			deq.push(DiscreteEvent(current_time()+dt, current_event()));
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
