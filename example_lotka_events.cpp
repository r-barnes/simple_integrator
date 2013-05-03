#include <iostream>
#include <array>
#include "arraystate.hpp"
#include "integrator.hpp"
using namespace std;

typedef arraystate<double, 2> state_type;

void df(const state_type &state, state_type &dxdt, double t){
  dxdt[0]= 1.5*state[0]-1*state[0]*state[1];
  dxdt[1]=-3  *state[1]+1*state[0]*state[1];
}

int main(){
  state_type x={{10,4}};
  event_integrator<state_type> stepper(x, df, 0.01, 1e-3);

  stepper.insert_event(10,"large_drought");
  stepper.insert_event(4,"recurring_drought",3);

  cout<<stepper.time()<<" "<<stepper.state()<<endl;
  while(stepper.time()<20){
    if(stepper.is_event()){
      if(stepper.event()=="large_drought")
        stepper.state()[0]*=0.3;
      else if(stepper.event()=="recurring_drought")
        stepper.state()[0]*=0.5;
    }
    stepper.step();
    cout<<stepper.time()<<" "<<stepper.state()<<endl;
  }
}

