/**
@file example_lotka_events.cpp
Demonstrates how the integrator could be used to implement a Lotka-Volterra
predator-prey model with discrete events.

Here one of the discrete events is a sudden single large drought which reduces
one of the species' populations to about one-third of its pre-drought size. The
other event is a smaller drought which reduces the same population to half its
pre-drought size every 3 years beginning on the 4th year of the simulation.
*/

#include <iostream>
#include "integrator.hpp"
using namespace std;

//Define state_type to be an double array of size 2
typedef si_lib::ArrayState<double, 2> state_type;

//First derivate of the state at a given time
void df(const state_type &state, state_type &dxdt, double t){
  dxdt[0]= 1.5*state[0]-1*state[0]*state[1];
  dxdt[1]=-3  *state[1]+1*state[0]*state[1];
}

int main(){
  //Set the initial state
  state_type x={{10,4}};

  //Create an EventIntegrator object with initial state **x**, first-derivative
  //**df**, minimum step-size 1e-3, and maximum step-size 0.01
  si_lib::EventIntegrator<state_type> stepper(x, df, 1e-3, 0.01);

  //Insert an event which occurs once at t=10
  stepper.insert_event(10,"large_drought");

  //Insert an event which occurs for the first time at t=4 and every 3 time
  //units afterwards
  stepper.insert_event(4,"recurring_drought",3);

  //Display the initial state
  cout<<stepper.time()<<" "<<stepper.state()<<endl;

  //Run the integrator until t>=20
  while(stepper.time()<20){

    //Has the time step has ended on an event
    if(stepper.is_event()){
      if(stepper.event()=="large_drought") //Look up the event
        stepper.state()[0]*=0.3;           //Perform appropriate actions
      else if(stepper.event()=="recurring_drought")
        stepper.state()[0]*=0.5;
    }

    //Progress the integrator by one step
    stepper.step();

    //Print out the current state
    cout<<stepper.time()<<" "<<stepper.state()<<endl;
  }
}

