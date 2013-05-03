/**
@file example_lotka_events.cpp
Demonstrates how the integrator could be used to implement a Lotka-Volterra
predator-prey model.
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
  si_lib::Integrator<state_type> stepper(x, df, 0.01, 1e-3);

  //Display the initial state
  cout<<stepper.time()<<" "<<stepper.state()<<endl;

  //Run the integrator for 1000 time steps
  while(stepper.steps()<1000){
    //Progress the integrator by one step
    stepper.step();

    //Print out the current state
    cout<<stepper.time()<<" "<<stepper.state()<<endl;
  }
}
