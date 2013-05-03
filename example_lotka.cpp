#include <iostream>
#include <array>
#include "integrator.hpp"
using namespace std;

typedef si_lib::arraystate<double, 2> state_type;

void df(const state_type &state, state_type &dxdt, double t){
  dxdt[0]= 1.5*state[0]-1*state[0]*state[1];
  dxdt[1]=-3  *state[1]+1*state[0]*state[1];
}

int main(){
  state_type x={{10,4}};
  si_lib::integrator<state_type> stepper(x, df, 0.01, 1e-3);

  cout<<stepper.time()<<" "<<stepper.state()<<endl;
  while(stepper.steps()<1000){
    stepper.step();
    cout<<stepper.time()<<" "<<stepper.state()<<endl;
  }
}
