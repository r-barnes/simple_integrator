#include <iostream>
#include "integrator.hpp"
using namespace std;

void df(const double &x, double &dxdt, double t){
  dxdt=2*t;
}

int main(){
  double x=0.5;
  integrator<double> stepper(x, df, 0.01, 1e-6);

  while(stepper.steps()<100){
    stepper.step();
    cout<<stepper.time()<<" "<<stepper.state()<<endl;
  }
}
