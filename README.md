SIMPLE INTEGRATOR
=================

How would you perform a numerical integration of an equation with the form

    y'(t)=f(t, y(t))

There are many integrators in the world. Many of the emphasise accuracy and use
high-order methods to estimate and minimize error. This kind of accuracy is
expensive. For systems that are generally well-behaved where the goal is
understanding of dynamics, rather than exact predictions, accuracy is no longer
paramount.

This C++11 library implements a simple adaptive step-size integrator based on
[Euler's method](https://en.wikipedia.org/wiki/Euler_method). The method is
both fast and simple, allowing for the quick simulation of many-dimensional
systems.

Additionally, a special class is provided which allows for discrete-time events
which instantaneously alter the state of a system. The integrator approaches
such events cautiously, using an exponentially-decreasing step-size. It
withdraws with equal caution.

The file **integrator.hpp** implements the integrator.

Examples
========

Two examples of the integrator in action are provided.

The first example **example\_lotka.cpp** simulates a Lotka-Volterra
predator-prey model. This code:

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
      si_lib::Integrator<state_type> stepper(x, df, 1e-3, 0.01);

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

Produces this result:

![Lotka-Volterra predator-prey model output](readme_imgs/lotka.png)

The second example **example\_lotka\_events.cpp** simulates a Lotka-Volterra
predator-prey model with discrete events. Here one of the discrete events is a
sudden single large drought which reduces one of the species' populations to
about one-third of its pre-drought size. The other event is a smaller drought
which reduces the same population to half its pre-drought size every 3 years
beginning on the 4th year of the simulation.

This code:

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

        //Has the time step ended on an event?
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

Produces this result:

![Lotka-Volterra predator-prey model output with discrete events]
(readme_imgs/lotka_events.png)

Building
========

To include the integrator in your own project, just add **integrator.hpp** to
your project and include it in your implementation, as demonstrated in the
examples.

To build, run

    make

and **Makefile** will be run. Or compile with a C++11 compiler.

This library is known to work with G++ 4.7.3.
