/*
 * simulate.cc
 *
 *  created: Sep 2020
 *   author: Ivan Fedotov
 */

#include <iostream>
#include <array>
#include <cmath>

/* SCOTS header */
#include "scots.hh"
/* ode solver */
#include "RungeKutta4.hh"

/* state space dim */
const int state_dim=2;
/* input space dim */
const int input_dim=2;

/* sampling time */
const double tau = 0.25;

using state_type = std::array<double,state_dim>;
using input_type = std::array<double,input_dim>;

/* we integrate the power system ode by 0.25 sec (the result is stored in x)  */
auto power_system_post = [] (state_type &x, const input_type &u) {
    /* the ode describing the aircraft */
    auto rhs =[] (state_type& xx,  const state_type &x, const input_type &u) {
        xx[0] = 1/L * (u[0] - (x[0] + x[1])*R);
        xx[1] = 1/L * (u[1] - (x[0] + x[1])*R);
    };
    /* use 10 intermediate steps */
    scots::runge_kutta_fixed4(rhs,x,u,state_dim,tau,5);
};

int main() {

  /* Overall current inside some bounds. */
  auto target = [](const state_type& x) {
      double const MAX_CURRENT = 17.6; //amps
      double const MIN_CURRENT = 14.4;

      return MIN_CURRENT <= x[0] + x[1] <= MAX_CURRENT;
  };

  /* read controller from file */
  scots::StaticController con;
  if(!read_from_file(con,"controller")) {
    std::cout << "Could not read controller from controller.scs\n";
    return 0;
  }
  
  std::cout << "\nSimulation:\n " << std::endl;

  // TODO: is it initial position for the state space?
  state_type x={{81, -1*M_PI/180, 55}};
  while(1) {
    std::vector<input_type> u = con.get_control<state_type,input_type>(x);
    std::cout << x[0] <<  " "  << x[1] << " " << x[2] << "\n";
    //std::cout << u[0][0] <<  " "  << u[0][1] << "\n";
    aircraft_post(x,u[0]);
    if(target(x)) {
      std::cout << "Arrived: " << x[0] <<  " "  << x[1] << " " << x[2] << std::endl;
      break;
    }
  }

  return 1;
}
