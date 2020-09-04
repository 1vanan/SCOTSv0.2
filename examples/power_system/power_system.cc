/*
 * power_system.cc
 *
 *  created: Sep 2020
 *  author: Ivan Fedotov
 *
 */

#include <iostream>
#include <array>

/* SCOTS header */
#include "scots.hh"
/* ode solver */
#include "RungeKutta4.hh"

/* time profiling */
#include "TicToc.hh"
/* memory profiling */
#include <sys/time.h>
#include <sys/resource.h>
struct rusage usage;

/* state space dim */
const int state_dim=3;
/* input space dim */
const int input_dim=2;
/* sampling time */
const double tau = 0.25;

/* data types of the state space elements and input
 * space elements used in uniform grid and ode solver */
using state_type = std::array<double,state_dim>;
using input_type = std::array<double,input_dim>;

// Values are taken for the average inductivity and resistance.
double L = 20 * 0.000001; // Henry
double R = 1000; // Ohm

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

/* we integrate the growth bound by 0.25 sec (the result is stored in r)  */
auto radius_post = [] (state_type &r, const state_type &, const input_type &u) {
    /* the ode for the growth bound */
    auto rhs =[] (state_type& rr,  const state_type &r, const input_type &u) {
        /* lipschitz matrix in this case can be computed manually.
         * ||  df_0/dx_0          df_0/dx_1  ||
         * ||  df_1/dx_0          df_1/dx_1  ||
         * */
        double L[1][1];
        L[0][0]=R/L;
        L[0][1]=R/L;
        L[1][0]=R/L;
        L[1][1]=R/L;
        /* to account for input disturbances TODO: why input? disturbance here is for inductivity*/
        const state_type w={{.108,0.002}};
        rr[0] = L[0][0]*r[0]+L[0][1]*r[1]+w[0];
        rr[1] = L[1][0]*r[0]+L[1][1]*r[1]+w[1];
    };
    /* use 10 intermediate steps */
    scots::runge_kutta_fixed4(rhs,r,u,state_dim,tau,5);
};

int main() {
    /* to measure time */
    TicToc tt;

    /* construct grid for the state space */
    /* setup the workspace of the synthesis problem and the uniform grid */
    /* grid node distance diameter */
    /* optimized values computed according to doi: 10.1109/CDC.2015.7403185 */
    // TODO: Come up which amount for diameter is better
    state_type s_eta={{25.0/362,3*M_PI/180/66}};
    /* lower bounds of the hyper rectangle */
    state_type s_lb={{13,13}};
    /* upper bounds of the hyper rectangle */
    state_type s_ub={{19,19}};
    scots::UniformGrid ss(state_dim,s_lb,s_ub,s_eta);
    std::cout << "Uniform grid details:" << std::endl;
    ss.print_info();

    /* construct grid for the input space */
    /* lower bounds of the hyper rectangle */
    input_type i_lb={{0,0}};
    /* upper bounds of the hyper rectangle.
     * Max voltage is 220 V and disturbance can be up to 10%.*/
    input_type i_ub={{242,242}};
    /* grid node distance diameter */
    // TODO: Come up which amount for diameter is better
    input_type i_eta={{32000,9.0/8.0*M_PI/180}};
    scots::UniformGrid is(input_dim,i_lb,i_ub,i_eta);
    is.print_info();

    /* transition function of symbolic model */
    scots::TransitionFunction tf;

    /* setup object to compute the transition function */
    scots::Abstraction<state_type,input_type> abs(ss,is);
    /* measurement disturbances  */
    state_type z={{0.2,0.2}};
    abs.set_measurement_error_bound(z);

    std::cout << "Computing the transition function: " << std::endl;
    tt.tic();
    abs.compute_gb(tf,power_system_post,radius_post);
    tt.toc();
    if(!getrusage(RUSAGE_SELF, &usage))
        std::cout << "Memory per transition: " << usage.ru_maxrss/(double)tf.get_no_transitions() << std::endl;
    std::cout << "Number of transitions: " << tf.get_no_transitions() << std::endl;

    /* define target set */
    // TODO: why this is different current than in @power_system?
    auto target = [&s_eta, &z, &ss](const scots::abs_type& abs_state) {
        state_type t_lb = {{63,-3*M_PI/180,0}};
        state_type t_ub = {{75,0,2.5}};
        state_type c_lb;
        state_type c_ub;
        /* center of cell associated with abs_state is stored in x */
        state_type x;
        ss.itox(abs_state,x);
        /* hyper-interval of the quantizer symbol with perturbation */
        for(int i=0; i<state_dim; i++) {
            c_lb[i] = x[i]-s_eta[i]/2.0-z[i];
            c_ub[i] = x[i]+s_eta[i]/2.0+z[i];
        }
        if( t_lb[0]<=c_lb[0] && c_ub[0]<=t_ub[0] &&
            t_lb[1]<=c_lb[1] && c_ub[1]<=t_ub[1] &&
            t_lb[2]<=c_lb[2] && c_ub[2]<=t_ub[2]) {
            if(-0.91<=(x[0]*std::sin(x[1])-s_eta[0]/2.0-z[0]-(c_ub[0])*(s_eta[1]/2.0-z[1]))) {
                return true;
            }
        }
        return false;
    };
    /* write grid point IDs with uniform grid information to file */
    write_to_file(ss,target,"target");

    std::cout << "\nSynthesis: " << std::endl;
    tt.tic();
    scots::WinningDomain win=scots::solve_reachability_game(tf,target);
    tt.toc();
    std::cout << "Winning domain size: " << win.get_size() << std::endl;

    std::cout << "\nWrite controller to controller.scs \n";
    if(write_to_file(scots::StaticController(ss,is,std::move(win)),"controller"))
        std::cout << "Done. \n";

    return 1;
}

