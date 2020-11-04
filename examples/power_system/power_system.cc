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
const int state_dim=2;
/* input space dim */
const int input_dim=2;
/* sampling time */
const double tau = 0.25;

/* data types of the state space elements and input
 * space elements used in uniform grid and ode solver.
 * 1 - value, 2 - to which dimension it corresponds
 */
using state_type = std::array<double,state_dim>;
using input_type = std::array<double,input_dim>;


/* we integrate the power system ode by 0.25 sec (the result is stored in x)  */
auto power_system_post = [] (state_type &x, const input_type &u) {
    const char *path="./resources/nominal.txt";
    std::ofstream outputFile; //open in constructor
    outputFile.open(path, std::ios_base::app);
    /* the ode describing the power system */
    auto rhs =[] (state_type& xx,  const state_type &x, const input_type &u) {
        // Values are taken for the average inductivity and resistance.
        double L = 5; // Henry
        double R = 1; // Ohm

        // f1, f2 definition
        xx[0] = 1/L * (u[0] - (x[0] + x[1])*R);
        xx[1] = 1/L * (u[1] - (x[0] + x[1])*R);
    };

    outputFile << "Writing nominal result for state:";
    for (auto state : x)
        outputFile << state << std::endl;

    outputFile << "And input:";
    for (auto input : u)
        outputFile << input << "  ";
      outputFile << std::endl;


    /* use 5 intermediate steps */
    scots::runge_kutta_fixed4(rhs,x,u,state_dim,tau,5);

        for (auto element : x)
        outputFile << "Result of the reachable set: " << element << std::endl;
    outputFile << std::endl;
};

// TODO: check this in case of error. May be matrix is unstable and gb will be countliniously increase
/* we integrate the growth bound by 0.25 sec (the result is stored in r)  */
auto radius_post = [] (state_type &r, const state_type &, const input_type &u) {
    std::string path="./resources/deviation.txt";
    std::ofstream outputFile; //open in constructor
    outputFile.open(path, std::ios_base::app);

    /* the ode for the growth bound */
    auto rhs =[] (state_type& rr,  const state_type &r, const input_type &u) {
        // Values are taken for the average inductivity and resistance.
        double L = 5; // Henry
        double R = 1; // Ohm
        /* Lipschitz matrix in this case can be computed manually.
         * ||  df_0/dx_0          df_0/dx_1  ||
         * ||  df_1/dx_0          df_1/dx_1  ||
         * */
        double Lp[2][2];
        Lp[0][0]= -R/L;
        Lp[0][1]= R/L;
        Lp[1][0]= R/L;
        Lp[1][1]= -R/L;
        /* to account for input disturbances */
        const state_type w={{0.5,0.5}}; //zonotope square ([0.5 0.5 0], [0.5 0 0.5])
        rr[0] = Lp[0][0]*r[0]+Lp[0][1]*r[1] + w[1];
        rr[1] = Lp[1][0]*r[0]+Lp[1][1]*r[1] + w[2];
    };

    /* logging the reachability set result with deviation.
     * The overall reachable set is: nominal result minus deviation result.
     */
    outputFile << "Writing the deviation result for the state:   ";
    for (auto state : r)
        outputFile << state  << "  ";
    outputFile << std::endl;

    outputFile << "And input:   ";
    for (auto input : u)
        outputFile << input << "  ";
    outputFile << std::endl;

    /* use 10 intermediate steps. Value of the growth bound in time \tau.
     * Result is stored in r, it gives deviation from the center in the final point.
     * Zonotope: c (nominal solution) + g (deviations).
     * */
    scots::runge_kutta_fixed4(rhs,r,u,state_dim,tau,5);

    outputFile << "Result of the deviation reachable set:   ";
    for (auto element : r)
        outputFile << element  << "  ";
    outputFile << std::endl;
};

// TODO: change disturbance in Latex (linear/non-linear).
int main() {
    std::string path_cora="./resources/reachability_res.json";
    std::string path_compare="./resources/compare.json";

    /* to measure time */
    TicToc tt;

    /* construct grid for the state space */
    /* setup the workspace of the synthesis problem and the uniform grid */
    /* grid node distance diameter */
    /* optimized values computed according to doi: 10.1109/CDC.2015.7403185
     * Node is less, precision is bigger, more probability to find the winning domain.*/
    // TODO: Come up which amount for diameter is better
    state_type s_eta={{0.5,0.5}};
    /* lower bounds of the hyper rectangle */
    state_type s_lb={10,10};
    /* upper bounds of the hyper rectangle */
    state_type s_ub={{11,11}};
    scots::UniformGrid ss(state_dim,s_lb,s_ub,s_eta);
    std::cout << "Uniform grid details:" << std::endl;
    ss.print_info();

    /* construct grid for the input space */
    /* lower bounds of the hyper rectangle */
    input_type i_lb={{205,205}};
    /* upper bounds of the hyper rectangle.
     * Max voltage is 220 V and disturbance can be up to 10%.*/
    input_type i_ub={{206,206}};
    /* grid node distance diameter */
    // TODO: Come up which amount for diameter is better
    input_type i_eta={{0.5,0.5}};
    scots::UniformGrid is(input_dim,i_lb,i_ub,i_eta);
    is.print_info();

    /* transition function of symbolic model */
    scots::TransitionFunction tf;

    /* setup object to compute the transition function */
    scots::Abstraction<state_type,input_type> abs(ss,is);
    /* Noise of the sensors. Disturbance in the state vars, not in the derivatives. */
    state_type z={{0,0}};
    abs.set_measurement_error_bound(z);

    std::cout << "Computing the transition function: " << std::endl;
    tt.tic();
//    /* compute tf with scots */
//    abs.compute_gb(tf, power_system_post,radius_post);
//    /* compute tf with cora from @path_cora */
//    abs.compute_transitions(tf, path_cora, power_system_post,radius_post);
//    /* compute tf with scots and compare with cora. result in @path_compare*/
    abs.compute_compare_gb(tf,path_cora, path_compare, power_system_post,radius_post);

    tt.toc();
    if(!getrusage(RUSAGE_SELF, &usage))
        std::cout << "Memory per transition: " << usage.ru_maxrss/(double)tf.get_no_transitions() << std::endl;
    std::cout << "Number of transitions: " << tf.get_no_transitions() << std::endl;

    /* define target set */
    auto target = [&is,&s_eta](const scots::abs_type& idx) {
        state_type x;
        is.itox(idx,x);

        double const MAX_CURRENT = 0; //amps
        double const MIN_CURRENT = 100;

        return MIN_CURRENT <= x[0] + x[1] - s_eta[0] && x[0] + x[1] + s_eta[0] <= MAX_CURRENT;
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

