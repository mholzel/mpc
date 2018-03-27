#include "mpc.h"
#include "cost.h"
#include "control.h"
#include "state.h"
#include "cost.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include <iostream>

using std::cout;
using std::endl;

void MPC::solve(const std::vector<T> &initial_state,
                const std::vector<T> &reference_polynomial,
                const std::vector<T> &previous_controls) {

    typedef CPPAD_TESTVECTOR(
    double) Dvector;

    /*
     * ---------------------------------------------------------------------
     *
     * Problem dimensions
     *
     * ---------------------------------------------------------------------
     */
    const size_t n_states = 6;
    const size_t n_controls = 2;
    const size_t n_vars = N * n_controls;
    const size_t n_constraints = 0;
    assert(initial_state.size() == n_states);

    /*
     * ---------------------------------------------------------------------
     *
     * Cost function
     *
     * ---------------------------------------------------------------------
     */
    Cost cost(initial_state, reference_polynomial, N, dt, reference_velocity);

    /*
     * ---------------------------------------------------------------------
     *
     * Initial Guess
     *
     * ---------------------------------------------------------------------
     */
    /* We assume that the steering angle is 0 for all time steps,
     * although the throttle remains constant and equal to the previous control value. */
    Dvector vars(n_vars);
    Control<T> control(vars.data());
    for (size_t t = 0; t < N; ++t) {
        control.steering_angle(t) = 0;
        control.throttle(t) = Control<T>::throttle(previous_controls.data());
    }

    /*
     * ---------------------------------------------------------------------
     *
     * Variable bounds
     *
     * ---------------------------------------------------------------------
     */
    /* The control throttles are constrained to the range +/- 1.
     * The control steering angles are constrained to the range +/- 25 degrees */
    Dvector vars_lowerbound(n_vars);
    Dvector vars_upperbound(n_vars);
    Control<T> control_lowerbounds(vars_lowerbound.data());
    Control<T> control_upperbounds(vars_upperbound.data());
    for (size_t t = 0; t < N; ++t) {
        control_lowerbounds.steering_angle(t) = -0.436332;
        control_upperbounds.steering_angle(t) = 0.436332;
        control_lowerbounds.throttle(t) = -.99;
        control_upperbounds.throttle(t) = .99;
    }

    /*
     * ---------------------------------------------------------------------
     *
     * Constraint Bounds
     *
     * ---------------------------------------------------------------------
     */
    /* We have no constraints */
    Dvector constraints_lowerbound(n_constraints);
    Dvector constraints_upperbound(n_constraints);

    /*
     * ---------------------------------------------------------------------
     *
     * IPOPT Options
     *
     * ---------------------------------------------------------------------
     */
    std::string options;
    options += "Integer print_level  0\n";
    // NOTE: Setting sparse to true allows the solver to take advantage
    // of sparse routines, this makes the computation MUCH FASTER. If you
    // can uncomment 1 of these and see if it makes a difference or not but
    // if you uncomment both the computation time should go up in orders of
    // magnitude.
    options += "Sparse  true        forward\n";
    options += "Sparse  true        reverse\n";
    // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
    // Change this as you see fit.
    //    options += "Numeric max_cpu_time          0.5\n";

    /*
     * ---------------------------------------------------------------------
     *
     * Solve
     *
     * ---------------------------------------------------------------------
     */
    /* Allocate stack space for the solution */
    CppAD::ipopt::solve_result <Dvector> solution;

    /* solve */
    CppAD::ipopt::solve<Dvector, Cost>(
            options,
            vars,
            vars_lowerbound,
            vars_upperbound,
            constraints_lowerbound,
            constraints_upperbound,
            cost,
            solution);

    /*
     * ---------------------------------------------------------------------
     *
     * Diagnostics and clean up
     *
     * ---------------------------------------------------------------------
     */
    bool ok = true;
    ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

    std::cout << "cost = " << solution.obj_value << ", status = " << solution.status << std::endl
              << std::endl;

    /* Now set the fields in the MPC controller so that they can be retrieved elsewhere. */
    x_vals.resize(N);
    y_vals.resize(N);

    /* The optimization will give us a steering angle in radians. The maximum value that this can be is +/- 25 degrees.
     * After that, we need to convert from that angle to the range [-1,1], where -1 denotes -25 degress and +1 denotes
     * +25 degrees.
     */
    T *t = solution.x.data();
    normalized_steering_angle = Control<T>::steering_angle(t) / 0.436332;
    normalized_throttle = Control<T>::throttle(t);

    /* Next, we want to calculate the predicted trajectory so that
     * we can display it */
    std::vector<T> states = cost.calculateTrajectory(t);
    T *state_data = states.data();
    for (size_t i = 0; i < N; ++i) {
        x_vals[i] = State<T>::x(state_data);
        y_vals[i] = State<T>::y(state_data);
        state_data += n_states;
    }

    const bool verbose = true;
    if (verbose) {
        /* Show the steering angles that we want to use */
        Eigen::Map<Eigen::MatrixXd> control(t, n_controls, N);
        cout << "Steering angles" << endl
             << control.topRows<1>() / 0.436332 << endl
             << endl;
        cout << "Throttles" << endl
             << control.bottomRows<1>() << endl
             << endl;
    }
}
