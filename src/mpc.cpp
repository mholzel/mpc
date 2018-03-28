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

void print_state(const std::vector<T> &state) {
    for (size_t i = 0; i < state.size(); ++i)
        cout << state[i] + 1e-4 << ", ";
    cout << endl;
}

void MPC::solve(const std::vector<T> &initial_state,
                const std::vector<T> &reference_polynomial,
                const std::vector<T> &previous_controls,
                const T global_psi) {

    typedef CPPAD_TESTVECTOR(
    double) Dvector;

    /*
     * ---------------------------------------------------------------------
     *
     * Problem dimensions
     *
     * ---------------------------------------------------------------------
     */
    const size_t n_states = 4;
    const size_t n_controls = 2;
    const size_t n_vars = N * n_controls;
    const size_t n_constraints = 0;
    assert(initial_state.size() == n_states);

    /*
     * ---------------------------------------------------------------------
     *
     * Time delay correction
     *
     * ---------------------------------------------------------------------
     */
    /* Here we update our time delay estimate if this is not the first time solving the MPC problem */
    T to = Control<T>::throttle(previous_controls.data()); // throttle_optimized
    T so = Control<T>::steering_angle(previous_controls.data());; // steering_optimized
    T vc = State<T>::v(initial_state.data()); // v_current
    T pc = global_psi -
           previous_global_psi; // psi_current. We need to do it like this because psi is continuously reset to 0
    if (initialized) {

        /* The instantaneous time delay is estimated to be the most likely
         * time delay which fits the previous samples.
         *
         * Specifically, in the previous iteration, we were given a
         * previous velocity, psi, throttle, and steering angle
         * (these are values from two time-steps ago relative to the present).
         * The assumption is that the previous throttle and steering angle
         * (from two timesteps ago) are being used until the simulator receives
         * our optimized throttle and steering angle from the last iteration
         * (one timestep ago).
         * It will then use these values until we give it a new command.
         * However, we don't know what the time delay is between when we send a command, and when
         * it will be implemented by the simulator. So we need to estimate it.
         * We also need to estimate the time between consecutive control commands.
         * In summary, our equations look like:
         *
            v_intermediate = v_previous + throttle_previous * time_delay;
            v_current = v_intermediate + throttle_optimized * time_until_next_control;

            psi_intermediate = psi_previous + v_previous * steering_previous / Lf * time_delay;
            psi_current = psi_intermediate + v_intermediate * steering_optimized / Lf * time_until_next_control;

            Solving this system of equations yields two possible solutions:
         */
        T sqrt_term = sqrt(pow(so, 2) * pow(tp, 2) * pow(vc, 2) + pow(sp, 2) * pow(to, 2) * pow(vp, 2) -
                           4 * so * sp * to * tp * pow(vp, 2) - 4 * Lf * pc * so * to * pow(tp, 2) +
                           4 * Lf * pp * so * to * pow(tp, 2) + 2 * so * sp * to * tp * vc * vp);
        T td1 = +(sqrt_term + so * tp * vc - 2 * so * tp * vp + sp * to * vp) / (2 * so * pow(tp, 2));
        T td2 = -(sqrt_term - so * tp * vc + 2 * so * tp * vp - sp * to * vp) / (2 * so * pow(tp, 2));

        /* We only permit time delay estimates in the range [ 0, .5 ].
         * If one of these values is inside that range, update the time delay estimate
         */
        if (td1 >= 0 && td1 <= 0.5) {
            time_delay = current_td_scale * td1 + past_td_scale * time_delay;
            cout << "Valid time delay estimate : " << td1 << endl;
        } else if (td2 >= 0 && td2 <= 0.5) {
            time_delay = current_td_scale * td2 + past_td_scale * time_delay;
            cout << "Valid time delay estimate : " << td2 << endl;
        } else {
            cout << "Both time delay estimates were invalid: " << td1 << ", " << td2 << endl;
        }

    } else {
        initialized = true;
    }

    /* Update the cached values */
    tp = to;
    sp = so;
    vp = vc;
    pp = 0;
    previous_global_psi = global_psi;

    /* Now that we have updated our time delay estimate, we will shift our initial state to
     * value that we think it will be at when the simulator actually receives our command. */
    // TODO
    time_delay = 0.1;
    Cost tmp_cost(Lf, initial_state, reference_polynomial, N, time_delay, reference_velocity);
    std::vector<T> initial_state_after_delay = initial_state;
    tmp_cost.dynamics(initial_state_after_delay.data(),
                      const_cast<T *>(initial_state.data()),
                      const_cast<T *>(previous_controls.data()));

    cout << time_delay << ", Initial state " << endl;
    print_state(initial_state);
    cout << "Initial state after time delay" << endl;
    print_state(initial_state_after_delay);

    /*
     * ---------------------------------------------------------------------
     *
     * Cost function
     *
     * ---------------------------------------------------------------------
     */
    // TODO Change to initial_state_after_delay
    Cost cost(Lf, initial_state_after_delay, reference_polynomial, N, dt, reference_velocity);

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

    std::cout << "cost = " << solution.obj_value << ", ";
    std::cout << "status = " << solution.status << std::endl << std::endl;

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
