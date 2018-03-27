#ifndef MPC_H
#define MPC_H

#include <vector>
#include "Eigen/Core"
#include "cost_ic.h"

using namespace std;

template <typename T>
class MPC
{
  public:
    /* Constants used to define the optimization problem */
    const size_t N;
    const T dt;
    const T reference_velocity;

    /* These values will be filled after solving */
    vector<T> x_vals;
    vector<T> y_vals;
    T normalized_steering_angle;
    T normalized_throttle;

    /* Constructor */
    MPC(size_t N,
        T dt,
        T reference_velocity)
            : N(N),
              dt(dt),
              reference_velocity(reference_velocity) {
    }

    virtual ~MPC() {}

    /* Given an initial state and reference polynomial to track, solve determine the optimal trajectory, and return
     * the first actuations */
    void solve(const std::vector<T> &initial_state,
               const std::vector<T> &reference_polynomial,
               const std::vector<T> &previous_controls)
    {

        typedef CPPAD_TESTVECTOR(T)
            Dvector;

        /*
         * ---------------------------------------------------------------------
         *
         * Cost function
         *
         * ---------------------------------------------------------------------
         */
        /* For the purposes of the cost function, the initial state needs to have type AD<T> */
        Cost<T> cost(reference_polynomial, N, dt, reference_velocity);

        /*
         * ---------------------------------------------------------------------
         *
         * Problem dimensions
         *
         * ---------------------------------------------------------------------
         */
        const size_t n_states = 6;
        const size_t n_controls = 2;
        const size_t n_vars = n_states + N * n_controls;
        const size_t n_constraints = n_states;
        assert(initial_state.size() == n_states);

        /*
         * ---------------------------------------------------------------------
         *
         * Initial Guess
         *
         * ---------------------------------------------------------------------
         */
        /* The first portion of the variables holds the initial state. 
         * After that come all of the controls.
         * 
         * We assume that the steering angle is 0 for all time steps,
         * although the throttle remains constant and equal to the previous control value. */
        Dvector vars(n_vars);
        for (size_t i = 0; i < n_states; ++i)
            vars[i] = initial_state[i];

        Control<T> control(vars.data() + n_states);
        for (size_t t = 0; t < N; ++t)
        {
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
        Control<T> control_lowerbounds(vars_lowerbound.data() + n_states);
        Control<T> control_upperbounds(vars_upperbound.data() + n_states);
        for (size_t t = 0; t < N; ++t)
        {
            control_lowerbounds.steering_angle(t) = -0.436332;
            control_upperbounds.steering_angle(t) = 0.436332;
            control_lowerbounds.throttle(t) = -1;
            control_upperbounds.throttle(t) = 1;
        }

        /*
         * ---------------------------------------------------------------------
         *
         * Constraint Bounds
         *
         * ---------------------------------------------------------------------
         */
        /* The constraints just ensure that the initial state is not optimized */
        Dvector constraints_lowerbound(n_constraints);
        Dvector constraints_upperbound(n_constraints);
        for (size_t i = 0; i < n_states; ++i)
        {
            constraints_lowerbound[i] = initial_state[i];
            constraints_upperbound[i] = initial_state[i];
        }

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
        CppAD::ipopt::solve_result<Dvector> solution;

        /* solve */
        CppAD::ipopt::solve<Dvector, Cost<T>>(
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

        T *controls = solution.x.data() + n_states;
        /* The optimization will give us a steering angle in radians. The maximum value that this can be is +/- 25 degrees.
         * After that, we need to convert from that angle to the range [-1,1], where -1 denotes -25 degress and +1 denotes
         * +25 degrees.
         */
        normalized_steering_angle = Control<T>::steering_angle(controls) / 0.436332;
        normalized_throttle = Control<T>::throttle(controls);

        // TODO
        // T *states = solution.x.data();
        // for (size_t i = 0; i < N; ++i)
        // {
        //     x_vals[i] = State<T>::x(states);
        //     y_vals[i] = State<T>::y(states);
        //     states += n_states;
        // }

        const bool verbose = true;
        if (verbose)
        {
            /* Show the steering angles that we want to use */
            Eigen::Map<Eigen::MatrixXd> c(controls, n_controls, N);
            cout << "Steering angles" << endl
                 << c.topRows<1>() / 0.436332 << endl
                 << endl;
            cout << "Throttles" << endl
                 << c.bottomRows<1>() << endl
                 << endl;
        }
    }
};

#endif /* MPC_H */
