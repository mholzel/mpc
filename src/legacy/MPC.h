#ifndef MPC_HEADER
#define MPC_HEADER

#include <vector>
#include "Eigen/Core"
#include "cost.h"
#include <cppad/ipopt/solve.hpp>

#include <iostream>

using std::cout;
using std::endl;

template<typename T>
class MPC {
public:

    const size_t N;
    const double dt;
    const double reference_velocity;

    std::vector<double> x_vals;
    std::vector<double> y_vals;
    double normalized_steering_angle;
    double normalized_throttle;

    MPC(size_t N,
        double dt,
        double reference_velocity)
            : N(N), dt(dt), reference_velocity(reference_velocity) {
    }

    virtual ~MPC() {}

    /* Given an initial state and reference polynomial to track, solve determine the optimal trajectory, and return
     * the first actuations */
    void Solve(const std::vector<T> &initial_state,
               const std::vector<T> &reference_polynomial,
               const std::vector<T> &previous_controls) {


        typedef CPPAD_TESTVECTOR(
        T) Dvector;

        /* This object will compute the cost function and constraints for our specific problem. */
        Cost<T> cost(initial_state, reference_polynomial, N, dt, reference_velocity);

        /* Specify the problem dimensions */
        size_t n_states = 6;
        size_t n_controls = 2;
        size_t n_vars = N * n_states + (N - 1) * n_controls;
        size_t n_constraints = N * n_states;


        /*---------------------------------------------------------------------*/
        /*
         * Generate the initial guess for the problem.
         * To do this, we will set the initial state equal to the passed-in value.
         *
         * Then we assume that the steering angle is 0 for all time steps,
         * although the throttle remains constant and equal to the previous control value.
         *
         * Finally, we compute the values of the other states assuming that these control values are correct.
         */
        /* The first step is to make sure all of the initial variables are zero. */
        Dvector vars(n_vars);
        for (size_t i = 0; i < n_vars; i++)
            vars[i] = 0.0;

        /* Set the desired control values */
        Control<T> control(vars.data() + N * n_states);
        for (size_t t = 0; t < N; ++t)
            control.throttle(t) = Control<T>::throttle(previous_controls.data());

        /* Now set the initial state */
        State<T> state(vars.data());
        state.x(0) = State<T>::x(initial_state.data());
        state.y(0) = State<T>::y(initial_state.data());
        state.psi(0) = State<T>::psi(initial_state.data());
        state.v(0) = State<T>::v(initial_state.data());
        state.crosstrack_error(0) = State<T>::crosstrack_error(initial_state.data());
        state.psi_error(0) = State<T>::psi_error(initial_state.data());

        /* Finally, set the remaining states using the previously-specified control values */
        for (size_t t = 1; t < N; ++t)
            cost.dynamics(state.at(t), state.at(t - 1), control.at(t - 1));

        /*---------------------------------------------------------------------*/
        /* Next, generate the upper and lower variable bounds.
         *
         * For the state, these are essentially use +/- infinity.
         *
         * The control throttles are constrained to the range +/- 1.
         * The control steering angles are constrained to the range +/- 25 degrees
         *
         * Note that we use the constraints to ensure that the initial state
         * does not change.
         */
        Dvector vars_lowerbound(n_vars);
        Dvector vars_upperbound(n_vars);

        /* First, set all of the bounds to +/- infinity */
        for (size_t i = 0; i < n_vars; i++) {
            vars_lowerbound[i] = -1.0e19;
            vars_upperbound[i] = 1.0e19;
        }

        /* Next, constrain the steering angles to +/- 25 degrees (in radians),
         * and the throttles to +/-1.
         *
         * TODO: Might want to use something smaller for the angles to make sure that the controller is not erratic. */
        T *control_lowerbounds = vars_lowerbound.data() + N * n_states;
        T *control_upperbounds = vars_upperbound.data() + N * n_states;
        for (size_t t = 0; t < N - 1; ++t) {

            Control<T>::steering_angle(control_lowerbounds) = -0.436332;
            Control<T>::steering_angle(control_upperbounds) = 0.436332;
            Control<T>::throttle(control_lowerbounds) = -1;
            Control<T>::throttle(control_upperbounds) = 1;

            control_lowerbounds += n_controls;
            control_upperbounds += n_controls;
        }

        /* Finally, create the constraint bounds.
         * All of these are 0 (to signify equality constraints),
         * except for the initial condition constraints.
         */
        Dvector constraints_lowerbound(n_constraints);
        Dvector constraints_upperbound(n_constraints);
        for (size_t i = 0; i < n_constraints; ++i) {
            constraints_lowerbound[i] = 0;
            constraints_upperbound[i] = 0;
        }

        T *constraint_lowerbounds = constraints_lowerbound.data();
        State<T>::x(constraint_lowerbounds) = State<T>::x(initial_state.data());
        State<T>::y(constraint_lowerbounds) = State<T>::y(initial_state.data());
        State<T>::psi(constraint_lowerbounds) = State<T>::psi(initial_state.data());
        State<T>::v(constraint_lowerbounds) = State<T>::v(initial_state.data());
        State<T>::crosstrack_error(constraint_lowerbounds) = State<T>::crosstrack_error(initial_state.data());
        State<T>::psi_error(constraint_lowerbounds) = State<T>::psi_error(initial_state.data());

        T *constraint_upperbounds = constraints_upperbound.data();
        State<T>::x(constraint_upperbounds) = State<T>::x(initial_state.data());
        State<T>::y(constraint_upperbounds) = State<T>::y(initial_state.data());
        State<T>::psi(constraint_upperbounds) = State<T>::psi(initial_state.data());
        State<T>::v(constraint_upperbounds) = State<T>::v(initial_state.data());
        State<T>::crosstrack_error(constraint_upperbounds) = State<T>::crosstrack_error(initial_state.data());
        State<T>::psi_error(constraint_upperbounds) = State<T>::psi_error(initial_state.data());

        // options
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
//        options += "Numeric max_cpu_time          0.5\n";

        // place to return solution
        CppAD::ipopt::solve_result <Dvector> solution;

        // solve the problem
        CppAD::ipopt::solve<Dvector, Cost<T>>(
                options,
                vars,
                vars_lowerbound,
                vars_upperbound,
                constraints_lowerbound,
                constraints_upperbound,
                cost,
                solution);

        //
        // Compute the solution and show some diagnostics
        //
        bool ok = true;
        ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

        cout << "cost = " << solution.obj_value << ", status = " << solution.status << endl << endl;

        /* Now set the fields in the MPC controller so that they can be retrieved elsewhere. */
        x_vals.resize(N);
        y_vals.resize(N);

        T *controls = solution.x.data() + N * n_states;
        /* The optimization will give us a steering angle in radians. The maximum value that this can be is +/- 25 degrees.
         * After that, we need to convert from that angle to the range [-1,1], where -1 denotes -25 degress and +1 denotes
         * +25 degrees.
         */
        normalized_steering_angle = Control<T>::steering_angle(controls) / 0.436332;
        normalized_throttle = Control<T>::throttle(controls);

        T *states = solution.x.data();
        for (size_t i = 0; i < N; ++i) {
            x_vals[i] = State<T>::x(states);
            y_vals[i] = State<T>::y(states);
            states += n_states;
        }

        const bool verbose = true;
        if (verbose) {
            /* Show the steering angles that we want to use */
            Eigen::Map <Eigen::MatrixXd> c(controls, n_controls, N);
            cout << "Steering angles" << endl << c.topRows<1>() / 0.436332 << endl << endl;
            cout << "Throttles" << endl << c.bottomRows<1>() << endl << endl;
        }
    }

};

#endif /* MPC_HEADER */
