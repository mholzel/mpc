#include "MPC_withstates.h"
#include <math.h>
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen/Core"
#include "Eigen/QR"
#include "state.h"
#include "control.h"
#include <iostream>

using CppAD::AD;
using Eigen::Dynamic;
using Eigen::Map;
using Eigen::Matrix;
using std::cout;
using std::endl;

/* Define the datatype of all of the CppAD variables */
using T = double;

class FG_eval
{

  public:
    using ADvector =
        CPPAD_TESTVECTOR(AD<T>);

    /* Constants */
    const T Lf = 2.67;
    const Eigen::Matrix<T, Eigen::Dynamic, 1> initial_state;
    const Eigen::Matrix<T, Eigen::Dynamic, 1> coeffs;
    const size_t N;
    const T dt;
    const T reference_velocity;

    FG_eval(Eigen::Matrix<T, Eigen::Dynamic, 1> initial_state,
            Eigen::Matrix<T, Eigen::Dynamic, 1> coeffs,
            size_t N,
            T dt,
            T reference_velocity)
        : coeffs(coeffs),
          N(N),
          dt(dt),
          reference_velocity(reference_velocity)
    {
    }

    /** This function contains our model of the dynamics, that is,
     *
     * State(t) = dynamics( State(t0), Control(t0) )
     *
     */
    template <typename V>
    void dynamics(V *t, V *t0, V *u_t0)
    {

        /* Evaluate the polynomial at the previous value of x */
        V f0 = 0.0;
        for (int i = 0; i < coeffs.size(); i++)
        {
            f0 += coeffs[i] * pow(State<V>::x(t0), i);
        }

        // TODO. No need to keep recalculating
        V psides0 = CppAD::atan(coeffs[1]);

        /* Now set the values in the state  */
        State<V>::x(t) = State<V>::x(t0) + State<V>::v(t0) * CppAD::cos(State<V>::psi(t0)) * dt;
        State<V>::y(t) = State<V>::y(t0) + State<V>::v(t0) * CppAD::sin(State<V>::psi(t0)) * dt;
        State<V>::psi(t) = State<V>::psi(t0) + State<V>::v(t0) * Control<V>::steering_angle(u_t0) / Lf * dt;
        State<V>::v(t) = State<V>::v(t0) + Control<V>::throttle(u_t0) * dt;
        State<V>::crosstrack_error(t) = (f0 - State<V>::y(t0)) + (State<V>::v(t0) * CppAD::sin(State<V>::psi_error(t0)) * dt);
        State<V>::psi_error(t) = (State<V>::psi(t0) - psides0) + (State<V>::v(t0) * Control<V>::steering_angle(u_t0) / Lf * dt);
    }

    // `fg` is a vector containing the cost and constraints.
    // `vars` is a vector containing the variable values (state & actuators).
    void operator()(ADvector &constraints, ADvector &vars)
    {

        State<AD<T>> constraint(constraints.data() + 1);
        State<AD<T>> state(vars.data());
        Control<AD<T>> control(vars.data() + 6 * N);

        /* First, we calculate the cost. This value should be put in the first entry of fg_vector */
        AD<T> &cost = constraints[0];

        /* The part of the cost based on the reference state. */
        for (size_t t = 0; t < N; ++t)
        {
            cost += CppAD::pow(state.crosstrack_error(t), 2);
            cost += CppAD::pow(state.psi_error(t), 2);
            cost += CppAD::pow(state.v(t) - reference_velocity, 2);
        }

        /* Minimize the use of actuators. */
        for (size_t t = 0; t < N - 1; ++t)
        {
            cost += CppAD::pow(control.steering_angle(t), 2);
            cost += CppAD::pow(control.throttle(t), 2);
        }

        /* Minimize the value gap between sequential actuations. */
        for (size_t t = 0; t < N - 2; ++t)
        {
            cost += CppAD::pow(control.steering_angle(t + 1) - control.steering_angle(t), 2);
            cost += CppAD::pow(control.throttle(t + 1) - control.throttle(t), 2);
        }

        /* Next, we calculate the constraints */
        /* The first constraint is just the initial values of the state. */
        constraint.x(0) = state.x(0);
        constraint.y(0) = state.y(0);
        constraint.psi(0) = state.psi(0);
        constraint.v(0) = state.v(0);
        constraint.crosstrack_error(0) = state.crosstrack_error(0);
        constraint.psi_error(0) = state.psi_error(0);

        /* The rest of the constraints just check that the dynamics are satisfied. */
        for (size_t t = 1; t < N; t++)
        {

            /* First, compute the predicted states at time t */
            dynamics(constraint.at(t), state.at(t - 1), control.at(t - 1));

            /* Now, subtract off the value of the actual variables. The constraints should be 0 after this step. */
            constraint.x(t) -= state.x(t);
            constraint.y(t) -= state.y(t);
            constraint.psi(t) -= state.psi(t);
            constraint.v(t) -= state.v(t);
            constraint.crosstrack_error(t) -= state.crosstrack_error(t);
            constraint.psi_error(t) -= state.psi_error(t);
        }
    }
};

template <typename T>
void MPC<T>::Solve(Eigen::Matrix<T, Eigen::Dynamic, 1> x0,
                   Eigen::Matrix<T, Eigen::Dynamic, 1> coeffs,
                   Eigen::Matrix<T, Eigen::Dynamic, 1> previous_controls)
{

    typedef CPPAD_TESTVECTOR(T)
        Dvector;

    /* This object will compute the cost function and constraints for our specific problem. */
    FG_eval fg(x0, coeffs, N, dt, reference_velocity);

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
    state.x(0) = State<T>::x(x0.data());
    state.y(0) = State<T>::y(x0.data());
    state.psi(0) = State<T>::psi(x0.data());
    state.v(0) = State<T>::v(x0.data());
    state.crosstrack_error(0) = State<T>::crosstrack_error(x0.data());
    state.psi_error(0) = State<T>::psi_error(x0.data());

    /* Finally, set the remaining states using the previously-specified control values */
    for (size_t t = 1; t < N; ++t)
        fg.dynamics(state.at(t), state.at(t - 1), control.at(t - 1));

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
    for (size_t i = 0; i < n_vars; i++)
    {
        vars_lowerbound[i] = -1.0e19;
        vars_upperbound[i] = 1.0e19;
    }

    /* Next, constrain the steering angles to +/- 25 degrees (in radians),
     * and the throttles to +/-1.
     *
     * TODO: Might want to use something smaller for the angles to make sure that the controller is not erratic. */
    T *control_lowerbounds = vars_lowerbound.data() + N * n_states;
    T *control_upperbounds = vars_upperbound.data() + N * n_states;
    for (size_t t = 0; t < N - 1; ++t)
    {

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
    for (size_t i = 0; i < n_constraints; ++i)
    {
        constraints_lowerbound[i] = 0;
        constraints_upperbound[i] = 0;
    }

    T *constraint_lowerbounds = constraints_lowerbound.data();
    State<T>::x(constraint_lowerbounds) = State<T>::x(x0.data());
    State<T>::y(constraint_lowerbounds) = State<T>::y(x0.data());
    State<T>::psi(constraint_lowerbounds) = State<T>::psi(x0.data());
    State<T>::v(constraint_lowerbounds) = State<T>::v(x0.data());
    State<T>::crosstrack_error(constraint_lowerbounds) = State<T>::crosstrack_error(x0.data());
    State<T>::psi_error(constraint_lowerbounds) = State<T>::psi_error(x0.data());

    T *constraint_upperbounds = constraints_upperbound.data();
    State<T>::x(constraint_upperbounds) = State<T>::x(x0.data());
    State<T>::y(constraint_upperbounds) = State<T>::y(x0.data());
    State<T>::psi(constraint_upperbounds) = State<T>::psi(x0.data());
    State<T>::v(constraint_upperbounds) = State<T>::v(x0.data());
    State<T>::crosstrack_error(constraint_upperbounds) = State<T>::crosstrack_error(x0.data());
    State<T>::psi_error(constraint_upperbounds) = State<T>::psi_error(x0.data());

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
    options += "Numeric max_cpu_time          0.5\n";

    // place to return solution
    CppAD::ipopt::solve_result<Dvector> solution;

    // solve the problem
    CppAD::ipopt::solve<Dvector, FG_eval>(
        options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
        constraints_upperbound, fg, solution);

    //
    // Compute the solution and show some diagnostics
    //
    bool ok = true;
    ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

    auto cost = solution.obj_value;
    std::cout << "cost = " << cost << ", status = " << solution.status << std::endl
              << std::endl;

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
    for (size_t i = 0; i < N; ++i)
    {
        x_vals[i] = State<T>::x(states);
        y_vals[i] = State<T>::y(states);
        states += n_states;
    }

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
