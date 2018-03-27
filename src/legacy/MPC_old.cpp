#include "MPC_ic.h"
#include <math.h>
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen/Core"
#include "Eigen/QR"

using CppAD::AD;

// This value assumes the model presented in the classroom is used.
//
// It was obtained by measuring the radius formed by running the vehicle in the
// simulator around in a circle with a constant steering angle and velocity on a
// flat terrain.
//
// Lf was tuned until the the radius formed by the simulating the model
// presented in the classroom matched the previous radius.
//
// This is the length from front to CoG that has a similar radius.
const double Lf = 2.67;

class FG_eval {

    using Get = VariableGetter<T, N>;

public:

    Eigen::VectorXd coeffs;
    const size_t N;
    const double dt;
    const double reference_velocity;
    const size_t x_start;
    const size_t y_start;
    const size_t psi_start;
    const size_t v_start;
    const size_t cte_start;
    const size_t epsi_start;
    const size_t delta_start;
    const size_t a_start;

    FG_eval(Eigen::VectorXd coeffs,
            size_t N,
            double dt,
            double reference_velocity)
            : coeffs(coeffs),
              N(N),
              dt(dt),
              reference_velocity(reference_velocity),
              x_start(0),
              y_start(x_start + N),
              psi_start(y_start + N),
              v_start(psi_start + N),
              cte_start(v_start + N),
              epsi_start(cte_start + N),
              delta_start(epsi_start + N),
              a_start(delta_start + N - 1) {
    }


    typedef CPPAD_TESTVECTOR(AD
    <double>)
    ADvector;

    // `fg` is a vector containing the cost and constraints.
    // `vars` is a vector containing the variable values (state & actuators).
    void operator()(ADvector &fg, const ADvector &vars) {
        // NOTE: You'll probably go back and forth between this function and
        // the Solver function below.
        // The cost is stored is the first element of `fg`.
        // Any additions to the cost should be added to `fg[0]`.
        fg[0] = 0;

        // The part of the cost based on the reference state.
        for (size_t t = 0; t < N; t++) {
            fg[0] += CppAD::pow(vars[cte_start + t], 2);
            fg[0] += CppAD::pow(vars[epsi_start + t], 2);
            fg[0] += CppAD::pow(vars[v_start + t] - reference_velocity, 2);
        }

        // Minimize the use of actuators.
        for (size_t t = 0; t < N - 1; t++) {
            fg[0] += CppAD::pow(vars[delta_start + t], 2);
            fg[0] += CppAD::pow(vars[a_start + t], 2);
        }

        // Minimize the value gap between sequential actuations.
        for (size_t t = 0; t < N - 2; t++) {
            fg[0] += CppAD::pow(vars[delta_start + t + 1] - vars[delta_start + t], 2);
            fg[0] += CppAD::pow(vars[a_start + t + 1] - vars[a_start + t], 2);
        }

        //
        // Setup Constraints
        //
        // NOTE: In this section you'll setup the model constraints.

        // Initial constraints
        //
        // We add 1 to each of the starting indices due to cost being located at
        // index 0 of `fg`.
        // This bumps up the position of all the other values.
        fg[1 + x_start] = vars[x_start];
        fg[1 + y_start] = vars[y_start];
        fg[1 + psi_start] = vars[psi_start];
        fg[1 + v_start] = vars[v_start];
        fg[1 + cte_start] = vars[cte_start];
        fg[1 + epsi_start] = vars[epsi_start];

        // The rest of the constraints
        for (size_t t = 1; t < N; t++) {
            // The state at time t+1 .
            AD<double> x1 = vars[x_start + t];
            AD<double> y1 = vars[y_start + t];
            AD<double> psi1 = vars[psi_start + t];
            AD<double> v1 = vars[v_start + t];
            AD<double> cte1 = vars[cte_start + t];
            AD<double> epsi1 = vars[epsi_start + t];

            // The state at time t.
            AD<double> x0 = vars[x_start + t - 1];
            AD<double> y0 = vars[y_start + t - 1];
            AD<double> psi0 = vars[psi_start + t - 1];
            AD<double> v0 = vars[v_start + t - 1];
            AD<double> cte0 = vars[cte_start + t - 1];
            AD<double> epsi0 = vars[epsi_start + t - 1];

            // Only consider the actuation at time t.
            AD<double> delta0 = vars[delta_start + t - 1];
            AD<double> a0 = vars[a_start + t - 1];

            /* Evaluate the polynomial at this value of x */
            AD<double> f0 = 0.0;
            for (int i = 0; i < coeffs.size(); i++) {
                f0 += coeffs[i] * pow(x0, i);
            }
            AD<double> psides0 = CppAD::atan(coeffs[1]);

            // Here's `x` to get you started.
            // The idea here is to constraint this value to be 0.
            //
            // Recall the equations for the model:
            // x_[t+1] = x[t] + v[t] * cos(psi[t]) * dt
            // y_[t+1] = y[t] + v[t] * sin(psi[t]) * dt
            // psi_[t+1] = psi[t] + v[t] / Lf * delta[t] * dt
            // v_[t+1] = v[t] + a[t] * dt
            // cte[t+1] = f(x[t]) - y[t] + v[t] * sin(epsi[t]) * dt
            // epsi[t+1] = psi[t] - psides[t] + v[t] * delta[t] / Lf * dt
            fg[1 + x_start + t] = x1 - (x0 + v0 * CppAD::cos(psi0) * dt);
            fg[1 + y_start + t] = y1 - (y0 + v0 * CppAD::sin(psi0) * dt);
            fg[1 + psi_start + t] = psi1 - (psi0 + v0 * delta0 / Lf * dt);
            fg[1 + v_start + t] = v1 - (v0 + a0 * dt);
            fg[1 + cte_start + t] =
                    cte1 - ((f0 - y0) + (v0 * CppAD::sin(epsi0) * dt));
            fg[1 + epsi_start + t] =
                    epsi1 - ((psi0 - psides0) + v0 * delta0 / Lf * dt);
        }
    }
};

void MPC::Solve(Eigen::VectorXd x0, Eigen::VectorXd coeffs) {
    typedef CPPAD_TESTVECTOR(
    double) Dvector;

    // object that computes objective and constraints
    FG_eval fg(coeffs, N, dt, reference_velocity);

    double x = x0[0];
    double y = x0[1];
    double psi = x0[2];
    double v = x0[3];
    double cte = x0[4];
    double epsi = x0[5];

    // number of independent variables
    // N timesteps == N - 1 actuations
    size_t n_vars = N * 6 + (N - 1) * 2;
    // Number of constraints
    size_t n_constraints = N * 6;

    // Initial value of the independent variables.
    // Should be 0 except for the initial values.
    Dvector vars(n_vars);
    for (size_t i = 0; i < n_vars; i++) {
        vars[i] = 0.0;
    }

    // Set the initial variable values
    for (size_t t = 0; t < N; t++) {
        vars[fg.x_start + t] = x;
        vars[fg.y_start + t] = y;
        vars[fg.psi_start + t] = psi;
        vars[fg.v_start + t] = v;
        vars[fg.cte_start + t] = cte;
        vars[fg.epsi_start + t] = epsi;
    }

    /* TODO: Probably would be faster to either use these initial values for all states,
     * TODO: or, even better, use the coefficients to interpolate and find good initial guesses. */

    // Lower and upper limits for x
    Dvector vars_lowerbound(n_vars);
    Dvector vars_upperbound(n_vars);

    // Set all non-actuators upper and lowerlimits
    // to the max negative and positive values.
    for (size_t i = 0; i < fg.delta_start; i++) {
        vars_lowerbound[i] = -1.0e19;
        vars_upperbound[i] = 1.0e19;
    }

    // The upper and lower limits of delta are set to -25 and 25
    // degrees (values in radians).
    // NOTE: Feel free to change this to something else.
    for (size_t i = fg.delta_start; i < fg.a_start; i++) {
        vars_lowerbound[i] = -0.436332;
        vars_upperbound[i] = 0.436332;
    }

    // Acceleration/decceleration upper and lower limits.
    // NOTE: Feel free to change this to something else.
    for (size_t i = fg.a_start; i < n_vars; i++) {
        vars_lowerbound[i] = -1.0;
        vars_upperbound[i] = 1.0;
    }

    // Lower and upper limits for constraints
    // All of these should be 0 except the initial
    // state indices.
    Dvector constraints_lowerbound(n_constraints);
    Dvector constraints_upperbound(n_constraints);
    for (size_t i = 0; i < n_constraints; i++) {
        constraints_lowerbound[i] = 0;
        constraints_upperbound[i] = 0;
    }
    constraints_lowerbound[fg.x_start] = x;
    constraints_lowerbound[fg.y_start] = y;
    constraints_lowerbound[fg.psi_start] = psi;
    constraints_lowerbound[fg.v_start] = v;
    constraints_lowerbound[fg.cte_start] = cte;
    constraints_lowerbound[fg.epsi_start] = epsi;

    constraints_upperbound[fg.x_start] = x;
    constraints_upperbound[fg.y_start] = y;
    constraints_upperbound[fg.psi_start] = psi;
    constraints_upperbound[fg.v_start] = v;
    constraints_upperbound[fg.cte_start] = cte;
    constraints_upperbound[fg.epsi_start] = epsi;

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
//    options += "Numeric max_cpu_time          0.5\n";

    // place to return solution
    CppAD::ipopt::solve_result <Dvector> solution;

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
    std::cout << "cost = " << cost << ", status = " << solution.status << std::endl << std::endl;

    /* Now set the fields in the MPC controller so that they can be retrieved elsewhere. */
    x_vals.resize(N);
    y_vals.resize(N);
    for (size_t i = 0; i < N; ++i) {
        x_vals[i] = solution.x[fg.x_start + i];
        y_vals[i] = solution.x[fg.y_start + i];
    }
    normalized_steering_angle = solution.x[fg.delta_start] / 0.436332;
    normalized_throttle = solution.x[fg.a_start];
}
