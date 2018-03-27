#ifndef COST_HEADER
#define COST_HEADER

#include <cppad/cppad.hpp>
#include "state.h"
#include "control.h"

using CppAD::AD;

template<typename T>
class Cost {

public:

    using ADvector =
    CPPAD_TESTVECTOR(AD<T>);

    /* Constants */
    const T Lf = 2.67;
    const std::vector<T> initial_state;
    const std::vector<T> reference_polynomial;
    const size_t N;
    const T dt;
    const T reference_velocity;

    Cost(const std::vector<T> &initial_state,
         const std::vector<T> &reference_polynomial,
         size_t N,
         T dt,
         T reference_velocity)
            : reference_polynomial(reference_polynomial),
              N(N),
              dt(dt),
              reference_velocity(reference_velocity) {
    }

    /** This function contains our model of the dynamics, that is,
     *
     * State(t) = dynamics( State(t0), Control(t0) )
     *
     */
    template<typename V>
    void dynamics(V *t, V *t0, V *u_t0) {

        /* Evaluate the polynomial at the previous value of x */
        V f0 = 0.0;
        for (size_t i = 0; i < reference_polynomial.size(); i++) {
            f0 += reference_polynomial[i] * pow(State<V>::x(t0), i);
        }

        // TODO. No need to keep recalculating
        V psides0 = CppAD::atan(reference_polynomial[1]);

        /* Now set the values in the state  */
        State<V>::x(t) = State<V>::x(t0) + State<V>::v(t0) * CppAD::cos(State<V>::psi(t0)) * dt;
        State<V>::y(t) = State<V>::y(t0) + State<V>::v(t0) * CppAD::sin(State<V>::psi(t0)) * dt;
        State<V>::psi(t) = State<V>::psi(t0) + State<V>::v(t0) * Control<V>::steering_angle(u_t0) / Lf * dt;
        State<V>::v(t) = State<V>::v(t0) + Control<V>::throttle(u_t0) * dt;
        State<V>::crosstrack_error(t) = (f0 - State<V>::y(t0))
                                        + (State<V>::v(t0) * CppAD::sin(State<V>::psi_error(t0)) * dt);
        State<V>::psi_error(t) = (State<V>::psi(t0) - psides0)
                                 + (State<V>::v(t0) * Control<V>::steering_angle(u_t0) / Lf * dt);
    }

    // `fg` is a vector containing the cost and constraints.
    // `vars` is a vector containing the variable values (state & actuators).
    void operator()(ADvector &constraints, ADvector &vars) {

        State<AD<T>> constraint(constraints.data() + 1);
        State<AD<T>> state(vars.data());
        Control<AD<T>> control(vars.data() + 6 * N);

        /* First, we calculate the cost. This value should be put in the first entry of fg_vector */
        AD<T> & cost = constraints[0];

        /* The part of the cost based on the reference state. */
        for (size_t t = 0; t < N; ++t) {
            cost += CppAD::pow(state.crosstrack_error(t), 2);
            cost += CppAD::pow(state.psi_error(t), 2);
            cost += CppAD::pow(state.v(t) - reference_velocity, 2);
        }

        /* Minimize the use of actuators. */
        for (size_t t = 0; t < N - 1; ++t) {
            cost += CppAD::pow(control.steering_angle(t), 2);
            cost += CppAD::pow(control.throttle(t), 2);
        }

        /* Minimize the value gap between sequential actuations. */
        for (size_t t = 0; t < N - 2; ++t) {
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
        for (size_t t = 1; t < N; t++) {

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

#endif /* COST_HEADER */