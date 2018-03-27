#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "state.h"
#include "control.h"

using CppAD::AD;

template<typename T>
class Cost {

public:
    using ADvector =
    CPPAD_TESTVECTOR(AD<T>);

    using States = State<AD < T>>;
    using Controls = Control<AD < T>>;

    /* Constants */
    const T Lf = 2.67;
    const std::vector<T> &initial_state;
    const std::vector<T> reference_polynomial;
    const size_t N;
    const T dt;
    const T reference_velocity;

    /* Declared to avoid repeated reallocation */
    std::vector<AD < T>> x;
    std::vector<AD < T>> x0;

    /* Constructor */
    Cost(const std::vector<T> &initial_state,
         const std::vector<T> &reference_polynomial,
         size_t N,
         T dt,
         T reference_velocity)
            : initial_state(initial_state),
              reference_polynomial(reference_polynomial),
              N(N),
              dt(dt),
              reference_velocity(reference_velocity),
              x(6),
              x0(6) {
    }

    /**
     * This function contains our model of the dynamics, that is,
     *
     * State(t) = dynamics( State(t0), Control(t0) )
     */
    void dynamics(AD <T> *t, AD <T> *t0, AD <T> *u_t0) {

        /* Evaluate the polynomial at the previous value of x */
        AD<T> f0 = 0.0;
        for (size_t i = 0; i < reference_polynomial.size(); i++)
            f0 += reference_polynomial[i] * pow(States::x(t0), i);

        // TODO. No need to keep recalculating
        T psides0 = CppAD::atan(reference_polynomial[1]);

        /* Now set the values in the state  */
        States::x(t) = States::x(t0) + States::v(t0) * CppAD::cos(States::psi(t0)) * dt;
        States::y(t) = States::y(t0) + States::v(t0) * CppAD::sin(States::psi(t0)) * dt;
        States::psi(t) = States::psi(t0) + States::v(t0) * Controls::steering_angle(u_t0) / Lf * dt;
        States::v(t) = States::v(t0) + Controls::throttle(u_t0) * dt;
        States::crosstrack_error(t) = (f0 - States::y(t0)) + (States::v(t0) * CppAD::sin(States::psi_error(t0)) * dt);
        States::psi_error(t) = (States::psi(t0) - psides0) + (States::v(t0) * Controls::steering_angle(u_t0) / Lf * dt);
    }

    /**
     * This is where the cost is actually computed.
     * The first entry in "constraints" is the cost, followed by the constraints.
     */
    void operator()(ADvector &constraints, ADvector &vars) {
        /* The first entry in the constraints vector is actually the cost. */
        constraints[0] = 0.0;
        AD<T> & cost = constraints[0];

        /* Compute the cost */
        Control<AD<T>> control(vars.data());
        for (size_t i = 0; i < x0.size(); ++i)
            x0[i] = initial_state[i];
        for (size_t t = 0; t < N; ++t) {
            /* Calculate x(t) */
            dynamics(x.data(), x0.data(), control.at(t));
            x0 = x;

            /* The part of the cost based on the reference state. */
            cost += CppAD::pow(States::crosstrack_error(x.data()), 2);
            cost += CppAD::pow(States::psi_error(x.data()), 2);
            cost += CppAD::pow(States::v(x.data()) - reference_velocity, 2);

            /* Minimize the use of actuators. */
            cost += CppAD::pow(control.steering_angle(t), 2);
            cost += CppAD::pow(control.throttle(t), 2);

            /* Minimize the value gap between sequential actuations. */
            if (t < N - 1) {
                cost += CppAD::pow(control.steering_angle(t + 1) - control.steering_angle(t), 2);
                cost += CppAD::pow(control.throttle(t + 1) - control.throttle(t), 2);
            }
        }
    }
};