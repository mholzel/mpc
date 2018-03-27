#include "cost.h"

/**
  * This is where the cost is actually computed.
  * The first entry in "constraints" is the cost, followed by the constraints.
  */
void Cost::operator()(ADvector &constraints, ADvector &vars) {

    /* The first entry in the constraints vector is actually the cost. */
    CppAD::AD <T> &cost = constraints[0];
    cost = 0.0;

    /* Compute the cost */
    Control<CppAD::AD < T>>
    control(vars.data());
    for (size_t i = 0; i < n_states; ++i)
        x0[i] = initial_state[i];
    for (size_t t = 0; t < N; ++t) {

        /* Calculate x(t) */
        dynamics(x.data(), x0.data(), control.at(t));
        x0 = x;

        /* The part of the cost based on the reference state. */
        T scale = 100 * pow(.99, t);
        cost += CppAD::pow(States::crosstrack_error(x.data()), 2) * scale;
        cost += CppAD::pow(States::psi_error(x.data()), 2) * scale;
        cost += CppAD::pow(States::v(x.data()) - reference_velocity, 2) * scale;

        /* Minimize the use of actuators. */
        // cost += CppAD::pow(control.steering_angle(t), 2);
        // cost += CppAD::pow(control.throttle(t), 2);

        /* Minimize the value gap between sequential actuations. */
//        if (t < N - 1) {
//            cost += CppAD::pow(control.steering_angle(t + 1) - control.steering_angle(t), 2) / 10;
//            cost += CppAD::pow(control.throttle(t + 1) - control.throttle(t), 2) / 10;
//        }
    }
}

using States = State<CppAD::AD < T>>;
using Controls = Control<CppAD::AD < T>>;

void Cost::dynamics(CppAD::AD <T> *t, CppAD::AD <T> *t0, CppAD::AD <T> *u_t0) {

    /* Evaluate the polynomial at the previous value of x */
    CppAD::AD <T> f0 = 0.0;
    for (size_t i = 0; i < reference_polynomial.size(); i++)
        f0 += reference_polynomial[i] * pow(States::x(t0), i);

    // TODO. Is this correct?
    CppAD::AD <T> psides0 = CppAD::atan(reference_polynomial[1]);

    /* Now set the values in the state  */
    States::x(t) = States::x(t0) + States::v(t0) * CppAD::cos(States::psi(t0)) * dt;
    States::y(t) = States::y(t0) + States::v(t0) * CppAD::sin(States::psi(t0)) * dt;
    States::psi(t) = States::psi(t0) + States::v(t0) * Controls::steering_angle(u_t0) / Lf * dt;
    States::v(t) = States::v(t0) + Controls::throttle(u_t0) * dt;
    States::crosstrack_error(t) =
            (f0 - States::y(t0)) + (States::v(t0) * CppAD::sin(States::psi_error(t0)) * dt);
    States::psi_error(t) =
            (States::psi(t0) - psides0) + (States::v(t0) * Controls::steering_angle(u_t0) / Lf * dt);
}

void Cost::dynamics(T *t, T *t0, T *u_t0) {

    /* Evaluate the polynomial at the previous value of x */
    T f0 = 0.0;
    for (size_t i = 0; i < reference_polynomial.size(); i++)
        f0 += reference_polynomial[i] * pow(State<T>::x(t0), i);

    // TODO. Is this correct?
    T psides0 = atan(reference_polynomial[1]);

    /* Now set the values in the state  */
    State<T>::x(t) = State<T>::x(t0) + State<T>::v(t0) * cos(State<T>::psi(t0)) * dt;
    State<T>::y(t) = State<T>::y(t0) + State<T>::v(t0) * sin(State<T>::psi(t0)) * dt;
    State<T>::psi(t) = State<T>::psi(t0) + State<T>::v(t0) * Control<T>::steering_angle(u_t0) / Lf * dt;
    State<T>::v(t) = State<T>::v(t0) + Control<T>::throttle(u_t0) * dt;
    State<T>::crosstrack_error(t) =
            (f0 - State<T>::y(t0)) + (State<T>::v(t0) * sin(State<T>::psi_error(t0)) * dt);
    State<T>::psi_error(t) =
            (State<T>::psi(t0) - psides0) + (State<T>::v(t0) * Control<T>::steering_angle(u_t0) / Lf * dt);

}

std::vector<T> Cost::calculateTrajectory(T *vars) {

    Control<T> control(vars);
    std::vector<T> states((N + 1) * n_states);
    T *x0 = states.data();
    T *x = states.data() + n_states;
    for (size_t i = 0; i < n_states; ++i)
        x0[i] = initial_state[i];

    for (size_t t = 0; t < N; ++t) {
        dynamics(x, x0, control.at(t));
        x0 += n_states;
        x += n_states;
    }
    return states;
}
