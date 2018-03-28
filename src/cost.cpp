#include "cost.h"

/**
  * This is where the cost is actually computed.
  * The first entry in "constraints" is the cost, followed by the constraints.
  */
void Cost::operator()(ADvector &constraints, ADvector &vars) {

    /* The first entry in the constraints vector is actually the cost. */
    CppAD::AD <T> &cost = constraints[0];
    cost = 0.0;

    /* Copy the initial state to x0. Note that we need to do it like this because they have different data types. */
    for (size_t i = 0; i < n_states; ++i)
        x0[i] = initial_state[i];

    /* Compute the cost */
    Control<CppAD::AD < T>>
    control(vars.data());
    for (size_t t = 0; t < N; ++t) {

        /* Calculate x(t) */
        dynamics(x.data(), x0.data(), control.at(t));
        x0 = x;

        /* Calculate the crosstrack error, which is the difference between y
         * and the reference polynomial evaluated at x */
        CppAD::AD <T> crosstrack_error = -States::y(x.data());
        for (size_t i = 0; i < reference_polynomial.size(); i++)
            crosstrack_error += reference_polynomial[i] * pow(States::x(x.data()), i);

        /* Obviously, we want the car to stay on the reference polynomial */
        cost += CppAD::pow(crosstrack_error, 2);

        /* And try to hit the reference velocity.
         * We are weighting the percent error because when the velocity is large,
         * this might otherwise be some ridiculously large number.
         * Hitting the velocity target is our least concern. Our main priority
         * is making sure that the car safely goes around the track */
        cost += CppAD::pow((States::v(x.data()) - reference_velocity), 2);

        /* Minimize the value gap between sequential actuations.
         * This depends on velocity because when we are going very fast
         * then we want to make sure that we have a very smooth trajectory.
         * This will ensure that consecutive control commands are similar,
         * and should therefore help us avoid instability.
         */
        if (t < N - 1)
            cost += CppAD::pow(control.steering_angle(t + 1) - control.steering_angle(t), 2) *
                    CppAD::pow(States::v(x.data()), 2);
    }
}

using States = State<CppAD::AD < T>>;
using Controls = Control<CppAD::AD < T>>;

void Cost::dynamics(CppAD::AD <T> *t, CppAD::AD <T> *t0, CppAD::AD <T> *u_t0) {
    States::x(t) = States::x(t0) + States::v(t0) * CppAD::cos(States::psi(t0)) * dt;
    States::y(t) = States::y(t0) + States::v(t0) * CppAD::sin(States::psi(t0)) * dt;
    States::psi(t) = States::psi(t0) + States::v(t0) * Controls::steering_angle(u_t0) / Lf * dt;
    States::v(t) = States::v(t0) + Controls::throttle(u_t0) * dt;
}

void Cost::dynamics(T *t, T *t0, T *u_t0) {
    State<T>::x(t) = State<T>::x(t0) + State<T>::v(t0) * cos(State<T>::psi(t0)) * dt;
    State<T>::y(t) = State<T>::y(t0) + State<T>::v(t0) * sin(State<T>::psi(t0)) * dt;
    State<T>::psi(t) = State<T>::psi(t0) + State<T>::v(t0) * Control<T>::steering_angle(u_t0) / Lf * dt;
    State<T>::v(t) = State<T>::v(t0) + Control<T>::throttle(u_t0) * dt;
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
