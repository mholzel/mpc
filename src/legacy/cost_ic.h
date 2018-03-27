#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "state.h"
#include "control.h"

using CppAD::AD;

template <typename T>
class Cost
{

  public:
    using ADvector =
        CPPAD_TESTVECTOR(AD<T>);

    using States = State<AD<T>>;
    using Controls = Control<AD<T>>;

    /* Constants */
    const T Lf = 2.67;
    const std::vector<T> reference_polynomial;
    const size_t N;
    const T dt;
    const T reference_velocity;
    const size_t n_states = 6;

    /* Declared to avoid repeated reallocation */
    std::vector<AD<T>> x;
    std::vector<AD<T>> x0;

    /* Constructor */
    Cost(const std::vector<T> &reference_polynomial,
         size_t N,
         T dt,
         T reference_velocity)
        : reference_polynomial(reference_polynomial),
          N(N),
          dt(dt),
          reference_velocity(reference_velocity),
          x(n_states),
          x0(n_states)
    {
    }

    /**
     * This function contains our model of the dynamics, that is,
     *
     * State(t) = dynamics( State(t0), Control(t0) )
     */
    void dynamics(AD<T> *t, AD<T> *t0, AD<T> *u_t0)
    {

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
     * This is where the cost and constraints are actually computed.
     * The first entry in "cost_constraints" is the cost. The rest holds the constraints.
     */
    void operator()(ADvector &cost_constraints, ADvector &vars)
    {
        /* The first entry in fg is the cost.  
        The entries after that are the constraints. */
        AD<T> &cost = cost_constraints[0];
        AD<T> *constraints = cost_constraints.data() + 1;

        /* The first n_states variables hold the initial state. */
        for (size_t i = 0; i < n_states; ++i)
            x0[i] = vars[i];

        /* The remaining terms are the controls */
        Control<AD<T>> control(vars.data() + n_states);

        /* Compute the cost */
        cost = 0.0;
        for (size_t t = 0; t < N; ++t)
        {
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
            if (t < N - 1)
            {
                cost += CppAD::pow(control.steering_angle(t + 1) - control.steering_angle(t), 2);
                cost += CppAD::pow(control.throttle(t + 1) - control.throttle(t), 2);
            }
        }

        /* The rest of the constraints are just there to ensure that the initial state is not optimized */
        for (size_t i = 0; i < n_states; ++i)
            constraints[i] = vars[i];
    }
};