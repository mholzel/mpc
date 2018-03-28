#ifndef COST_HEADER
#define COST_HEADER

#include "control.h"
#include "state.h"
#include "math.h"
#include "base_type.h"
#include <vector>
#include <cppad/cppad.hpp>

class Cost {
public:

    /* Types */
    using States = State<CppAD::AD < T>>;
    using Controls = Control<CppAD::AD < T>>;

    /* Constants */
    const T Lf;
    const std::vector<T> initial_state;
    const std::vector<T> reference_polynomial;
    const size_t N;
    const T dt;
    const T reference_velocity;
    const size_t n_states = 6;

    /* Declared to avoid repeated reallocation */
    std::vector<CppAD::AD < T>> x;
    std::vector<CppAD::AD < T>> x0;

    Cost(const T Lf,
         const std::vector<T> &initial_state,
         const std::vector<T> &reference_polynomial,
         size_t N,
         T dt,
         T reference_velocity)
            : Lf(Lf),
              initial_state(initial_state),
              reference_polynomial(reference_polynomial),
              N(N),
              dt(dt),
              reference_velocity(reference_velocity),
              x(n_states),
              x0(n_states) {
        assert(initial_state.size() == n_states);
    }

    using ADvector =
    CPPAD_TESTVECTOR(CppAD::AD<double>);

    /**
     * This function contains our model of the dynamics, that is,
     *
     * State(t) = dynamics( State(t0), Control(t0) )
     *
     */
    void dynamics(CppAD::AD <T> *t, CppAD::AD <T> *t0, CppAD::AD <T> *u_t0);

    void dynamics(T *t, T *t0, T *u_t0);

    void operator()(ADvector &fg, ADvector &vars);

    std::vector<T> calculateTrajectory(T *vars);
};

#endif /* COST_HEADER */
