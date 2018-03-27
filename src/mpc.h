#ifndef MPC_HEADER
#define MPC_HEADER

#include <vector>
#include "Eigen/Core"
#include "base_type.h"

class MPC {

    /* Constants used to define the optimization problem */
    const size_t N;
    const T dt;
    const T reference_velocity;

public:

    /* These values will be filled after solving */
    std::vector<T> x_vals;
    std::vector<T> y_vals;
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

    /* Virtual destructor in case we extend this class */
    virtual ~MPC() {
    }

    /**
     * Given an initial state and reference polynomial to track,
     * solve determine the optimal trajectory,
     * and return the first actuations
     */
    void solve(const std::vector<T> &initial_state,
               const std::vector<T> &reference_polynomial,
               const std::vector<T> &previous_controls);
};

#endif /* MPC_HEADER */
