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
    const T Lf = 2.67;

    /* A flag indicating that we have solved the MPC problem at least once */
    bool initialized = false;

    /* To solve for the time delay, we will cache some values */
    T tp; // throttle_previous
    T sp; // steering_previous
    T vp; // v_previous
    T pp; // psi_previous
    T previous_global_psi;

public:

    /* These values will be filled after solving */
    std::vector<T> x_vals;
    std::vector<T> y_vals;
    T normalized_steering_angle;
    T normalized_throttle;

    /* This is our estimate of the time delay */
    T time_delay;

    /* The time delay is simply an exponential moving average of the mean.
     * Gamma is a constant in the range [0,1].
     * A value of 0 means, "just use the instantaneous time delay estimate".
     * A values of 1 means, "weight all samples equally" */
    const T gamma = .95;
    const T current_td_scale = 1 / (1 + gamma);
    const T past_td_scale = gamma / (1 + gamma);

    /* Constructor */
    MPC(size_t N,
        T dt,
        T reference_velocity,
        T initial_time_delay)
            : N(N),
              dt(dt),
              reference_velocity(reference_velocity),
              time_delay(initial_time_delay) {
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
               const std::vector<T> &previous_controls,
               const T global_psi);
};

#endif /* MPC_HEADER */
