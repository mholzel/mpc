#ifndef CONTROL_HEADER
#define CONTROL_HEADER

#include <stdlib.h>

/**
 * This class is a simply interface that interfaces with the underlying memory,
 * so that your code doesn't have to handle all of the nasty indexing over and over again.
 */
template<typename T, typename Index = size_t>
class Control {

    static const Index n = 2;

public:

    T *data;

    Control(T *data) : data(data) {}

    /** Return a pointer to all of the variables starting at the specified time index. */
    T *at(Index i) {
        return data + i * n;
    }

    ////////////////////////////////////////////
    template<typename U = T>
    static U &steering_angle(U *data) {
        return data[0];
    }

    template<typename U = T>
    static const U &steering_angle(const U *data) {
        return data[0];
    }

    template<typename U = T>
    static U &throttle(U *data) {
        return data[1];
    }

    template <typename U = T>
    static const U &throttle(const U *data)
    {
        return data[1];
    }

    ////////////////////////////////////////////

    T &steering_angle(Index i) {
        return steering_angle(at(i));
    }

    T &throttle(Index i) {
        return throttle(at(i));
    }
};

#endif /* CONTROL_HEADER */
