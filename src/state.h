#ifndef STATE_HEADER
#define STATE_HEADER

#include <stdlib.h>

/**
 * This class is a simply interface that interfaces with the underlying memory,
 * so that your code doesn't have to handle all of the nasty indexing over and over again.
 */
template <typename T, typename Index = size_t>
class State
{

    static const Index n = 6;

  public:
    T *data;

    State(T *data) : data(data) {}

    /** Return a pointer to all of the variables starting at the specified time index.*/
    T *at(Index i)
    {
        return data + i * n;
    }

    ////////////////////////////////////////////
    template <typename U = T>
    static U &x(U *data)
    {
        return data[0];
    }

    template <typename U = T>
    static U &y(U *data)
    {
        return data[1];
    }

    template <typename U = T>
    static U &psi(U *data)
    {
        return data[2];
    }

    template <typename U = T>
    static U &v(U *data)
    {
        return data[3];
    }

    template <typename U = T>
    static U &crosstrack_error(U *data)
    {
        return data[4];
    }

    template <typename U = T>
    static U &psi_error(U *data)
    {
        return data[5];
    }

    ////////////////////////////////////////////

    T &x(Index i)
    {
        return x(at(i));
    }

    T &y(Index i)
    {
        return y(at(i));
    }

    T &psi(Index i)
    {
        return psi(at(i));
    }

    T &v(Index i)
    {
        return v(at(i));
    }

    T &crosstrack_error(Index i)
    {
        return crosstrack_error(at(i));
    }

    T &psi_error(Index i)
    {
        return psi_error(at(i));
    }
};

#endif /* STATE_HEADER */
