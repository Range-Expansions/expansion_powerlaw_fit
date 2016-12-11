import numpy as np
cimport numpy as np

from cython_gsl cimport *

cimport libc.math as math

cpdef extract_powerlaw_oskar(double[:] x, double[:] y):
    """Extract a powerlaw from a x, y trajectory using the procedure
    described in the 2007 PNAS paper."""

    #### Declaring variables for the main code ####

    cdef int traj_length = x.shape[0]
    cdef int window_size
    cdef int i

    cdef double[:] cur_x
    cdef double[:] cur_y
    cdef double theta, cos_theta, sin_theta
    cdef double xo, yo, delta_x, delta_y
    cdef double[:] x_prime, y_prime
    cdef double spacing
    cdef double mean_msq
    cdef double L

    cdef int dd

    #### Creating the output array to speed up computation ####

    # Create the output array for simplicity
    cdef double[:, :] output
    # Find the output dimensions.
    cdef int num_rows = 0
    for window_size in range(3, traj_length):
        for i in range(0, traj_length - window_size):
            num_rows += 1
    output = np.zeros((num_rows, 4), dtype=np.double)

    # Loop and fit the slopes

    # Create variables for the linear fit
    cdef double c0, c1, cov00, cov01, cov11, sumsq

    cdef int row_count = 0

    for window_size in range(3, traj_length):
        # Loop over every possible window of size "window_size"
        for i in range(0, traj_length - window_size):
            cur_x = x[i:(i+window_size)]
            cur_y = y[i:(i+window_size)]

            # Linear fit cur_x and cur_y. Get a pointer to the first element
            gsl_fit_linear(&cur_x[0], 1, &cur_y[0], 1, window_size,
                          &c0, &c1, &cov00, &cov01, &cov11, &sumsq)

            # Rotate the coordinate system clockwise
            theta = -math.atan2(c1, 1) # You rise c1 over a distance of 1

            # Rotate clockwise
            cos_theta = math.cos(theta)
            sin_theta = math.sin(theta)

            x_prime = cur_x.copy()
            y_prime = cur_y.copy()
            for dd in range(window_size):
                x_temp = cur_x[dd]
                y_temp = cur_y[dd]

                x_prime[dd] = x_temp*cos_theta - y_temp*sin_theta
                y_prime[dd] = x_temp*sin_theta + y_temp*cos_theta

            # Calculate average msq distance by doing the appropriate integral
            mean_msq = 0
            for dd in range(window_size - 1):
                spacing = x_prime[dd+1] - x_prime[dd]
                mean_msq += 0.5*spacing*(y_prime[dd+1]**2 + y_prime[dd]**2)

            # x' doesn't go between zero and one actually
            L = x_prime[window_size - 1] - x_prime[0]
            if L < 0:
                print 'L is less than zero...oh dear'
                print L
                return cur_x, cur_y, x_prime, y_prime, i, window_size, c1, c0
            mean_msq /= L
            # Get the average msq distance

            output[row_count, 0] = i
            output[row_count, 1] = window_size
            output[row_count, 2] = mean_msq
            output[row_count, 3] = L
            row_count += 1

    return output