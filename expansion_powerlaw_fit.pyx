import numpy as np
cimport numpy as np

from cython_gsl cimport *

cimport libc.math as math

cpdef double[:, :] extract_powerlaw_oskar(double[:] x, double[:] y):
    """Extract a powerlaw from a x, y trajectory using the procedure
    described in the 2007 PNAS paper."""

    #### Declaring variables for the main code ####

    cdef int traj_length = x.shape[0]
    cdef int window_size
    cdef int i

    cdef double[:] cur_x
    cdef double[:] cur_y
    cdef double[:] y_linear_fit
    cdef double theta, cos_theta
    cdef double[:] y_perp
    cdef double mean_msq

    cdef int dd

    #### Creating the output array to speed up computation ####

    # Create the output array for simplicity
    cdef double[:, :] output
    # Find the output dimensions.
    cdef int num_rows = 0
    for window_size in range(3, traj_length):
        for i in range(0, traj_length - window_size):
            num_rows += 1
    output = np.zeros((num_rows, 3), dtype=np.double)

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

            # Rotate and translate the coordinate system to the beginning of the contour
            theta = math.atan(c1)
            xo = cur_x[0]
            yo = c0 + c1*xo # Starting point of the linear fit

            cos_theta = math.cos(theta)
            sin_theta = math.sin(theta)

            x_prime = cur_x.copy()
            y_prime = cur_y.copy()
            for dd in range(window_size):
                delta_x = cur_x[dd] - xo
                delta_y = cur_y[dd] - yo

                x_prime[dd] = delta_x*cos_theta + delta_y*sin_theta
                y_prime[dd] = -delta_x*sin_theta + delta_y*cos_theta

            # Calculate average msq distance, do it in place for speed
            mean_msq = 0
            for dd in range(window_size):
                mean_msq += y_perp[dd]**2
            mean_msq /= window_size
            # Get the average msq distance

            output[row_count, 0] = i
            output[row_count, 1] = window_size
            output[row_count, 2] = mean_msq
            row_count += 1

    return output