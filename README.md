# Expansion_Powerlaw_Fit

Based on the algorithm discussed in the materials and methods of [this paper](http://www.pnas.org/cgi/doi/10.1073/pnas.0710150104),
I created a `cython` script that takes a time trace as input and then:

1. Creates a window over a given time interval
2. Uses linear regression (using the `cython_gsl` package) to fit the deterministic motion of data over the window.
3. Determines the average mean squared displacement from the linear regression over the window.
4. Loop over all possible windows and apply steps 1-3. 

See the [IPython Notebook](https://github.com/Range-Expansions/expansion_powerlaw_fit/blob/master/doc/powerlaw_fit_example_for_website.ipynb) 
illustrating how to  use the package in the **doc** folder.