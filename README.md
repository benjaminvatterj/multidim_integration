# Multidimensional Integration routines and Fixed point finder.

This file contains two simple extensions of numerical integrations written in C.
The first is a n-dimensional quadrature integration, which extends the GSL integration
routine "gsl_integration_qags".
The following is a simple Monte Carlo integration. The two have matching signatures to
allow for dynamically exchanging the two according the number of dimensions being used
and the precision required.

Additionally the fixed_point file has a very simple fixed point finder, that uses the Aitken's
delta method for acceleration. It is a copy of the one found in the SciPy modules.

Note: this isn't novel in any sense. I basically rewrote routines I normally use in Python,
and couldn't find in C. Given that I use these very often, maybe other's might benefit from
them too.
