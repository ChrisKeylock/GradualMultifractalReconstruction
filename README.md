# GradualMultifractalReconstruction
# Gradual Multifractal Reconstruction in two dimensions.
# The base algorithm (the IAAWT) was developed in:
#
# Keylock, C.J. 2017. Multifractal surrogate-data generation algorithm that 
# preserves pointwise Hölder regularity structure, with initial applications
# to turbulence, Physical Review E 95, 032123, https://doi.org/10.1103/PhysRevE.95.032123.
#
# GMR was introduced in:
# Keylock, C.J. 2018. Gradual multifractal reconstruction of time-series: 
# Formulation of the method and an application to the coupling between stock
# market indices and their Hölder exponents, Physica D 368, 1-9, https://doi.org/10.1016/j.physd.2017.11.011
#
# An example in two dimensions appeared in
# Keylock, C.J. 2019. Hypothesis testing for nonlinear phenomena in the 
# geosciences using synthetic, surrogate data, Earth and Space Science 6, doi: 10.1029/2018EA000435
# 
# And was first used properly in a geomorphological application:
# Keylock, C.J., Singh, A., Passalacqua, P., and Foufoula-Georgiou, E. 2021.
# ''Dissecting'' Landscapes with Holder Exponents to Reconcile Process and Form

# The code here is adopted to employ the new functions that exist in the
# Matlab Wavelet Toolbox. Our original code used the original toolbox
# developed by Nick Kingsbury who invented the DTCWT transform:
#
# Kingsbury, N. 2001. Complex wavelets for shift invariant analysis and
# filtering of signals, Appl. Comput. Harmon. Anal. 10, 234-253.

# The user inputs a 2-D signal, which is square with sides 2^N where N is an
# integer. The user states the number of surrogates to generate and the
# threhold (the GMR parameter, which defaults to 0 - the IAAWT algorithm). 

# The code returns a structure containing the surrogates 
# as well as a convergence error measure.
