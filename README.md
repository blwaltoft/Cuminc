# Cuminc
A SAS-macro for estimation of the cumulative incidence using Poisson regression

In survival analyses, we often estimate the hazard rate of a specific cause. Sometimes the
main focus is not the hazard rates but the cumulative incidences, i.e., the probability of
having failed from a specific cause prior to a given time. The cumulative incidences may
be calculated using the hazard rates, and the hazard rates are often estimated by the Cox
regression. This procedure may not be suitable for large studies due to limited computer
resources. Instead one uses Poisson regression, which approximates the Cox regression.
Rosthøj et al. presented a SAS-macro for the estimation of the cumulative incidences based
on the Cox regression. I present the functional form of the probabilities and variances when
using piecewise constant hazard rates and a SAS-macro for the estimation using Poisson
regression. The use of the macro is demonstrated through examples and compared to the
macro presented by Rosthøj et al.

The paper for this software can be found here:

http://www.cmpbjournal.com/article/S0169-2607(08)00210-1/fulltext
