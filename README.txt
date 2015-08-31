Estimation of phase-amplitude-coupling / cross-frequency coupling / nested oscillation

We implement a collection of methods for quantifying cross-frequency coupling.

## Significance via surrogate shuffling

If multiple trials are available, shuffling across trial to generate surrogate is ideal. This destroyes any joint structure between the two signals that are NOT due to fixed trial structure.

The main problem of surrogates would be insufficient decoherence of the phases. A perfect oscillation, results in infinite range correlations, and thus shuffling is very tricky. Make sure to check the autocorrelation function of the phase time series.
