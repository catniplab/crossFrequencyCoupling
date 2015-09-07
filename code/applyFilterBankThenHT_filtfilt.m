function [x_filtered, x_analytic] = applyFilterBankThenHT_filtfilt(b, a, x, nBoundaryRemoval)
% Use filtfilt to do zero-phase distortion, acausal bandpass filtering
% followed by Hilbert transform. If the bandpass is not narrow enough,
% Hilbert transform might suffer to extract envelop and phases correctly.
%
% [x_filtered, x_analytic] = applyFilterBankThenHT_filtfilt(b, a, x,
%                                                           nBoundaryRemoval);
%
% Input
%    b,a: {N x 1} transfer function coefficients
%    x: [T+2*nBoundaryRemoval x 1] broadband signal to be filtered
%    nBoundaryRemoval: [1] number of time bins to trim off from both ends
%                      for removing filter artifacts on the boundary
% 
% Output
%    x_filtered: [T x N; real] bandpass filtered output
%    x_analytic: [T x N; complex] result of hilbert transform on output

T = size(x, 1);
x_filtered = zeros(T - 2*nBoundaryRemoval, numel(b));
x_analytic = complex(x_filtered, 0);

for kCenter = 1:numel(b)
    xTemp = filtfilt(b{kCenter}, a{kCenter}, x);
    x_filtered(:, kCenter) = xTemp(nBoundaryRemoval+1:end-nBoundaryRemoval);
    x_analytic(:, kCenter) = hilbert(x_filtered(:, kCenter));
end