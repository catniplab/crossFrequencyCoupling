function [CFC, extra] = estimateCFC_template(fLow, fHigh, xLow, xHigh, aLow, phiLow, aHigh, phiHigh, varargin)
% Computes a cross-frequency-couping (CFC) index from 2 band-passed
% signals, and their corresponding Hilbert transforms (or equivalent).
% Same operation is applied to each column.
%
% Input
%    fLow, fHigh: [1] center frequency for low and high band-pass filters
%                 The value is [0, 1] where 1 denotes the sampling
%                 frequency (twice the Nyquist freqeuncy)
%                 condition: fLow < fHigh
%    xLow, xHigh: [T x N] independent length-T band-passed time series
%    aLow, aHigh: [T x N] amplitude signals
%    phiLow, phiHigh: [T x N] phase signals
%
% Output
%    CFC: [N x 1] single index which is higher when CFC is (supposedly) stronger
%    extra: anything else that needs to be returned

assert(fLow < fHigh);
T = size(phiLow, 1);
N = size(phiLow, 2);
assert(all(size(xLow) == [T,N]));
assert(all(size(xHigh) == [T,N]));
assert(all(size(aLow) == [T,N]));
assert(all(size(aHigh) == [T,N]));
assert(all(size(phiHigh) == [T,N]));
CFC = nan(N, 1);
extra = [];

error('I''m just a template');

end
