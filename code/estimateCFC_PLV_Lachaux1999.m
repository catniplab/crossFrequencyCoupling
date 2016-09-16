function CFC = estimateCFC_PLV_Lachaux1999(fLow, fHigh, xLow, xHigh, aLow, phiLow, aHigh, phiHigh)
% Phase-locking value (PLV) 
%    phase of the amplitude of the high frequency should be locked 
%    with the phase of the low frequency signal.
%
% Weakness: can't detect if coupling is at integer multiples of the lower frequency
% Weakness: when amplitude is small, and phase slip happens, you know.

% We need a second Hilbert transform to do it.
phiAmpHigh = angle(hilbert(aHigh));

sz = size(xLow);
CFC = abs(mean(exp(1i*(phiLow - phiAmpHigh))));
CFC = permute(CFC, [2:numel(sz), 1]);
