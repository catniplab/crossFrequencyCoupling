function CFC = estimateCFC_PLV_Lachaux1999(fLow, fHigh, xLow, xHigh, aLow, phiLow, aHigh, phiHigh)
% Phase-locking value (PLV)
%
% Weakness: can't detect if coupling is at integer multiples of the lower frequency
% Weakness: when amplitude is small, and phase slip happens, you know.

CFC = abs(mean(exp(1i*(phiLow - phiHigh))));