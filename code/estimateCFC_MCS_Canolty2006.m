function CFC = estimateCFC_MCS_Canolty2006(fLow, fHigh, xLow, xHigh, aLow, phiLow, aHigh, phiHigh)
% Mean of the Compound Signal (MCS: high-freq amplitude and low-freq phase)
%
% Weakness: can't detect if coupling is at double the frequency
% Weakness: can't detect if amplitude has symmetric deviation (for example
%           increased variance)

CFC = abs(mean(aHigh .* exp(1i*phiLow)));
