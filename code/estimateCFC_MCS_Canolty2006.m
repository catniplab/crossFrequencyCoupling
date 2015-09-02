function CFC = estimateCFC_Canolty2006(fLow, fHigh, xLow, xHigh, aLow, phiLow, aHigh, phiHigh)
% Weakness: can't detect if coupling is at double the frequency
% Weakness: can't detect if amplitude has symmetric deviation (for example
%           increased variance)

CFC = abs(mean(aHigh .* exp(1i*phiLow)));