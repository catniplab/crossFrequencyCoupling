function CFC = estimateCFC_ESC_Bruns2004(fLow, fHigh, xLow, xHigh, aLow, phiLow, aHigh, phiHigh)
% Envolop-to-Signal Correlation
%
% Weakness: can't detect if coupling is at double the frequency
% Weakness: can't detect if amplitude has symmetric deviation (for example
%           increased variance)
% Weakness: can't discriminate between common amplidue gain and
%           phase-correlation

sz = size(xLow); N = prod(sz(2:end));
CFC = zeros(N, 1);
for k = 1:N
    CFC(k) = abs(corr(xLow(:,k), aHigh(:,k)));
end
if numel(sz) > 2; CFC = reshape(CFC, sz(2:end)); end
