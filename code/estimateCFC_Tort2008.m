function [CFC, extra] = estimateCFC_Tort2008(fLow, fHigh, xLow, xHigh, aLow, phiLow, aHigh, phiHigh, param)
% CFC = estimateCFC_Tort2008(fLow, fHigh, xLow, xHigh, aLow, phiLow, aHigh, phiHigh, param)
% Entropy of average amplitude conditioned on phase distribution method
%
% Extra parameter: # of phase bins (default 18).

if nargin > 8
    nBin = param;
    assert(numel(nBin) == 1, 'single param for number of bins')
    assert(res(nBin, 1) == 0, 'number of bins should be integer')
    assert(nBin > 1, 'number of bins should be greater than 2');
else
    nBin = 18;
end

N = size(xLow, 2);
H = @(P) -sum(P(P~=0).*log(P(P~=0)));

phaseBins = linspace(-pi, pi, nBin + 1);
phaseBinCenters = (phaseBins(1:end-1) + phaseBins(2:end))/2;

aacpd = zeros(numel(phaseBins)-1, N);
for k = 1:N
    for kBin = 1:numel(phaseBins)-1
        bIdx = phiLow(:,k) >= phaseBins(kBin) & phiLow(:,k) < phaseBins(kBin+1);
        aacpd(kBin,:) = mean(aHigh(bIdx,k));
    end
    
    % normalize to make it a distribution
    aacpd = aacpd / sum(aacpd);
    
    Hmax = -log(1/nBin);
    CFC(k) = (Hmax - H(aacpd))/Hmax;
end

if nargout > 1
    extra.nBin = nBin;
    extra.phaseBins = phaseBins;
    extra.phaseBinCenters = phaseBinCenters;
    extra.aacpd = aacpd;
end