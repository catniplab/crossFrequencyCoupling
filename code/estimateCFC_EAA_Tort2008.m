function [CFC, extra] = estimateCFC_EAA_Tort2008(fLow, fHigh, xLow, xHigh, aLow, phiLow, aHigh, phiHigh, param)
% CFC = estimateCFC_EAA_Tort2008(fLow, fHigh, xLow, xHigh, aLow, phiLow, aHigh, phiHigh, param)
% Entropy of Average Amplitude (EAA) conditioned on phase distribution method
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

sz = size(xLow); N = prod(sz(2:end));
H = @(P) -sum(P(P~=0).*log(P(P~=0))); % entropy function

% The original paper requires z-scoring the filtered signals. It changes
% the normalization of mean amplitude distribution. Subtracting the mean
% can't be emulated, but normalizing by the standard deviation can be done
% even after Hilbert transform.
aHigh = bsxfun(@times, aHigh, 1 ./ std(xHigh));

phaseBins = linspace(-pi, pi, nBin + 1);
phaseBinCenters = (phaseBins(1:end-1) + phaseBins(2:end))/2;

CFC = zeros(N, 1);
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

if numel(sz) > 2; CFC = reshape(CFC, sz(2:end)); end

if nargout > 1
    extra.nBin = nBin;
    extra.phaseBins = phaseBins;
    extra.phaseBinCenters = phaseBinCenters;
    extra.aacpd = aacpd;
end
