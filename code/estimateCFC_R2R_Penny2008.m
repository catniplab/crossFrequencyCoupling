function CFC = estimateCFC_R2R_Penny2008(fLow, fHigh, xLow, xHigh, aLow, phiLow, aHigh, phiHigh)
% r^2 of the residual from GLM regression of amplitude from sin and cos of phase

CFC = zeros(size(xLow, 2), 1);
T = size(phiLow, 1);
for k = 1:size(xLow, 2) 
    Xtemp = [cos(phiLow(:,k)), sin(phiLow(:,k)), ones(T, 1)];
    res = Xtemp * (Xtemp \ aHigh(:,k)) - aHigh(:,k);
    RSS = sum(res.^2); % residual sum of squares
    TSS = sum(aHigh(:,k).^2); % total sum of squares
    CFC(k) = 1 - RSS/TSS; % r^2 is the metric
end
