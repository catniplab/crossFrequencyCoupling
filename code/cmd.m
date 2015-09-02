% rand('seed', 20150827); randn('seed', 20150827);

%% Setup frequencies of interest and corresponding filter banks
fs = 1000; % sampling frequency
f1 = 10 / fs;
f2 = 50 / fs;
df = 4 / fs; % pass band band-width

nTap = ceil(2/f1); %ceil(2/df); % For an FIR, nTap should be longer than one cycle of lowest freq
removeBoundaryHandle = @(z) z(nTap+1:end-nTap);

designBPhandle = @(f) firls(ceil(2/f), [0, f-df, f-df/2, f+df/2, f+df, 1/2] * 2, [0 0 1 1 0 0]);
b1 = designBPhandle(f1);
b2 = designBPhandle(f2);

%% Simulate some random signals
T = ceil(10 * fs);
TT = T + 2 * nTap;

switch(2)
case 1
e = randn(T + 2 * nTap, 1);
e1 = filtfilt(b1, 1, e);
e2 = filtfilt(b2, 1, e);
x = e1 .* sin((1:numel(e))'/T*8) + e2 .* e1;

case 2
e1 = cos(2 * pi * f1 * (1:TT) + 0.2 * cumsum(randn(1, TT)));
e2 = exp(-e1) .* cos(2 * pi * f2 * (1:TT));
x = e1 + e2;
x = x(:);
x = zscore(x) + 0.01 * randn(TT, 1);
end

figRaw = figure(4017); clf; plot(removeBoundaryHandle(x)); title('x(t)');
%% band-pass filter x
xLow = filtfilt(b1, 1, x);
xHigh = filtfilt(b2, 1, x);

%% Hilbert transform
X1 = hilbert(xLow);
X2 = hilbert(xHigh);

%% remove boundaries
x = removeBoundaryHandle(x);
xLow = removeBoundaryHandle(xLow);
xHigh = removeBoundaryHandle(xHigh);
X1 = removeBoundaryHandle(X1);
X2 = removeBoundaryHandle(X2);

%% Extract amplitude and phase from the analytical signal
aLow  = abs(X1); phiLow  = angle(X1);
aHigh = abs(X2); phiHigh = angle(X2);

%% Plot them
fig = figure(62789); clf;
subplot(2,2,1); plot(xLow);
hold all; plot(aLow);
subplot(2,2,3); plot(phiLow);
subplot(2,2,2); plot(xHigh);
hold all; plot(aHigh)
subplot(2,2,4); plot(phiHigh);
linkaxes(get(fig, 'Children'), 'x');

% Parameters for the surrogate generation
nSurrogate = 999; % number of shuffled surrogates requested
minTimeShift = 1 * fs; % minimum time to be separated to have phase decoherence
rpIdxAll = generateSurrogateIndices(T, minTimeShift, nSurrogate);

%% Inspect decoherence
% the phase autocorrelation should decay significantly by minTimeShift
[xcLow, lags]  = xcorr(phiLow,  ceil(min(T/4, minTimeShift*8)), 'coeff');
[xcHigh, lags] = xcorr(phiHigh, ceil(min(T/4, minTimeShift*8)), 'coeff');
figure(4125); clf; hold all; title('phase autocorrelation');
plot(lags(lags>=0), xcLow (lags>=0));
plot(lags(lags>=0), xcHigh(lags>=0));
line(minTimeShift * [1, 1], [-1, 1], 'Color', 'k');

%% Generate surrogates
xx1 = repmat(xLow(:), 1, nSurrogate+1);
xx2 = xHigh(rpIdxAll);
xaLow = repmat(aLow(:), 1, nSurrogate+1);
xphiLow = repmat(phiLow(:), 1, nSurrogate+1);
xaHigh = aHigh(rpIdxAll);
xphiHigh = phiHigh(rpIdxAll);

%%
CFC = cell(5, 1);
CFC{1} = estimateCFC_Canolty2006(f1, f2, xx1, xx2, xaLow, xphiLow, xaHigh, xphiHigh);
CFC{2} = estimateCFC_PLV_Lachaux1999(f1, f2, xx1, xx2, xaLow, xphiLow, xaHigh, xphiHigh);
CFC{3} = estimateCFC_ESC_Bruns2004(f1, f2, xx1, xx2, xaLow, xphiLow, xaHigh, xphiHigh);
CFC{4} = estimateCFC_GLM_Penny2008(f1, f2, xx1, xx2, xaLow, xphiLow, xaHigh, xphiHigh);
CFC{5} = estimateCFC_Tort2008(f1, f2, xx1, xx2, xaLow, xphiLow, xaHigh, xphiHigh);
CFCstr = {'MCS', 'PLV', 'ESC', 'R2R', 'EAA'};
% mean compound signal
% phase locking value
% envelop-to-signal correlation
% r^2 of regression
% entropy of average amplitude

%%
for kCFC = 1:numel(CFC)
    threshold = quantile(CFC{kCFC}(2:end), 1 - 0.05);
    pValue = (sum(CFC{kCFC}(2:end) > CFC{kCFC}(1)) + 1)/nSurrogate;
    m = mean(CFC{kCFC}(2:end));
    s = std(CFC{kCFC}(2:end));
    dev = (CFC{kCFC}(1) - m)/s;
    fprintf('CFC [%s]: %g\tp-value [%.4f]\tdeviation [%f]\n', CFCstr{kCFC}, CFC{kCFC}(1), pValue, dev);
end

figure(1408); clf; hold all;
hist(CFC{kCFC}(2:end), 20)
line(CFC{kCFC}(1) * [1, 1], [0, nSurrogate/10], 'Color', 'r');
line(threshold * [1, 1], [0, nSurrogate/10], 'Color', 'b');

return

%% Estimate simple PACs
rNESC = corr(cos(angle(X1)), abs(X2))
rAEC = corr(abs(X1), abs(X2))

%% Entropy of average amplitude conditioned on phase distribution method
% 
amp2 = abs(X2);
nBin = 18;
phaseBins = linspace(-pi, pi, nBin + 1);
phaseBinCenters = (phaseBins(1:end-1) + phaseBins(2:end))/2;
aacpd = zeros(numel(phaseBins)-1, 1);
for kBin = 1:numel(phaseBins)-1
    bIdx = phi >= phaseBins(kBin) & phi < phaseBins(kBin+1);
    aacpd(kBin) = mean(amp2(bIdx));
end
figure(957); clf;
bar(phaseBinCenters, aacpd);
aacpd = aacpd / sum(aacpd);
H = @(P) -sum(P(P~=0).*log(P(P~=0)));
Hmax = -log(1/nBin);
(Hmax - H(aacpd))/Hmax

%% dependence between the filtered signal and the amplitude signal
% I don't like phase signal because of phase slips
% also, don't need to be concerned with phase being circular variable
% the main problem is that overall amplitude fluctuations creates
% spurious PAC estimates.

% dep(xLow, abs(X2))
