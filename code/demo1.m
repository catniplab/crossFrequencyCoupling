%% Basic parameters for the demo
fs = 1000; % sampling frequency
plotDelay = 0.01; % delay in sequence designed filter plots (in seconds)
outDir = 'demo_output';

warning('off', 'MATLAB:MKDIR:DirectoryExists')
mkdir(outDir);

%% Design filter bank
switch 7
    case 1
        [b, a, fCenterList, nTap] = filterBank_uniform(8, 100, 5, 5, fs, 'cheby2', 3);
        FB_prefix = 'cb2_bw5_40';
    case 2
        [b, a, fCenterList, nTap] = filterBank_uniform(8, 100, 5, 4, fs, 'firls', 3);
        FB_prefix = 'firls_bw4';
    case 3
        [b, a, fCenterList, nTap] = filterBank_uniform(8, 100, 5, 4, fs, 'fir1', 3);
        FB_prefix = 'fir1_bw4';
    case 4
        [b, a, fCenterList, nTap] = filterBank_uniform(8, 100, 5, 4, fs, 'firls', 6);
        FB_prefix = 'firls_bw4_nC6';
    case 5
        [b, a, fCenterList, nTap, fEdges] = filterBank_prop(5, 100, 3, 5, 1000, 'firls');
        FB_prefix = 'pLS_TBP3';
    case 6
        [b, a, fCenterList, nTap, fEdges] = filterBank_prop(4, 100, 2, 7, 1000, 'fir1');
        FB_prefix = 'pfir1_TBP2';
    case 7
        [b, a, fCenterList, nTap, fEdges] = filterBank_prop(4, 150, 1, 3, 1000, 'fir1', 3);
        % uniformly spaced frequency bands
        FB_prefix = 'pfir1_TBP1';
    case 8
        [b, a, fCenterList, nTap, fEdges] = filterBank_prop(4, 150, 1, 3, 1000, 'cheby2', 3);
        FB_prefix = 'pchb2_TBP1';
    case 9
        [b, a, fCenterList, nTap, fEdges] = filterBank_bal(4:3:150, 1, 3, 1000, 'fir1');
        FB_prefix = 'bfir1_1';
end
figVis = visualizeFilterBank(b, a, fCenterList, fs, outDir, FB_prefix);

%% Make some fake signal with cross-frequency coupling
TT = 5000;
T = TT + 2 * nTap;

f1 = 10 / fs;
f2 = 75 / fs;
e1 = cos(2 * pi * f1 * (1:T) + 0.1 * cumsum(randn(1, T)));
e2 = exp(-2*e1) .* cos(2 * pi * f2 * (1:T)) / 4;
x = e1 + e2;
x = x(:);

%% Add spectrally matching noise
% We want to make the Fourier transform magnitude 1
fx = fft(x);
ap = abs(fx);
aq = quantile(ap, 0.95);
ap(ap > aq) = aq;
ap = (aq - ap) / aq;

% introduce random phase, unit power noise
fnoise = fft(randn(size(x)));
fnoise = fnoise .* ap;
noise = 2 * ifft(fnoise);
x = x + noise;

clear fx ap aq fnoise

%%
figure(924); clf;
subplot(4,1,1); hold all;
pwelch(x, [], [], [], fs);
line(fs * f1 * [1, 1], [-30, -20], 'Color', 'r');
line(fs * f2 * [1, 1], [-30, -20], 'Color', 'r');
xlim([0 200]);

subplot(4,1,2); hold all
tRange = nTap + (100:500);
plot(e1, 'LineWidth', 2);
plot(e2, 'LineWidth', 2);
plot(x, 'k');
xlim([tRange(1), tRange(end)]);
lh = legend(sprintf('slow oscillation (%d Hz)', f1), sprintf('fast oscillation (%d Hz)', f2), 'noisy observation');

tRange = tRange - nTap;

%%
[x_filtered, x_analytic] = applyFilterBankThenHT_filtfilt(b, a, x, nTap);
amplitude = abs(x_analytic);
phase = angle(x_analytic);

%% Visualize some of the filtered signals
fig = figure(924);
nBand = numel(b); ph = [];
for kBand = 1:nBand
    subplot(4,1,1);
    [h, w] = freqz(b{kBand}, a{kBand}, 1024);
    if ~isempty(ph)
        delete(ph);
    end
    ph = plot(w * fs / 2 / pi, 10*log10(abs(h)));
    ylim([-60, 0]);
    
    subplot(4,1,3); cla; hold all;
    plot(tRange, x_filtered(tRange, kBand));
    plot(tRange, amplitude(tRange, kBand));
    xlim([tRange(1), tRange(end)]);
    subplot(4,1,4); cla; hold all;
    plot(tRange, phase(tRange, kBand));
    title(num2str(fEdges(:, kBand)'))
    drawnow
    pause(plotDelay);
end
% r = input('Go? ', 's'); if lower(r(1)) ~= 'y'; disp('Aborting'); return; end

ts = datestr(now,30);

set(fig, 'PaperSize', [5 8], 'PaperPosition', [0 0 5 8]);
saveas(fig, sprintf('%s/%s_%s_sample.pdf', outDir, ts, FB_prefix));

set(figVis, 'PaperSize', [8 5], 'PaperPosition', [0 0 8 5]);
saveas(figVis, sprintf('%s/%s_%s_filters.pdf', outDir, ts, FB_prefix));

%%
fLowRange = fCenterList(fCenterList <= 35); nLow = numel(fLowRange);
fHighRange = fCenterList(fCenterList >= 40); nHigh = numel(fHighRange);
highIdxOffset = find(fCenterList == fHighRange(1)) - 1;
assert(nLow > 0);
assert(nHigh > 0);

%% Get some estimators
estimators = CFCestimatorFactory('all');
nEstimator = numel(estimators);

% Parameters for the surrogate generation
nSurrogate = 999; % number of shuffled surrogates requested
minTimeShift = 0.1 * fs; % minimum time to be separated to have phase decoherence
rpIdxAll = generateSurrogateIndices(T - 2*nTap, minTimeShift, nSurrogate);

CFC = {}; clear CFC
CFC(nEstimator, nLow, nHigh) = surrogateStats(); % initialize results structure
for kLow = 1:nLow
    fLow = fLowRange(kLow);
    xLow = x_filtered(:, kLow);
    aLow = amplitude(:, kLow);
    pLow = phase(:, kLow);
    for kHigh = 1:nHigh
        fHigh = fHighRange(kHigh);
        xHigh = x_filtered(:, kHigh + highIdxOffset);
        aHigh = amplitude(:, kHigh + highIdxOffset);
        pHigh = phase(:, kHigh + highIdxOffset);

        %% Generate surrogates
        xxLow  = repmat(xLow(:), 1, nSurrogate+1);
        xxHigh = xHigh(rpIdxAll);
        xaLow  = repmat(aLow(:), 1, nSurrogate+1);
        xpLow  = repmat(pLow(:), 1, nSurrogate+1);
        xaHigh = aHigh(rpIdxAll);
        xpHigh = pHigh(rpIdxAll);

        for kEstim = 1:nEstimator
            CFCtemp = estimators(kEstim).handle(fLow, fHigh, ...
                            xxLow, xxHigh, xaLow, xpLow, xaHigh, xpHigh);

            stat = surrogateStats(CFCtemp(1), CFCtemp(2:end));
            CFC(kEstim, kLow, kHigh) = stat;
        end
    end
end

%%
for kEstim = 1:nEstimator
fig = figure(5877+kEstim); clf;

subplot(1,3,1);
imagesc(fLowRange, fHighRange, reshape([CFC(kEstim,:,:).value], nLow, nHigh)'); axis xy; colorbar; colormap('jet')
title(estimators(kEstim).desc);

subplot(1,3,2);
imagesc(fLowRange, fHighRange, reshape([CFC(kEstim,:,:).deviation], nLow, nHigh)'); axis xy; colorbar; colormap('jet')
title('deviation');
xlabel('Freq (Hz)'); ylabel('Freq (Hz)');
hold on
plot(f1*fs, f2*fs, 'ro');

subplot(1,3,3);
pValueMap = reshape([CFC(kEstim,:,:).pValue], nLow, nHigh)'; % returns rounded up p-value...
significanceMap = (pValueMap <= 0.1) + (pValueMap <= 0.05) + (pValueMap <= 0.001);
imagesc(fLowRange, fHighRange, significanceMap); axis xy;
cbh = colorbar; colormap('jet'); caxis([0 3]);
set(cbh, 'Ticks', [0 1 2 3], 'TickLabels', {'n.s.', 'p < 0.1', 'p < 0.05', 'p < 0.001'});
title('significance (bootstrap, unadjusted)')
xlabel('Freq (Hz)'); ylabel('Freq (Hz)');
hold on
plot(f1*fs, f2*fs, 'ro');
drawnow

set(fig, 'PaperSize', [8 3], 'PaperPosition', [0 0 8 3]);
saveas(fig, sprintf('%s/%s_%s_%s.pdf', outDir, ts, FB_prefix, estimators(kEstim).ID));
end
