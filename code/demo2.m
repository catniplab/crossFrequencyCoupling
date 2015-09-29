%% Demonstrate Phase-Amplitude coupling with a clean sample signal
% Purpose is to visualize the concept and provide intuition for the readers

rand('seed', 20150928); randn('seed', 20150928);

%% Generate the two signals that are interacting
fs = 1000;
f1 = 10;
f2 = 75;
T = 2000;
sigma1 = 0.1;
sigman = 0.15;
e1 = cos(2 * pi * f1/fs * (1:T)' + sigma1 * cumsum(randn(T, 1)));
e2 = exp(-3*e1) .* cos(2 * pi * f2/fs * (1:T)') / 10;
tr = ((1:T)-1)*1e-3;
x = e1 + e2 + sigman * randn(T, 1);

%%
[b, a, fCenterList, nTap, fEdges] = filterBank_prop(f1, f2, 1, 3, 1000, 'fir1', f2-f1);
[x_filtered, x_analytic] = applyFilterBankThenHT_filtfilt(b, a, x, nTap);
amplitude = abs(x_analytic);
phase = angle(x_analytic);

%%
tshift = nTap/fs;
tr2 = tr(1:end-2*nTap);
tr = tr - tshift;

%%
fig = figure(924); clf;
ax = [];
ax(1) = subplot(4,1,1); cla; hold on
plot(tr, x, 'k', 'LineWidth', 2);
plot(tr, e1, 'b--');
plot(tr, e2, 'r--');
xl = [0, 0.3];
xlim(xl);
title('synthetic nonlinear LFP');

%%
ax(2) = subplot(4,1,2); cla; hold on
plot(tr2, x_filtered(:, 1));
xlim(xl);
plot(tr2, amplitude(:, 1));
title(sprintf('bandpass at %d Hz', f1));

pos = get(ax(2), 'Position');
pos(4) = pos(4)/20; % 10th of the height of subplot
ax(3) = axes('Position', pos);
imagesc(tr2, 1, phase(:, 1)');
xlim(xl);
set(ax(3), 'Xtick', [], 'YTick', [], 'box', 'off');
cmap = colormap('hsv');

%%

ax(4) = subplot(4,1,3); cla; hold on
plot(tr2, x_filtered(:, 2));
xlim(xl);
plot(tr2, amplitude(:, 2));
title(sprintf('bandpass at %d Hz', f2));

pos = get(ax(4), 'Position');
pos(4) = pos(4)/20; % 10th of the height of subplot
ax(5) = axes('Position', pos);
imagesc(tr2, 1, phase(:, 2)');
xlim(xl);
set(ax(5), 'Xtick', [], 'YTick', [], 'box', 'off');
cmap = colormap('hsv');

%%
set(ax(1:5), 'YTick', []);
set(ax(2:5), 'XTickLabel', []);

%%
compound_signal = amplitude(:, 2) .* exp(1i*phase(:, 1));
ax(6) = subplot(4,2,7); cla; hold on
plot(exp(1i*linspace(0, 2*pi)), '-', 'Color', 0.7 * [1, 1, 1], 'LineWidth', 2);
for k = 1:numel(compound_signal)
    kColor = ceil((1 + angle(compound_signal(k)) / pi) * size(cmap, 1) / 2);
    plot(compound_signal(k), '.', 'Color', cmap(kColor,:))
end
plot(mean(compound_signal), 'xk', 'MarkerSize', 10)
axis equal; grid on
axis([-2,2,-1.2,1.2])
title('compound signal')
xlabel('real');
ylabel('imaginary');

%%
rpa = amplitude(randperm(size(amplitude, 1)), 2);
compound_signal_shuffled = rpa .* exp(1i*phase(:, 1));

ax(7) = subplot(4,2,8); cla; hold on
plot(exp(1i*linspace(0, 2*pi)), '-', 'Color', 0.7 * [1, 1, 1], 'LineWidth', 2);
for k = 1:numel(compound_signal_shuffled)
    kColor = ceil((1 + angle(compound_signal_shuffled(k)) / pi) * size(cmap, 1) / 2);
    plot(compound_signal_shuffled(k), '.', 'Color', cmap(kColor,:))
end
plot(mean(compound_signal_shuffled), 'xk', 'MarkerSize', 10);
axis equal; grid on
axis([-2,2,-1.2,1.2])
title('shuffled surrogate')

%%
set(fig, 'PaperUnit', 'inches', 'PaperSize', [5, 10], 'PaperPosition', [0, 0, 5, 10]);
saveas(fig, 'demo2.png');
saveas(fig, 'demo2.pdf');