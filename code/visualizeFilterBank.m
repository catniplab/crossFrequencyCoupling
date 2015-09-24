function fig = visualizeFilterBank(b, a, fCenterList, fs, outDir, prefix, groupSize)

nFilter = numel(b);
assert(numel(a) == nFilter);
nPoint = 1024; % frequency resolution
warning('off', 'MATLAB:MKDIR:DirectoryExists')
mkdir(outDir)

if nargin < 7 || isempty(groupSize)
    groupSize = nFilter;
end

for kCenter = 1:nFilter
    if mod(kCenter, groupSize) == 1
        fig = figure(43896 + ceil(kCenter/groupSize));
        set(fig, 'Units', 'inches');
        set(fig, 'PaperSize', [10 4], 'PaperPosition', [0 0 10 4]);
        clf; hold all
    end

    if groupSize ~= nFilter
        subplot(groupSize, 2, 2*(mod(kCenter-1, groupSize))+1);
        zplane(b{kCenter}, a{kCenter});
        stem(b{kCenter}); axis tight
        hold all; stem(a{kCenter} * max(b{kCenter}));
        if mod(kCenter, groupSize) == 0 || kCenter == nFilter
            xlabel('taps'); ylabel('filter coeff');
        end
        subplot(groupSize, 2, 2*(mod(kCenter-1, groupSize))+2);
    end
    
    fCenter = fCenterList(kCenter);
    hold on;
    line(fCenter * [1, 1], [-60, 0], 'Color', 'r');
    
    [h, w] = freqz(b{kCenter}, a{kCenter}, nPoint);
    plot(w * fs / 2 / pi, 10*log10(abs(h)));
    ylim([-60, 0]);
    xlim([0, 200]);
    
    if kCenter == nFilter % mod(kCenter, groupSize) == 0 || kCenter == nFilter
        xlabel('Frequency (Hz)');
        ylabel('gain (dB)');
        title('Filter bank');
        drawnow;
        saveas(fig, sprintf('%s/filterBank_%s_%d.pdf', outDir, prefix, ...
            ceil(kCenter/groupSize)));
    end
end
