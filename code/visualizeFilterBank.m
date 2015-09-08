function visualizeFilterBank(b, a, fCenterList, fs, prefix)
%%
nFilter = numel(b);
assert(numel(a) == nFilter);
nPoint = 1024; % frequency resolution

for kCenter = 1:nFilter
    if kCenter == 1 %mod(kCenter, 5) == 1
        fig = figure(43896 + ceil(kCenter/5)); set(fig, 'Units', 'inches');
        set(fig, 'PaperSize', [10 4], 'PaperPosition', [0 0 10 4]); clf;
        hold all
    end
    %subplot(5, 2, 2*(mod(kCenter-1, 5))+1);
    %zplane(b{kCenter}, a{kCenter});
    %stem(b{kCenter}); axis tight
    %hold all; stem(a{kCenter} * max(b{kCenter}));
    %if mod(kCenter, 5) == 0 || kCenter == nFilter
    %    xlabel('taps'); ylabel('filter coeff');
    %end

    %subplot(5, 2, 2*(mod(kCenter-1, 5))+2);
    
    fCenter = fCenterList(kCenter);
    hold on;
    line(fCenter * [1, 1], [-80, 0], 'Color', 'r');
    
    [h, w] = freqz(b{kCenter}, a{kCenter}, nPoint);
    plot(w * fs / 2 / pi, 10*log10(abs(h)));
    ylim([-60, 0]);
    xlim([0, 200]);
    
    if kCenter == nFilter % mod(kCenter, 5) == 0 || kCenter == nFilter
        xlabel('Frequency (Hz)');
        ylabel('gain (dB)');
        saveas(fig, sprintf('filterBank_%s_%d.png', prefix, ceil(kCenter/5)));
    end
end