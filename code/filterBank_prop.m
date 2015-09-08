function [b, a, fCenterList, suggestedBoundaryRemoval, fEdges] = filterBank_prop(fMin, fMax, TBP, nCycle, fs, designType, fStep)
%
% Goal of this design is to produce filter bank whose output captures
% modulation of amplitude in narrow-band oscillations. Due to uncertainty
% principle, we can't really have both narrow band and short-time scale.
% This design focuses on short-time scale, and assigns larger bandwidth
% for higher frequencies, namely, bandwidth = (center freq) * c.
%
% For an FIR, we end up with (similarly to wavelet transform):
%
%    Frequency bin: (center freq) +/- (TBP/nCycle) * (center freq)
%    Time bin: nCycle / (center freq)
% 
% Input
%    fMin: [1] minimum center frequency
%    fMax: [1] maximum center frequency
%    TBP: [1] time-bandwidth product
%    nCycle: number of cycles to include for the FIR filters (default: 3).
%            length of FIR filter is nCycle * (sampling freq)/(center freq)
%    fs: [1] sampling frequency
%    designType: 'fir1' [default], 'firls', or 'cheby2' for filter design
%    fStep: [1] step-size between center frequencies (uniform grid)
%               if fStep == 0, proportional grid is automatically generated
%               (default)
%
% Output
%    b, a: {nFilter x 1} filter coefficients
%    fCenterList: [nFilter x 1] center pass-band frequencies
%    suggestedBoundaryRemoval: [1] FIR boundary effect size
%    fEdges: [2 x nFilter] bandpass filter specifications
%
% [b, a, fCenterList, suggestedBoundaryRemoval, fEdges] = filterBank_prop(8, 100, 1, 3, 1000, 'fir1')

if nargin < 6
    designType = 'fir1';
end

if nargin < 7
    fStep = 0;
end

% If we want the bandwidth to touch each other, then
% f1 + f1 * TBP/nCycle/2 = f2 - f2 * TBP/nCycle/2
% therefore, f2 = f1 * freqIncRatio
freqIncRatio = (1 + TBP/nCycle/2) / (1 - TBP/nCycle/2);

if fStep == 0
    % freqIncRatio^n = fMax/fMin
    nBands = ceil(log(fMax/fMin) / log(freqIncRatio));
    fCenterList = fMin * (freqIncRatio.^(0:nBands));
else
    fCenterList = fMin:fStep:fMax;
end

fLenHandle = @(f) ceil(fs * nCycle / f);
fBWHandle = @(f) TBP * f / nCycle / 2;

switch lower(designType)
    case {'fir1'}
        fSpecHandle = @(f, df) [f-df/2, f+df/2];
        designBPhandle = @(f) fir1(fLenHandle(f), fSpecHandle(f, fBWHandle(f)) * 2 / fs, 'bandpass');
        isFIR = true;
        
    case {'firls'}
        % Linear phase with least squares over side lobes of the shape:
        %
        %              ------------------
        %
        %
        %    -------??                    ??----------
        %
        %    0  f-df    f-df/2  f   f+df/2   f+df fs/2
        fSpecHandle = @(f, df) [0, f-df, f-df/2, f+df/2, f+df, fs/2];
        designBPhandle = @(f) firls(fLenHandle(f), fSpecHandle(f, fBWHandle(f)) * 2 / fs, [0 0 1 1 0 0]);
        isFIR = true;
        
    case {'cheby2'}
        fSpecHandle = @(f, df) [f-df/2, f+df/2];
        designBPhandle = @(f) cheby2(4, 40, fSpecHandle(f, fBWHandle(f)) * 2 / fs);
        isFIR = false;
        
    otherwise
        error('Unknown filter design [%s]', designType);
end

b = cell(numel(fCenterList),1);
a = cell(numel(fCenterList),1);
fEdges = zeros(2, numel(fCenterList));
for kCenter = 1:numel(fCenterList)
    fCenter = fCenterList(kCenter);
    fEdges(:, kCenter) = fCenter + fBWHandle(fCenter) * [-1, 1]/2;
    if isFIR
        b{kCenter} = designBPhandle(fCenter);
        a{kCenter} = 1;
    else
        [b{kCenter}, a{kCenter}] = designBPhandle(fCenter);
    end
end

% The longest filter side-effect
% boundary effect of filtering to be removed for all signals
suggestedBoundaryRemoval = fLenHandle(fMin);