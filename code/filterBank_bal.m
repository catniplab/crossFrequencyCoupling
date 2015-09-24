function [b, a, fCenterList, suggestedBoundaryRemoval, fEdges] = filterBank_bal(fCenterList, TBP, nCycleBegin, fs, designType)
% Design a balanced filter bank with equal temporal and bandwidth features
%
% filterBank_prop has terrible bandwidth for higher frequencies. We try to
% balance it out by taking more cycles for the higher frequencies.
% Number of cycles is linearly scaled as center frequency increases.
% The resulting filters are of same length, and hence the temporal resolution
% is kept constant, and the bandwidth is kept constant.
%
% Input
%    fCenterList: list of center frequencies in ascending order
%    TBP: [1] time-bandwidth product
%    nCycleBegin: number of cycles to include for the FIR filters (default: 3).
%            for the lowest center frequency.
%    fs: [1] sampling frequency
%    designType: 'fir1' [default], 'firls', or 'cheby2' for filter design
%
% Output
%    b, a: {nFilter x 1} filter coefficients
%    fCenterList: [nFilter x 1] center pass-band frequencies
%    suggestedBoundaryRemoval: [1] FIR boundary effect size
%    fEdges: [2 x nFilter] bandpass filter specifications
%
% [b, a, fCenterList, suggestedBoundaryRemoval, fEdges] = filterBank_bal(8:2:100, 1, 3, 1000, 'fir1')

if nargin < 6
    designType = 'fir1';
end

assert(fCenterList(1) == min(fCenterList), 'need to be ascending order');

nCycleList = nCycleBegin * (fCenterList / fCenterList(1));
fLenHandle = @(f,nCycle) ceil(fs * nCycle / f);
fBWHandle = @(f,nCycle) TBP * f / nCycle / 2;

switch lower(designType)
    case {'fir1'}
        fSpecHandle = @(f, df) [f-df/2, f+df/2];
        designBPhandle = @(f, nCycle) fir1(fLenHandle(f, nCycle), fSpecHandle(f, fBWHandle(f, nCycle)) * 2 / fs, 'bandpass');
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
        designBPhandle = @(f, nCycle) firls(fLenHandle(f), fSpecHandle(f, fBWHandle(f, nCycle)) * 2 / fs, [0 0 1 1 0 0]);
        isFIR = true;
        
    case {'cheby2'}
        fSpecHandle = @(f, df) [f-df/2, f+df/2];
        designBPhandle = @(f, nCycle) cheby2(4, 40, fSpecHandle(f, fBWHandle(f, nCycle)) * 2 / fs);
        isFIR = false;
        
    otherwise
        error('Unknown filter design [%s]', designType);
end

b = cell(numel(fCenterList),1);
a = cell(numel(fCenterList),1);
fEdges = zeros(2, numel(fCenterList));
fLenList = zeros(numel(fCenterList), 1);
for kCenter = 1:numel(fCenterList)
    fCenter = fCenterList(kCenter);
    nCycle = nCycleList(kCenter);
    fLenList(kCenter) = fLenHandle(fCenter, nCycle);
    fEdges(:, kCenter) = fCenter + fBWHandle(fCenter, nCycle) * [-1, 1]/2;
    if isFIR
        b{kCenter} = designBPhandle(fCenter, nCycle);
        a{kCenter} = 1;
    else
        [b{kCenter}, a{kCenter}] = designBPhandle(fCenter, nCycle);
    end
end

% The longest filter side-effect
% boundary effect of filtering to be removed for all signals
suggestedBoundaryRemoval = max(fLenList);
