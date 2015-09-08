function [b, a, fCenterList, suggestedBoundaryRemoval] = filterBank_uniform(fMin, fMax, fStep, df, fs, designType, nCycle)
% Generate a filter bank comprised of bandpass digital filters.
% This creates equal-bandwidth filters. If FIR with time resolution
% inversely propotional to frequency is used (for Hilber transform),
% uncertainty principle tells us, we need to sacrifice something...
% If IIR is used, then for higher frequencies we will have longer
% time scale which is not appropriate for envelop & phase extraction.
%
% Input
%    fMin: [1] minimum center frequency
%    fMax: [1] maximum center frequency
%    fs: [1] sampling frequency
%    fStep: [1] interval between center frequencies
%    df: [1] bandwidth around each center frequency
%    designType: 'fir1' [default], 'firls', or 'cheby2' for filter design
%    nCycle: number of cycles to include for the FIR filters (default: 3).
%            length of FIR filter is nCycle * (sampling freq)/(center freq)
%
% [b, a, fCenterList, nTap] = filterBank_simple(fMin, fMax, fStep, df, fs, designType, nCycle)
% e.g.
% [b, a, fCenterList, nTap] = filterBank_simple(8, 100, 5, 4, 1000, 'fir1', 3);
fCenterList = fMin:fStep:fMax;

if nargin < 6
    designType = 'fir1';
end

if nargin < 7
    nCycle = 3;
end

switch lower(designType)
    case {'fir1'}
        designBPhandle = @(f) fir1(ceil(nCycle/f*fs), [f-df/2, f+df/2] * 2 / fs, 'bandpass');
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
        designBPhandle = @(f) firls(ceil(nCycle/f*fs), [0, f-df, f-df/2, f+df/2, f+df, fs/2] * 2 / fs, [0 0 1 1 0 0]);
        isFIR = true;
        
    case {'cheby2'}
        %[n,Ws] = cheb2ord(Wp,Ws,Rp,Rs)
        designBPhandle = @(f) cheby2(4, 40, [f-df/2, f+df/2] * 2 / fs);
        isFIR = false;
        
    otherwise
        error('Unknown filter design [%s]', designType);
end

b = cell(numel(fCenterList),1);
a = cell(numel(fCenterList),1);
for kCenter = 1:numel(fCenterList)
    fCenter = fCenterList(kCenter);
    if isFIR
        b{kCenter} = designBPhandle(fCenter);
        a{kCenter} = 1;
    else
        [b{kCenter}, a{kCenter}] = designBPhandle(fCenter);
    end
end

% The longest filter side-effect
% boundary effect of filtering to be removed for all signals
suggestedBoundaryRemoval = ceil(nCycle/fMin*fs);