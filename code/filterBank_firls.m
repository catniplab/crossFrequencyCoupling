% Generate a filter bank

fs = 1000; % sampling frequency
fMin = 8;
fMax = 100; %200;
fStep = 5; %2;
df = 4; % bandwidth
fCenterList = fMin:fStep:fMax;

% Linear phase with least squares over side lobes of the shape:
%
%              ------------------
%
%
%    -------??                    ??----------
%
%    0  f-df    f-df/2  f   f+df/2   f+df fs/2
designBPhandle = @(f) firls(ceil(2/f*fs), [0, f-df, f-df/2, f+df/2, f+df, fs/2] * 2 / fs, [0 0 1 1 0 0]);

for kCenter = 1:numel(fCenterList)
    fCenter = fCenterList(kCenter);
    b{kCenter} = designBPhandle(fCenter);
end

% The longest filter side-effect
% boundary effect of filtering to be removed for all signals
nTap = ceil(2/fMin*fs);
removeBoundaryHandle = @(z) z(nTap+1:end-nTap);
