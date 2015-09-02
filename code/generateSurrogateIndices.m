function rpIdxAll = generateSurrogateIndices(T, minTimeShift, nSurrogate)
% Generate random block permutation as surrogates.
%
% rpIdxAll = generateSurrogateIndices(T, minTimeShift, nSurrogate)
%
% Generate surrogate shuffles from one continuously recorded time series
% of length T.
%
% If you have enough trials it is probably better to shuffle across trials
% given that inter-trial intervals are long enough.
%
% Input
%    T: [1] length of time series
%    minTimeShift: [1] block size to be shuffled
%    nSurrogate: [1] number of surrogates (default: 999)
%
% Output
%    rpIdxAll: [floor(T / minTimeShift) x nSurrogate+1)]
%             first column is the original order
%             the rest of the columns are randomly shuffled

nChunk = floor(T / minTimeShift); % number of chunks to be shuffled
if nChunk < 6 % 6! = 720, 7! = 5040
    error('Not enough data');
end

rpIdxAll = zeros(nChunk * minTimeShift, nSurrogate+1);
rpIdxAll(:, 1) = 1:size(rpIdxAll, 1); % first column is the original ordering
for k = 2:nSurrogate+1
    rp = randperm(nChunk);
    originalIdx = reshape(1:nChunk * minTimeShift, minTimeShift, nChunk);
    rpIdx = originalIdx(:, rp);
    rpIdxAll(:, k) = rpIdx(:);
end
