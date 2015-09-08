function [stat] = surrogateStats(observedSample, surrogateSamples)
% Extract statistics from surrogate distribution
% stat = surrogateStats(observedSample, surrogateSamples)
%
% Input
%   observedSample: [1] observed value
%   surrogateSamples: [nSurrogate x 1] surrogate samples from null hypothesis
%
% Output
%   stat.value: [1] observedSample itself is returned here
%   stat.pValue: [1] estimated p-value
%   stat.deviation: [1] deviation from the mean
%   stat.smean: [1] mean of the surrogates
%   stat.sstd: [1] standard deviation of the surrogates
%   stat.thresholds: [4x1] estimated quantiles of the surrogate distribution
%                    at alpha = 0.01, 0.05, 0.1, 0.5

if nargin < 1
    stat = struct('value', 0, 'pValue', 0, 'deviation', 0, 'smean', 0, 'sstd', 0, 'thresholds', []);
    return
end

nSurrogate = numel(surrogateSamples);
pValue = (sum(surrogateSamples > observedSample) + 1)/(nSurrogate+1); % small bias away from 0
m = mean(surrogateSamples);
s = std(surrogateSamples);
dev = (observedSample - m)/s;

stat.value = observedSample;
stat.pValue = pValue;
stat.deviation = dev;
stat.smean = m;
stat.sstd = s;
stat.thresholds = quantile(surrogateSamples, 1 - [0.01, 0.05, 0.1, 0.5]);
