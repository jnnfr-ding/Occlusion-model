function [onset, duration, type] = getOnsetDuration(onsetTimes,corrOnsetDuration,N)
%{
Sample spike response onset time from a Pearson distribution. Get spike 
response duration from correlation between spike response onset and 
duration.
%}
% INPUTS:
%   onsetTimes        =  Array of spike response onset times from 
%                        experiment (sec)
%   corrOnsetDuration =  Correlation between spike response onset and 
%                        duration
%   N                 =  Number of data points to generate
% OUTPUTS:
%   onset    =  Nx1 Sampled spike response onset times (sec)
%   duration =  Nx1 Spike response durations corresponding to sampled onset
%                   times according to linear correlation (sec)
%   type     =  Nx1 Type of distribution in the Pearson system

% Get moments and correlation
muOnset = nanmean(onsetTimes,'all');       % Mean onset time
sigmaOnset = nanstd(onsetTimes,0,'all');   % Standard deviation
skewOnset = skewness(onsetTimes,0,'all');  % Skewness (3rd moment)
kurtOnset = kurtosis(onsetTimes,0,'all');  % Kurtosis (4th moment)
slope = corrOnsetDuration(1);
intercept = corrOnsetDuration(2);

% Spike response onset times
[onset, type] = pearsrnd(muOnset,sigmaOnset,skewOnset,kurtOnset,N,1);
% Spike response durations
duration = slope.*onset+intercept;

end