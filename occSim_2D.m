function [DSGCpopResp_subType onr_idx] = occSim_2D(PD,onr,Ncells,Tbins,Ktrials,dt,...
    radius_movingField,barSpeed,barPos_y,spikeRespSD,...
    somaPos_x,somaPos_y,occPos_x,occPos_y,prefSide_width,RFwidths,...
    ff_cells_onset,ff_onset_duration_corr,ff_peakFR_PD,...
    onr_cells_onset,onr_onset_duration_corr,onr_peakFR_ND,...
    alignPD,alignMargin,tuningWidth_gaussFit_param,gaussFit_collapsedTuningCurve_SD,...
    baselineFR)
%{
INPUTS:
    PD      = T/F   Bar moving in preferred direction (PD)
    onr     = T/F   Simulate with occlusion?
    Ncells  = 1x1
    Tbins   = 1x1
    Ktrials = 1x1
    dt      = 1x1
    
    radius_movingField  = 1x1
    barSpeed            = 1x1
    barPos_y            = 1x1
    spikeRespSD         = 
    
    somaPos_x           = Nx1
    somaPos_y           = Nx1
    occPos_x            = 
    occPos_y            =
    prefSide_width      =
    RFwidths            = Nx1
    
    ff_cells_onset          = 
    ff_onset_duration_corr  =
    ff_peakFR_PD            =
    onr_cells_onset         =
    onr_onset_duration_corr =
    onr_peakFR_ND           =
    
    alignPD                             =
    alignMargin                         = 
    tuningWidth_gaussFit_param          = 
    gaussFit_collapsedTuningCurve_SD    = 
    
    baselineFR = 1x1 Baseline firing rate (spikes/sec)

OUTPUTS:
    DSGCpopResp_subtype     = NxTxK
%}

% (0) Initialize population responses for one sub-type
DSGCpopResp_subType = zeros(Ncells,Tbins,Ktrials);
% (0) Baseline firing rate
%DSGCpopResp_subType = spikeRespSD*randn(Ncells,Tbins,Ktrials);
%DSGCpopResp_subType(DSGCpopResp_subType<0) = 0;

% Generate uniform random numbers between 0 and 1
firingProb = rand(Ncells,Tbins,Ktrials);
DSGCpopResp_subType(firingProb<baselineFR*dt) = 1/dt;

% (1) Sample N spike response onset times in a 2-second trial (Nx1)
% (2) Get corresponding spike response durations (Nx1)
[spikeRespOnset spikeRespDuration] = ...
    getOnsetDuration(ff_cells_onset,ff_onset_duration_corr,Ncells);
% (3) Convert onset times to receptive field (RF) radii (Nx1)
radius_RF = radius_movingField - spikeRespOnset*barSpeed;

if onr  % Occlusion null response
    % Indices of cells that will give ONR
    onr_x_idx = find(occPos_x > somaPos_x &...
        somaPos_x > occPos_x-prefSide_width);
    onr_y_idx = find(occPos_y-100 < somaPos_y &...
        somaPos_y < occPos_y+100);

    onr_idx = intersect(onr_x_idx,onr_y_idx);
    somaPos_x(onr_idx) = occPos_x - prefSide_width; 
    rcells = length(onr_idx);          % r occluded cells give ONR
    
    % (1) and (2) for ONR
    [spikeRespOnsetONR spikeRespDurationONR] = ... % (rx1)
        getOnsetDuration(onr_cells_onset,onr_onset_duration_corr,rcells);
    radius_RF(onr_idx) = spikeRespOnsetONR*barSpeed - radius_movingField;
    spikeRespOnset(onr_idx) = spikeRespOnsetONR;
    spikeRespDuration(onr_idx) = spikeRespDurationONR;
end

% (4) Convert RF radii to time bins at which spike responses start (1xN)
tbin_respStarts = ceil((somaPos_x - radius_RF)'./barSpeed./dt);

% (5) Sample N peak firing rates (Nx1)
if PD   % PD-preferring cells
    spikeRespPeakFR = datasample(ff_peakFR_PD,Ncells);
else    % ND-preferring cells, zero mean
    spikeRespPeakFR = zeros(1,Ncells);
end

if onr
    % For occlusion null response, sample from ONR-specific distribution
    % and adjust start time bin
    spikeRespPeakFR(onr_idx) = datasample(onr_peakFR_ND,rcells);
    tbin_respStarts(onr_idx) = ceil((somaPos_x(onr_idx) + radius_RF(onr_idx))'./barSpeed./dt);
end

if ~alignPD     % Jitter in PD axis alignment?
    spikeRespPeakFR = jitterPDalignment(Ncells,alignMargin,spikeRespPeakFR,...
        tuningWidth_gaussFit_param,gaussFit_collapsedTuningCurve_SD);
end

% (6)* Change peak FR according to somaPos_y!!!!!
% Nx1
attenuate_y_offset = normpdf(barPos_y,somaPos_y,RFwidths)...    % Gaussian RF
    ./normpdf(0,0,RFwidths);                                    % Normalize                                      
spikeRespPeakFR = spikeRespPeakFR.*attenuate_y_offset';         % Scale peak FR

% Parameters for sine fit
sinFit_param = NaN(Ncells,3);               % Initialize
sinFit_param(:,1) = spikeRespPeakFR';       % Amplitude
sinFit_param(:,2) = spikeRespDuration.*2;   % Period
sinFit_param(:,3) = spikeRespOnset;         % Phase
for ktrial = 1:Ktrials
    % Introduce trial-to-trial variability
    for icell = 1:Ncells
        % (7) Create sine fit function handle
        sinFit = createSinFitHandle(sinFit_param(icell,:));
        % (8) Time bins with nonzero response
        tbins_nonzeroResp = spikeRespOnset(icell):dt:...
            spikeRespOnset(icell)+spikeRespDuration(icell)+dt;
        % (9) Evaluate sine fit at those time bins to get mean response
        DSGCmeanResp = sinFit(tbins_nonzeroResp);
        % (10) Add sub-Poisson noise to sinusoidal mean response    
        DSGCresp = spikeRespSD.*randn(1,length(DSGCmeanResp))+DSGCmeanResp;
        DSGCresp(DSGCresp<0) = 0;     % Rectify
        % (11) Paste in response  
        DSGCpopResp_subType(icell,tbin_respStarts(icell):...
            tbin_respStarts(icell)+length(tbins_nonzeroResp)-1,...
            ktrial) = DSGCresp;
    end
end

end




function [onset,duration,type] = getOnsetDuration(onsetTimes,corrOnsetDuration,N)
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
[onset,type] = pearsrnd(muOnset,sigmaOnset,skewOnset,kurtOnset,N,1);
% Spike response durations
duration = slope.*onset+intercept;

end

function fr = jitterPDalignment(Ncells,alignMargin,spikeRespPeakFR,...
    tuningWidth_gaussFit_param,gaussFit_collapsedTuningCurve_SD)
% Deviations from PD (1xN)
alignJitter = 2*alignMargin*rand(1,Ncells)-alignMargin; % Sample from uniform distribution
% Sample circular variances (1xN)
alignCircVar = normrnd(tuningWidth_gaussFit_param(1),...
    tuningWidth_gaussFit_param(2),1,Ncells);
% Convert circular variances to angular deviations (1xN)
alignAngDev = sqrt(2*alignCircVar);
% Tuning curve widths (1xN)
tuningWidths = alignAngDev.*gaussFit_collapsedTuningCurve_SD; % SD of direction tuning curve
% New peak firing rates after jitter (1xN)
fr = gaussTuningCurve(alignJitter,spikeRespPeakFR,tuningWidths);

end

function fr = gaussTuningCurve(x,a,sigma)
%{
Find modulated/attenuated peak firing rates due to jittered PD alignment
from von Mises direction tuning curve approximated as Gaussian.
%}
% INPUTS:
%   x     = 1xN Jitters from preferred direction (deg)
%   a     = 1xN Peak firing rates, which are the amplitudes of the tuning
%               curves (Hz)
%   sigma = 1xN Standard deviations of the tuning curves (deg)
% OUTPUT:
%   fr     = 1xN Firing rates at jittered directions (Hz)

fr = a.*exp(-0.5*(x./sigma).^2);

end

function y = createSinFitHandle(p)
%{
Create function handle for sine fit to spike response. Amplitude is given 
by peak firing rate; period is given by twice the spike response duration; 
phase is given by spike response onset time.
%}
% INPUTS:
%   p(1) = amplitude
%   p(2) = period
%   p(3) = phase
% OUTPUT:
%   y    = function handle for sine fit

y = @(t) p(1).*sin(2.*pi./p(2).*(t - p(3)));

end

