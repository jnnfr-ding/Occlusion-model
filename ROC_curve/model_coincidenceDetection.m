%{
PD = preferred direction
ND = null direction
ONR = occlusion null response
%}

function [irun_truePositiveRate, irun_falsePositiveRate] = model_coincidenceDetection(threshold, filter, filterLength, nrows, mcells, somaDistance_mean_PD, somaDistance_mean_ND, somaDistance_sd, radius_movingField, dt, barSpeed, Ktrials, baselineFR, alignPD, alignMargin, spikeResp_sd, occ_xpos, occ_tbin, prefSide_width, occ_window, noise_only, ff_cells_onset, onr_cells_onset, ff_onset_duration_corr, onr_onset_duration_corr, ff_peakFR_PD, onr_peakFR_ND, tuningWidth_gaussFit_param, gaussFit_collapsedTuningCurve_SD)

%% Create dependent parameters.
Ncells = nrows*mcells; % N cells for each subtype: N = n x m

%% Soma positions
%{ 
Set up model by creating soma positions. Right- and left-preferring 
populations are independently distributed across space, i.e. no position 
correlation between populations. The inter-soma (nearest neighbor) distance
distribution is Gaussian.
%}

% Generate soma positions for right- and left-preferring populations
somaPositions_rightPop = somaPositions(nrows,mcells,somaDistance_mean_PD,somaDistance_sd,radius_movingField*2);
somaPositions_leftPop = somaPositions(nrows,mcells,somaDistance_mean_ND,somaDistance_sd,radius_movingField*2);

% Get soma positions for each subtype as a single column vector
somaPositions_rightPop = somaPositions_rightPop';   % PD-preferring pop
somaPositions_rightPop = somaPositions_rightPop(:); % Reshape
somaPositions_leftPop = somaPositions_leftPop';     % ND-preferring pop
somaPositions_leftPop = somaPositions_leftPop(:);   % Reshape

% Total number of time bins in model
% i.e. Amount of time (in bins) for leading edge of bar to reach all somas 
% plus 1320-micron buffer distance
Tbins = ceil(1/dt * 1/barSpeed * ...
    (max([somaPositions_leftPop somaPositions_rightPop],... % Rightmost soma position
    [],'all') + radius_movingField*4/barSpeed/dt));         % 1320-micron buffer distance on far right end

%% Simulation
DSGCpopResp = sim_coincidenceDetection(1, noise_only, Tbins, somaPositions_rightPop, somaPositions_leftPop, Ncells, Ktrials, baselineFR, dt, barSpeed, alignPD, alignMargin, occ_xpos, spikeResp_sd, prefSide_width, radius_movingField, ff_cells_onset, onr_cells_onset, ff_onset_duration_corr, onr_onset_duration_corr, ff_peakFR_PD, onr_peakFR_ND, tuningWidth_gaussFit_param, gaussFit_collapsedTuningCurve_SD);

%% Decoder
%{
... threshold... detect...discriminate... ROC curve
%}

irun_truePositiveRate = zeros(1,Ktrials);
irun_falsePositiveRate = zeros(1,Ktrials);

for ktrial = 1:Ktrials                  % Loop through K trials
    % Extract firing rates from k-th trial
    DSGCpopResp_ktrial = squeeze(DSGCpopResp(:,:,ktrial));
    % Sum firing rates across 2N cells
    DSGCpopResp_ktrial = sum(DSGCpopResp_ktrial);
    % Convolve firing-rate response with temporal filter
    % Detect time bins when convolution output crosses threshold
    
    threshold_cross = decoder_coincidenceDetection(DSGCpopResp_ktrial, filter, filterLength, Tbins, threshold);
    
    irun_truePositiveRate(1,ktrial) = sum(threshold_cross(occ_tbin + occ_window)) / length(occ_window);
    
    irun_falsePositiveRate(1,ktrial) = (sum(threshold_cross(1:occ_tbin+occ_window(1)-1)) + sum(threshold_cross(occ_tbin+occ_window(end)+1:end))) / (length(threshold_cross) - length(occ_window));
end

end



