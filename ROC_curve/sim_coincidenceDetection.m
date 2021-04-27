% Simulate an experiment with a single long vertical occlusion while the
% bar moves from left to right.

function DSGCpopResp = sim_coincidenceDetection(occlusion, noise_only, Tbins, somaPositions_rightPop, somaPositions_leftPop, Ncells, Ktrials, baselineFR, dt, barSpeed, alignPD, alignMargin, occ_xpos, spikeResp_sd, prefSide_width, radius_movingField, ff_cells_onset, onr_cells_onset, ff_onset_duration_corr, onr_onset_duration_corr, ff_peakFR_PD, onr_peakFR_ND, tuningWidth_gaussFit_param, gaussFit_collapsedTuningCurve_SD)

% Initialize population responses
DSGCpopResp_PD = zeros(Ncells,Tbins,Ktrials); % PD-preferring pop
DSGCpopResp_ND = zeros(Ncells,Tbins,Ktrials); % ND-preferring pop
% DSGCpopResp_PD_noise = zeros(Ncells,Tbins,Ktrials); % PD-preferring pop noise
% DSGCpopResp_ND_noise = zeros(Ncells,Tbins,Ktrials); % ND-preferring pop noise

% if baselineFR > 0
%     % Generate uniform random numbers between 0 and 1
%     firingProb = rand(Ncells,Tbins,Ktrials);            % Generate probabilities
%     % Probability of firing is approximated by rate x time
%     DSGCpopResp_PD_noise(firingProb<baselineFR*dt) = 1/dt;    %spike/sec
% 
%     firingProb = rand(Ncells,Tbins,Ktrials);            % Generate probabilities
%     DSGCpopResp_ND_noise(firingProb<baselineFR*dt) = 1/dt;    %spike/sec
% end

DSGCpopResp_PD_noise = 0.5*baselineFR + baselineFR*rand(Ncells,Tbins,Ktrials);
DSGCpopResp_ND_noise = 0.5*baselineFR + baselineFR*rand(Ncells,Tbins,Ktrials);


%*% Right-preferring pop (PD)
if ~noise_only
    % (1) Sample N spike response onset times in a 2-second trial (Nx1)
    % (2) Get corresponding spike response durations (Nx1)
    [spikeRespOnset, spikeRespDuration] = getOnsetDuration(ff_cells_onset,ff_onset_duration_corr,Ncells);
    % (3) Convert onset times to receptive field (RF) radii (Nx1)
    radius_RF = radius_movingField - spikeRespOnset*barSpeed;
    % (4) Convert RF radii to time bins at which spike responses start (1xN)
    tbin_respStarts_rightPop = ceil((somaPositions_rightPop - radius_RF)' ./ barSpeed ./ dt);
    % (5) Sample N peak firing rates (Nx1)
    spikeRespPeakFR = datasample(ff_peakFR_PD,Ncells);
    if ~alignPD     % Jitter in PD axis alignment?
        % Deviations from PD (1xN)
        alignJitter = 2*alignMargin*rand(1,Ncells)-alignMargin; % Sample from uniform distribution
        % Sample circular variances (1xN)
        alignCircVar = normrnd(tuningWidth_gaussFit_param(1),tuningWidth_gaussFit_param(2),1,Ncells);
        % Convert circular variances to angular deviations (1xN)
        alignAngDev = sqrt(2*alignCircVar);
        % Tuning curve widths (1xN)
        tuningWidths = alignAngDev.*gaussFit_collapsedTuningCurve_SD; % SD of direction tuning curve
        % New peak firing rates after jitter (1xN)
        spikeRespPeakFR = gaussTuningCurve(alignJitter,spikeRespPeakFR,tuningWidths);
    end
    % Parameters for sine fit
    sinFit_param = NaN(Ncells,3);               % Initialize
    sinFit_param(:,1) = spikeRespPeakFR';       % Amplitude
    sinFit_param(:,2) = spikeRespDuration.*2;   % Period
    sinFit_param(:,3) = spikeRespOnset;         % Phase
    for ktrial = 1:Ktrials
        % Introduce trial-to-trial variability
        for icell = 1:Ncells
            % (6) Create sine fit function handle
            sinFit = createSinFitHandle(sinFit_param(icell,:));
            % (7) Time bins with nonzero response
            tbins_nonzeroResp = spikeRespOnset(icell):dt:spikeRespOnset(icell)+spikeRespDuration(icell)+dt;
            % (8) Evaluate sine fit at those time bins to get mean response
            DSGCmeanResp_fullfield_PD = sinFit(tbins_nonzeroResp);
            % (9) Add sub-Poisson noise to sinusoidal mean response    
            DSGCresp_fullfield_PD = spikeResp_sd.*randn(1,length(DSGCmeanResp_fullfield_PD))+DSGCmeanResp_fullfield_PD;
            DSGCresp_fullfield_PD(DSGCresp_fullfield_PD<0) = 0;     % Rectify
            % (10) Paste in response  
            DSGCpopResp_PD(icell, tbin_respStarts_rightPop(icell):tbin_respStarts_rightPop(icell)+length(tbins_nonzeroResp)-1, ktrial) = DSGCresp_fullfield_PD;
        end
    end



    %*% Left-preferring pop (ND)
    if occlusion
        %*% Occlusion null
        % Indices of cells that will give ONR
        onr_cell_idx = find(occ_xpos > somaPositions_leftPop &...
            somaPositions_leftPop > occ_xpos-prefSide_width);

        % Synchronize the responses of all the occluded cells
        % Make their soma positions identical
        % because spike response onset (latency) depends on soma position
        somaPositions_leftPop(onr_cell_idx) = occ_xpos - prefSide_width; 

        rcells = length(onr_cell_idx);          % r occluded cells give ONR
        % (0) Initialize
%          tbin_respStarts_leftPop = zeros(1,rcells);  % Time bin at which spike responses start (temporary)!!!
        % (1) Sample a spike response onset times in a 2-second trial (rx1)
        % (2) Get corresponding spike response duration (rx1)
        [spikeRespOnset, spikeRespDuration] = ...
            getOnsetDuration(onr_cells_onset,onr_onset_duration_corr,rcells);
        % (3) Convert onset times to receptive field (RF) radius (rx1)
        % Note: Different for occlusion null response!!!
        radius_RF = spikeRespOnset*barSpeed - radius_movingField; 
        % (4) Convert RF radius to time bin at which spike response start (rx1)
        tbin_respStarts_leftPop = ceil((somaPositions_leftPop(onr_cell_idx)...
            + radius_RF)'./barSpeed./dt);
        % (5) Sample a peak firing rate (rx1)
        spikeRespPeakFR = datasample(onr_peakFR_ND,rcells);
        % Parameters for sine fit
        sinFit_param = NaN(rcells,3);               % Initialize
        sinFit_param(:,1) = spikeRespPeakFR';       % Amplitude
        sinFit_param(:,2) = spikeRespDuration.*2;   % Period
        sinFit_param(:,3) = spikeRespOnset;         % Phase
        for ktrial = 1:Ktrials
            % Introduce trial-to-trial variability
            for icell = 1:rcells
                idx = onr_cell_idx(icell);          % Index of occluded cell
                % (6) Create sine fit function handle
                sinFit = createSinFitHandle(sinFit_param(icell,:));
                % (7) Time bins with nonzero response
                tbins_nonzeroResp = spikeRespOnset(icell):dt:...
                   spikeRespOnset(icell)+spikeRespDuration(icell)+dt;
                % (8) Evaluate sine fit at those time bins to get mean response
                DSGCmeanResp_occlusion_ND = sinFit(tbins_nonzeroResp);
                % (9) Add sub-Poisson noise to sinusoidal mean response    
                DSGCresp_occlusion_ND = spikeResp_sd.*randn(1,length(DSGCmeanResp_occlusion_ND))...
                    +DSGCmeanResp_occlusion_ND;
                DSGCresp_occlusion_ND(DSGCresp_occlusion_ND<0) = 0;     % Rectify
                % (10) Paste in response  
                DSGCpopResp_ND(idx,tbin_respStarts_leftPop(icell):...
                tbin_respStarts_leftPop(icell)+length(tbins_nonzeroResp)-1,...
                ktrial) = DSGCresp_occlusion_ND;
            end
        end
    else
        %*% Full-field null
        % (1) Sample N spike response onset times in a 2-second trial (Nx1)
        % (2) Get corresponding spike response durations (Nx1)
        [spikeRespOnset, spikeRespDuration] = ...
            getOnsetDuration(ff_cells_onset,ff_onset_duration_corr,Ncells);
        % (3) Convert onset times to receptive field (RF) radii (Nx1)
        radius_RF = radius_movingField - spikeRespOnset*barSpeed;
        % (4) Convert RF radii to time bins at which spike responses start (1xN)
        tbin_respStarts_leftPop = ceil((somaPositions_leftPop - radius_RF)' ./ barSpeed ./ dt);
        for ktrial = 1:Ktrials
            % Introduce trial-to-trial variability
            for icell = 1:Ncells
                % (7) Time bins with nonzero response
                tbins_nonzeroResp = spikeRespOnset(icell):dt:...
                    spikeRespOnset(icell)+spikeRespDuration(icell)+dt;
                % (9) Add sub-Poisson noise to mean zero response    
                DSGCresp_fullfield_ND = spikeResp_sd.*randn(1,length(tbins_nonzeroResp));
                DSGCresp_fullfield_ND(DSGCresp_fullfield_ND<0) = 0;     % Rectify
            % (10) Paste in response  
                DSGCpopResp_ND(icell,tbin_respStarts_leftPop(icell):tbin_respStarts_leftPop(icell)+length(tbins_nonzeroResp)-1,ktrial) = DSGCresp_fullfield_ND;
            end
        end
    end
end

DSGCpopResp_PD = DSGCpopResp_PD + DSGCpopResp_PD_noise;
DSGCpopResp_ND = DSGCpopResp_ND + DSGCpopResp_ND_noise;

% Population response of all cells (2N x T x K)
DSGCpopResp = [DSGCpopResp_PD; DSGCpopResp_ND];


end


