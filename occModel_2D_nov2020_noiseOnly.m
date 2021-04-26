%{ 
Change filename fn
Change        onr = 0 (Occlusion null response) or 1 (Full-field) (l.15)
Change barSpeed   = [330,660,1320,1980,2640]                      (l.53)
Change baselineFR = [0,0.025,0.1]                                 (l.75)
Change occluder position occPos_x = [1100,1800,2500,3200]         (l.26)
Change variable name on l.235 to match fn
%}

%fn = 'ff_speed0330_bkgdFR0000';             % Full-field
%fn = 'occ_speed2640_bkgdFR0000_occPos3200'; % Occlusion
plotRMSE = 1;           % Plot RMSE?
saveRMSE = 0;           % Save RMSE?

onr = 0;        % USE FULL-FIELD FOR NOISE ONLY SIMULATIONS

Rruns = 100;           % L different soma arrangements i.e. retinas
Ktrials = 10;           % K trials or repetitions for each soma arrangement

nrows = 10;             % n 1D arrays
mcells = 100;           % m cells per array for each subtype
Ncells = nrows*mcells;  % N cells for each subtype: N = n x m

barPos_y = 800;
occPos_y = barPos_y;
occPos_x = 1800;        % Fix this magic number (microns)!!!
                        % OR specify after setting up soma positions

dt = 0.010;             % Time step step size (sec)
t_axis = 0:dt:50;       % time axis (sec)

somaDist_mean_PD = 39;  % Inter-soma distance mean for PD-preferring (microns)
somaDist_mean_ND = 39;  % Inter-soma distance mean for ND-preferring (microns)
somaDist_SD = 16;       % Inter-soma distance standard deviation (microns)

% Spatial set-up
dx = 1;                 % Spatial discretization/resolution (microns)
x_axis = 0:dx:50000;     % x-axis

prefSide_width = 110;   % Width of region in which cell is sensitive to
                        % emerging edge moving in ND (microns)
                        
radius_movingField = 330;   % Moving field radius (micron)
radius_dendriticField_mean = 88;    % Dendritic field radius mean (microns)
radius_dendriticField_SD = 14.8;    % Dendritic field radius standard deviation
radius_occStim = 110;               % Occlusion stimulus radius (micron)
dendriticField2RF = 1.25;           % Scale dendritic field radius to decoder SD
dendriticField2onrRF = 0.25;        % Scale dendritic field radius to decoder SD for occluded cell
integrationTime = 0.070;            % Integration time (sec)
RFcenterOffset = 150;       % Spatial offset between soma and functional RF center (microns)

%*% Bar speed
barSpeed = 330;                % micron/sec
barSpeed_baseline = 330;        % micron/sec

onset_speed_slope_ff  = 0.095;  % sec
onset_speed_slope_onr = 0.083;  % sec

%*% Position of leading edge of bar moving rightward from 0 to 6000 microns
%*% Note there is at least 330 microns before bar reaches first soma center
barPos_x = 0:dt*barSpeed:x_axis(end);
speed_duration_scalingFactor_onr = -0.561; % Power-law exponent for ONR

%*% Jitter along PD-ND axis
alignPD = 0;        % Align preferred direction (PD) axes of all DSGCs?
alignMargin = 14.1; % PD +/- ___ degrees

%*% Sub-Poisson noise model 
spikeRespVar = 0.4;                  % Constant variance
spikeRespSD = sqrt(spikeRespVar);    % Standard deviation

%*% Baseline firing rate
baselineFR = 0.025;          % spikes/sec
%baselineFR = 4;             % Hz/bin, i.e. spikes/sec/dt
decodeRange = 400;          % microns

optical_conversion = 33;    % Conversion from visual angle subtended on
                            % retina to microns (microns per deg)

                            
%% Load data files
cd('C:\Users\alber\Dropbox\Occlusion Project\sampling_distributions')
% Load spike response onset times
load('ff_cells_onset')
load('occ_cells_onset')
load('onr_cells_onset')
% Load spike response onset-duration correlation
load('ff_onset_duration_corr')
load('occ_onset_duration_corr')
load('onr_onset_duration_corr')
% Load spike response peak firing rates
load('ff_peakFR_PD')
load('occ_peakFR_PD')
load('onr_peakFR_ND')
% Load direction tuning curve parameters
load('tuningWidth_gaussFit_param')
load('gaussFit_collapsedTuningCurve_SD')
cd('C:\Users\alber\Dropbox\Occlusion Project\occlusion_model_nov2020')


%% Model
posEstRMSE = NaN(Rruns,length(t_axis));     % Initialize
for rrun = 1:Rruns
    %*% SOMA POSITIONS
    % Generate soma positions
    somaPos_x_rightPop = somaPositions(...  % PD-preferring pop x-coordinates
        nrows,mcells,somaDist_mean_PD,somaDist_SD,2*radius_movingField);
    somaPos_y_rightPop = somaPositions(...  % PD-preferring pop y-coordinates
        mcells,nrows,somaDist_mean_PD,somaDist_SD,2*radius_movingField)';
    somaPos_x_leftPop = somaPositions(...   % ND-preferring pop x-coordinates
        nrows,mcells,somaDist_mean_ND,somaDist_SD,2*radius_movingField);
    somaPos_y_leftPop = somaPositions(...   % ND-preferring pop y-coordinates
        mcells,nrows,somaDist_mean_ND,somaDist_SD,2*radius_movingField)';
    % Reshape into single column vector (Nx1)
    somaPos_x_rightPop = somaPos_x_rightPop';   
    somaPos_x_rightPop = somaPos_x_rightPop(:); % PD-preferring pop x-coordinates
    somaPos_y_rightPop = somaPos_y_rightPop';
    somaPos_y_rightPop = somaPos_y_rightPop(:); % PD-preferring pop y-coordinates
    somaPos_x_leftPop = somaPos_x_leftPop';
    somaPos_x_leftPop = somaPos_x_leftPop(:);   % ND-preferring pop x-coordinates
    somaPos_y_leftPop = somaPos_y_leftPop';
    somaPos_y_leftPop = somaPos_y_leftPop(:);   % ND-preferring pop y-coordinates
    
    % Total number of time bins in model
    % i.e. Amount of time (in bins) for leading edge of bar to reach all somas 
    % plus 1320-micron buffer distance
    Tbins = ceil(1/dt * 1/barSpeed * (max([somaPos_x_rightPop; somaPos_x_leftPop]))...         % Rightmost soma position
        + radius_movingField*4/barSpeed/dt);   % 1320-micron buffer distance
    
    %*% SIMULATIONS   
    if onr  % Occlusion
        PD = 1;     % Right-preferring pop
        RFwidths_PD = receptiveFieldWidth(Ncells,...
            radius_dendriticField_mean,radius_dendriticField_SD,dendriticField2RF);
%        DSGCpopResp_PD = 0.5*baselineFR + baselineFR*rand(Ncells,Tbins,Ktrials);
        PD = 0;     % Left-preferring pop
        RFwidths_ND = receptiveFieldWidth(Ncells,...
            radius_dendriticField_mean,radius_dendriticField_SD,dendriticField2RF);
%        [DSGCpopResp_ND,onr_idx] = ...;
    else    % Full-field
        PD = 1;     % Right-preferring pop
        RFwidths_PD = receptiveFieldWidth(Ncells,...
            radius_dendriticField_mean,radius_dendriticField_SD,dendriticField2RF);
        
        DSGCpopResp_PD = zeros(Ncells,Tbins,Ktrials);
        firingProb = rand(Ncells,Tbins,Ktrials);
        DSGCpopResp_PD(firingProb<baselineFR*dt) = 1/dt;
%        DSGCpopResp_PD = 0.5*baselineFR + baselineFR*rand(Ncells,Tbins,Ktrials);
        PD = 0;     % Left-preferring pop
        RFwidths_ND = receptiveFieldWidth(Ncells,...
            radius_dendriticField_mean,radius_dendriticField_SD,dendriticField2RF);
        
        DSGCpopResp_ND = zeros(Ncells,Tbins,Ktrials);
        firingProb = rand(Ncells,Tbins,Ktrials);
        DSGCpopResp_ND(firingProb<baselineFR*dt) = 1/dt;
%        DSGCpopResp_ND = 0.5*baselineFR + baselineFR*rand(Ncells,Tbins,Ktrials);
    end
    % Population response of all cells (2N x T x K)
    DSGCpopResp = [DSGCpopResp_PD; DSGCpopResp_ND];

    
    %*% DECODING
    somaPos_x_rightPop = somaPos_x_rightPop - RFcenterOffset ...
        + integrationTime*barSpeed_baseline;
    somaPos_x_leftPop  = somaPos_x_leftPop + RFcenterOffset... 
        + integrationTime*barSpeed_baseline;
    if onr
        somaPos_x_leftPop(onr_idx) = occPos_x - prefSide_width;
        somaPos_x_leftPop(onr_idx) = somaPos_x_leftPop(onr_idx)... 
            + RFcenterOffset + integrationTime*barSpeed_baseline;        % Delay in ONR
        RFwidths_ND(onr_idx) = dendriticField2onrRF.*RFwidths_ND(onr_idx)./dendriticField2RF;
    end
    
    [posEst_x,somaPos_noisy_x] = popVectorDecoder_nov2020(...
        DSGCpopResp,[somaPos_x_rightPop; somaPos_x_leftPop],[RFwidths_PD; RFwidths_ND]);

    [posEst_y,somaPos_noisy_y] = popVectorDecoder_nov2020(...
        DSGCpopResp,[somaPos_y_rightPop; somaPos_y_leftPop],[RFwidths_PD; RFwidths_ND]);

%*% Decoder with moving window
%     [posEst_x,somaPos_noisy_x] = decoder_backgroundNoise(...
%         DSGCpopResp,[somaPos_x_rightPop; somaPos_x_leftPop],[RFwidths_PD; RFwidths_ND],...
%         barPos_x,decodeRange);
% 
%     [posEst_y,somaPos_noisy_y] = decoder_backgroundNoise(...
%         DSGCpopResp,[somaPos_y_rightPop; somaPos_y_leftPop],[RFwidths_PD; RFwidths_ND],...
%         repmat(barPos_y,1,Tbins),decodeRange); 
    
    posEstRMSE(rrun,1:Tbins) = sqrt(nanmean(...
        (posEst_x - repmat(barPos_x(1:Tbins),[Ktrials,1])).^2 + (posEst_y - barPos_y).^2));
    
    disp(rrun)
end

posEstRMSE_visualAngle = posEstRMSE./optical_conversion;

occPos_x_t = occPos_x/barSpeed;     % sec

if plotRMSE
    figure; hold on
    plot(t_axis,nanmean(posEstRMSE_visualAngle))
    line([occPos_x_t occPos_x_t],[-5 20],'Color','k','LineStyle','-.')
    xlabel('Time (s)')
    ylabel('RMSE (deg)')
    xlim([0 30])
    ylim([0 100])
end

% Make variable name same as filename fn at the top
ff_speed0330_bkgdFR0050 = nanmean(posEstRMSE_visualAngle);

if saveRMSE
    cd('C:\Users\alber\Dropbox\Occlusion Project\occlusion_model_nov2020\FigsForFred')
    save(fn,fn)
end

% figure; hold on
% viscircles([reshape(somaPos_x_rightPop,[],1)...
%     reshape(somaPos_y_rightPop,[],1)],...
%     8*ones(Ncells,1),'Color','r')  % PD is red
% viscircles([reshape(somaPos_x_leftPop,[],1)...
%     reshape(somaPos_y_leftPop,[],1)],...
%     88*ones(Ncells,1),'Color','b')  % ND is blue



