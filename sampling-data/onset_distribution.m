%{
(1) Find onset times for each cell and build distribution of onset times
(2) Find offset times for each cell
(3) Find spike response duration for each cell
(4) Plot correlation between onset time and spike response duration
%}

%% Parameters
% List of full-field cell indices
ff_cells_idx = [1:18 20:25 32:47 49:52 54:60 63:64 66:79 82:83 87:90];
                                        % Exclude Cell#19 !!!
ff_Ncells = length(ff_cells_idx);       % Number of full-field cells
% List of occlusion cell indices
occ_cells_idx = [1:16 18 21 23 26:45 47:48 53 55 58:64 66 68 71:73 76 79:91];
                                        % Exclude Cell#46 !!!
occ_Ncells = length(occ_cells_idx);     % Number of occlusion cells
Ncells = max(max(ff_cells_idx),max(occ_cells_idx)); % Total number of cells

num_dir = 8;           % Number of directions

binSize = 0.025;       %sec
startTime = 0;         %sec
endTime = 2;           %sec
edges = startTime:binSize:endTime; % Create time bins
Nbins = length(edges); % Number of bins

threshold_param = 4;        % Threshold is ___ s.d. above baseline
barspeed = 300;             %pixels/sec
movingField_diameter = 600; %pixels
movingField_radius = movingField_diameter/2; %pixels

binranges_onset = startTime:binSize:endTime;  % Edges for binning spike onset times
%binranges_onset = edges;
binranges_RF = 0:20:movingField_radius+5;   % Edges for binning receptive field radii

%{
Use data from bar speed 300 pixels/sec to get baseline firing rate
There is no firing (zero response) when no stimulus is present
%}
noStimStartTime = 5;        %sec
noStimEndTime = 6;          %sec

manualAdjust = 1;   % Manually adjust? (Usually yes)

% Plotting
plot_ff = 0;        % Plot figures for full-field?
plot_occ_pd = 0;    % Plot figures for occlusion PD?
plot_onr = 0;       % Plot figures for occlusion null response?

% Save onset times, offset times, and onset-duration correlation?
save_ff = 0;
save_occ_pd = 0;
save_onr = 0;

% Fitting
smoothSpline = 0;   % Smooth using cubic spline?
smoothParam = 0.9995;       % Smoothing parameter for PD and PD-adjacent directions
smoothParam_onr = 0.99998;  % Smoothing parameter for occlusion null response
fitGaussian = 0;    % Fit Gaussian to histograms?
fitDist = 1;        % (1) Fit beta distribution to histogram for PD and PD-adjacent directions
                    % (2) Fit gamma distribution to histogram for null and null-adjacent directions


%% Load file
% Direction index of each cell's PD
cd('C:\Users\alber\Dropbox\Occlusion Project\sampling_distributions') % Specify directory
load('all_cells_prefdir')   % Load PD index file
% First column for full-field
% Second column for occlusion


%% Full-field
% Direction index of each cell's PD
dir_idx_PD = all_cells_prefdir(:,1)/45+1; % Convert directions in degrees to indices 

% Initialize
% All 8 directions
ff_cells_onset_alldir = NaN(num_dir,Ncells);  % Onset
ff_cells_offset_alldir = NaN(num_dir,Ncells); % Offset
% PD and PD-adjacent directions
ff_cells_onset = NaN(3,Ncells);               % Onset
ff_cells_offset = NaN(3,Ncells);              % Offset
for icell = ff_cells_idx
    %---------------------------------------------------------------------%
    % PD and adjacent directions of i-th cell
    icell_dir_idx_PD = dir_idx_PD(icell);   % Direction index of i-th cell's PD
    PD_leftNeighbor = icell_dir_idx_PD-1;   % Neighbor left to PD
    PD_rightNeighbor =  icell_dir_idx_PD+1; % Neighbor right to PD
    if PD_leftNeighbor == 0
        % If PD = 0 degrees, set left neighbor to be 315 degrees
        PD_leftNeighbor = num_dir;
    end
    if PD_rightNeighbor > num_dir
        % If PD = 315 degrees, set right neighbor to be 0 degrees
        PD_rightNeighbor = 1;
    end
    % Array of direction indices for PD and adjacent directions
    PD_neighbors = [PD_leftNeighbor icell_dir_idx_PD PD_rightNeighbor];
    
    %---------------------------------------------------------------------%
    % Load data
    fn = strcat('fullfield_cell',num2str(icell)); % Filename
    filepath = strcat('C:\Users\alber\Dropbox\Occlusion Project\sampling_distributions\',fn);
    load(filepath) % Load data
    
    data = SortedDSData.total;                % Get data structure
    ktrials = size(data.spikeCounts,2);       %repetitions
    spikeTimes = data.spikeTimes;             % Spike times for all directions
    
    %---------------------------------------------------------------------%
    % Get baseline firing rate and threshold
    noStimFR = NaN(ktrials,num_dir);       % Initialize    
    for mdir = 1:num_dir
        % Loop through all directions
        dirSpikeTimes = spikeTimes{1,mdir}; % Get spike times for one direction
        for ktrial = 1:ktrials
            % Loop through all trials
            % Spike times for k-th trial
            ktrialSpikeTimes = dirSpikeTimes{ktrial,1};
            % Firing rate when NO stimulus is present
            noStimFR(ktrial,mdir) = max(histc(ktrialSpikeTimes,[noStimStartTime noStimEndTime]));
        end
    end
    % Threshold for change-point detection
    baseline = mean(noStimFR,'all'); % Baseline firing rate 
                                     % Averaged across all trials and all directions
    threshold = threshold_param*std(noStimFR,0,'all') + baseline;

    %---------------------------------------------------------------------%
    % Find onset and offset times
    ff_icell_onset_alldir = NaN(num_dir,1);         % Initialize onset array
    ff_icell_offset_alldir = NaN(num_dir,1);        % Initialize offset array
    for mdir = 1:num_dir
        % Loop through all directions
        dirSpikeTimes = spikeTimes{1,mdir};         % Get spike times for one direction
        dirSpikeCount = NaN(ktrials,length(edges)); % Initialize
        for ktrial = 1:ktrials
            % Loop through all trials
            % Spike times for k-th trial
            ktrialSpikeTimes = dirSpikeTimes{ktrial,1};
            % Bin spike times for k-th trial
            dirSpikeCount(ktrial,:) = histc(ktrialSpikeTimes,edges);
        end
        % Trial-averaged firing rate across bins for one direction
        dirFR = mean(dirSpikeCount)/binSize;  %spikes/sec
        % Find time bin at which response starts or ends
        % Onset
        for ibin = 1:Nbins
            % Find time bin at which response starts
            % Onset defined as above-threshold firing for 2 consecutive time bins
            if dirFR(ibin)>threshold && dirFR(ibin+1)>threshold
                ff_icell_onset_alldir(mdir) = edges(ibin);
                break % Exit for-loop once onset has been found
            end
        end
        % Offset
        for ibin = Nbins:-1:3                   % Loop backwards
            % Offset defined as above-threshold firing in 2 out of 3
            % consecutive time bins
            if dirFR(ibin)>threshold            % Detect above-threshold firing in first time bin
                if dirFR(ibin-1)>threshold      % Detect above-threshold firing in second time bin
                    ff_icell_offset_alldir(mdir) = edges(ibin);
                    break       % Exit for-loop once offset has been found
                else % If second time bin has zero response
                    if dirFR(ibin-2)>threshold  % Detect above-threshold firing in third time bin
                        ff_icell_offset_alldir(mdir) = edges(ibin);
                        break   % Exit for-loop once offset has been found
                    end
                end
            end
        end
    end
    % Store spike onset and offset times
    % Onset
    ff_cells_onset_alldir(:,icell) = ff_icell_onset_alldir;          % All 8 directions
    ff_cells_onset(:,icell) = ff_icell_onset_alldir(PD_neighbors);   % PD and PD-adjacent directions
    % Offset
    ff_cells_offset_alldir(:,icell) = ff_icell_offset_alldir;        % All 8 directions
    ff_cells_offset(:,icell) = ff_icell_offset_alldir(PD_neighbors); % PD and PD-adjacent directions
end
% Manual adjustment
if manualAdjust
    ff_cells_onset(2,25) = 0.8;     % Cell 25 PD
    ff_cells_onset(3,63) = 0.775;   % Cell 63 PD+1 
    ff_cells_onset(1,89) = 0.625;   % Cell 89 PD-1
end
% Convert onset times to receptive field radii
ff_cells_RF_alldir = movingField_radius - ff_cells_onset_alldir*barspeed; % All 8 directions
ff_cells_RF = movingField_radius - ff_cells_onset*barspeed;               % PD and PD-adjacent directions

% Build onset time and RF size distributions
% Initialize
ff_onset_distribution = NaN(length(binranges_onset),3); % Spike onset times
ff_RF_distribution = NaN(length(binranges_RF),3);       % Receptive field radii
for mdir = 1:3
    % Loop through PD and PD-adjacent directions
    ff_onset_distribution(:,mdir) = histc(ff_cells_onset(mdir,:),binranges_onset);
    ff_RF_distribution(:,mdir) = histc(ff_cells_RF(mdir,:),binranges_RF);
end
% Create histogram using PD and PD-adjacent directions
ff_onset_distribution_PD_PDadj = histc(ff_cells_onset(:),binranges_onset);

% Cubic smoothing spline
if smoothSpline
    ff_onset_dist_smoothed = csaps(binranges_onset,...
        ff_onset_distribution_PD_PDadj,smoothParam,edges);
    figure
    hold on
    bar(binranges_onset,ff_onset_distribution_PD_PDadj,'histc')
    %plot(edges,ff_onset_dist_smoothed,'r','LineWidth',1.2)
    xlim([startTime endTime])
    ylim([0 20])
    xlabel('Spike response onset time (sec)')
    ylabel('Cell count')
    title('Full-field PD and PD-adjacent directions')
end

% Fit beta distribution to histogram
if fitDist
    % Use PD and PD-adjacent directions
    ff_betaFit_onsetDist_PD_PDadj = fitdist(ff_cells_onset(:),'Beta')
    figure
    ff_histfit = histfit(ff_cells_onset(:),34,'beta');
    xlabel('Spike response onset time (sec)')
    xlim([startTime endTime])
    ylabel('Cell count')
    title('Full-field PD and PD-adjacent directions')
    legend(ff_histfit(2),'Beta')
end

% Fit Gaussian
% Use only PD
if fitGaussian
    ff_gaussianFit_onsetDist_PDonly = fitdist(ff_cells_onset(2,:)','normal')
    figure
    histfit(ff_cells_onset(2,:),34)
    xlabel('Spike response onset time (sec)')
    ylabel('Cell count')
    title('Full-field Preferred direction (PD) only')
    % Combine PD and PD-adjacent directions
    ff_gaussianFit_onsetDist_PD_PDadj = fitdist(ff_cells_onset(:),'normal')
    figure
    histfit(ff_cells_onset(:),34)
    xlabel('Spike response onset time (sec)')
    ylabel('Cell count')
    title('Full-field PD and PD-adjacent directions')
end

% Spike response duration
ff_resp_duration_alldir = ff_cells_offset_alldir - ff_cells_onset_alldir;
ff_resp_duration = ff_cells_offset - ff_cells_onset;
% Correlation between spike response onset and duration

% PD only
R = corrcoef(ff_cells_onset(2,ff_cells_idx),ff_resp_duration(2,ff_cells_idx));
rho = R(2,1)   % Pearson correlation coefficient
onset_duration_corr = polyfit(ff_cells_onset(2,ff_cells_idx),...  
    ff_resp_duration(2,ff_cells_idx),1);            % Find parameters of linear fit
if save_ff
    save('onset_duration_corr','onset_duration_corr');  % Save parameters of linear fit
end

% PD and PD-adjacent
R = corrcoef(ff_cells_onset(:,ff_cells_idx),ff_resp_duration(:,ff_cells_idx));
rho = R(2,1);   % Pearson correlation coefficient
ff_onset_duration_corr = polyfit(ff_cells_onset(:,ff_cells_idx),...
    ff_resp_duration(:,ff_cells_idx),1); % Parameters of linear fit
if save_ff
    save('ff_onset_duration_corr','ff_onset_duration_corr') % Save parameters
    save('ff_cells_onset','ff_cells_onset')                 % Save onset times
    save('ff_cells_offset','ff_cells_offset')               % Save offset times
end

% Visualization
if plot_ff
    figure  % Plot PD spike onset time distribution
    bar(binranges_onset,ff_onset_distribution(:,2),'histc')
    title('Spike Onset Time Distribution (PD)')
    xlabel('Time (sec)')
    ylabel('Cell Count')

    figure  % Plot PD receptive field size distribution
    bar(binranges_RF,ff_RF_distribution(:,2),'histc')
    title('Receptive Field Radius Distribution (PD)')
    xlabel('Distance from soma (pixels)')
    ylabel('Cell Count')
    
    figure  % Plot correlation between PD onset time and response duration
    hold on
    % PD-1 (left neighbor)
    scatter(ff_cells_onset(1,ff_cells_idx),ff_resp_duration(1,ff_cells_idx),...
        'filled','MarkerFaceColor','k')
    % PD+1 (right neighbor)
    ff_PDadj = scatter(ff_cells_onset(3,ff_cells_idx),ff_resp_duration(3,ff_cells_idx),...
        'filled','MarkerFaceColor','k');
    % PD
    ff_PD = scatter(ff_cells_onset(2,ff_cells_idx),ff_resp_duration(2,ff_cells_idx),...
        'filled','MarkerFaceColor',[0 0.4470 0.7410]);
    xaxis = linspace(startTime,endTime);            % Create x-values of linear fit
    linFit = polyval(ff_onset_duration_corr,xaxis); % Create y-values of linear fit
    linFit_PD = polyval(onset_duration_corr,xaxis); % Create y-values of linear fit (PD only)
    plotLinFit = plot(xaxis,linFit,'k');            % Plot linear fit
    plotLinFit_PD = plot(xaxis,linFit_PD,...
        'Color',[0 0.4470 0.7410]);                 % Plot linear fit (PD only)
    legend([ff_PDadj,ff_PD,plotLinFit],'PD-adjacent','PD','Fit')
    % Display Pearson correlation coefficient
    text(0.75,1.2,strcat('\rho =',{' '},num2str(rho)),'FontSize',11)
    xlabel('Onset time (sec)')
    ylabel('Duration (sec)')
    xlim([0 1.2])
    ylim([0 1.6])
    title('Full-field spike response onset-duration correlation')
end

%% Occlusion PD
% Direction index of each cell's PD
dir_idx_PD = all_cells_prefdir(:,2)/45+1; % Convert directions in degrees to indices 

% Initialize
% All 8 directions
occ_cells_onset_alldir = NaN(num_dir,Ncells);   % Onset
occ_cells_offset_alldir = NaN(num_dir,Ncells);  % Offset
% PD and PD-adjacent directions
occ_cells_onset = NaN(3,Ncells);                % Onset          
occ_cells_offset = NaN(3,Ncells);               % Offset 
for icell = occ_cells_idx
    %---------------------------------------------------------------------%
    % PD and adjacent directions of i-th cell
    icell_dir_idx_PD = dir_idx_PD(icell);   % Direction index of i-th cell's PD
    PD_leftNeighbor = icell_dir_idx_PD-1;   % Neighbor left to PD
    PD_rightNeighbor =  icell_dir_idx_PD+1; % Neighbor right to PD
    if PD_leftNeighbor == 0
        % If PD = 0 degrees, set left neighbor to be 315 degrees
        PD_leftNeighbor = num_dir;
    end
    if PD_rightNeighbor > num_dir
        % If PD = 315 degrees, set right neighbor to be 0 degrees
        PD_rightNeighbor = 1;
    end
    % Array of direction indices for PD and adjacent directions
    PD_neighbors = [PD_leftNeighbor icell_dir_idx_PD PD_rightNeighbor];
    
    %---------------------------------------------------------------------%
    % Load data
    fn = strcat('occlusion_cell',num2str(icell)); % Filename
    filepath = strcat('C:\Users\alber\Dropbox\Occlusion Project\sampling_distributions\',fn);
    load(filepath) % Load data
       
    %---------------------------------------------------------------------%
    data = SortedDSData.total;          % Use 'total' tab for threshold
    ktrials = size(data.spikeCounts,2); %repetitions
    spikeTimes = data.spikeTimes;       % Spike times for all directions
    
    % Get baseline firing rate and threshold
    noStimFR = NaN(ktrials,num_dir);       % Initialize    
    for mdir = 1:num_dir
        % Loop through all directions
        dirSpikeTimes = spikeTimes{1,mdir}; % Get spike times for one direction
        for ktrial = 1:ktrials
            % Loop through all trials
            % Spike times for k-th trial
            ktrialSpikeTimes = dirSpikeTimes{ktrial,1};
            % Firing rate when NO stimulus is present
            noStimFR(ktrial,mdir) = max(histc(ktrialSpikeTimes,[noStimStartTime noStimEndTime]));
        end
    end
    % Threshold for change-point detection
    baseline = mean(noStimFR,'all'); % Baseline firing rate 
                                     % Averaged across all trials and all directions
    threshold = threshold_param*std(noStimFR,0,'all') + baseline;
    
    %---------------------------------------------------------------------%
    % Find onset times
    data = SortedDSData.On;             % Use 'On' tab for occlusion PD
    ktrials = size(data.spikeCounts,2); %repetitions
    spikeTimes = data.spikeTimes;       % Spike times for all directions
    
    occ_icell_onset_alldir = NaN(num_dir,1);    % Initialize onset array
    occ_icell_offset_alldir = NaN(num_dir,1);   % Initialize offset array
    for mdir = 1:num_dir
        % Loop through all directions
        dirSpikeTimes = spikeTimes{1,mdir};     % Get spike times for one direction
        dirSpikeCount = NaN(ktrials,Nbins);     % Initialize
        for ktrial = 1:ktrials
            % Loop through all trials
            % Spike times for k-th trial
            ktrialSpikeTimes = dirSpikeTimes{ktrial,1};
            % Bin spike times for k-th trial
            dirSpikeCount(ktrial,:) = histc(ktrialSpikeTimes,edges);
            
            % Moving mean
            %dirSpikeCount(ktrial,:) = movmean(dirSpikeCount(ktrial,:),3);
            
        end
        % Trial-averaged firing rate across bins for one direction
        dirFR = mean(dirSpikeCount)/binSize;  %spikes/sec
        % Find time bin at which response starts or ends
        % Onset
        for ibin = 1:Nbins
            % Find time bin at which response starts
            if dirFR(ibin)>threshold && dirFR(ibin+1)>threshold
                occ_icell_onset_alldir(mdir) = edges(ibin);
                break % Exit for-loop once onset has been found
            end
        end
        % Offset
        for ibin = Nbins:-1:3                   % Loop backwards
            % Offset defined as above-baseline firing in 2 out of 3
            % consecutive time bins
            if dirFR(ibin)>threshold            % Detect above-baseline firing in first time bin
                if dirFR(ibin-1)>threshold      % Detect above-baseline firing in second time bin
                    occ_icell_offset_alldir(mdir) = edges(ibin);
                    break       % Exit for-loop once offset has been found
                else % If second time bin has zero response
                    if dirFR(ibin-2)>threshold  % Detect above-baseline firing in third time bin
                        occ_icell_offset_alldir(mdir) = edges(ibin);
                        break   % Exit for-loop once offset has been found
                    end
                end
            end
        end
    end
    % Store spike onset times
    occ_cells_onset_alldir(:,icell) = occ_icell_onset_alldir;        % All 8 directions
    occ_cells_onset(:,icell) = occ_icell_onset_alldir(PD_neighbors); % PD and PD-adjacent directions
    % Store spike offset times
    occ_cells_offset_alldir(:,icell) = occ_icell_offset_alldir;        % All 8 directions
    occ_cells_offset(:,icell) = occ_icell_offset_alldir(PD_neighbors); % PD and PD-adjacent directions
end
% Manual adjustment
if manualAdjust
    occ_cells_onset(2,12) = 0.425;      % Cell 12 PD
    occ_cells_onset(2,28) = 0.475;      % Cell 28 PD
    occ_cells_onset(1,30) = 0.45;       % Cell 30 PD-1
    occ_cells_onset(2,39) = 0.6;        % Cell 39 PD
end
% Convert onset times to receptive field radii
occ_cells_RF_alldir = movingField_radius - occ_cells_onset_alldir*barspeed; % All 8 directions
occ_cells_RF = movingField_radius - occ_cells_onset*barspeed;               % PD and PD-adjacent directions

% Build distributions
% Initialize
occ_onset_distribution = NaN(length(binranges_onset),3); % Spike onset times
occ_RF_distribution = NaN(length(binranges_RF),3);       % Receptive field radii
for mdir = 1:3
    % Loop through PD and PD-adjacent directions
    occ_onset_distribution(:,mdir) = histc(occ_cells_onset(mdir,:),binranges_onset);
    occ_RF_distribution(:,mdir) = histc(occ_cells_RF(mdir,:),binranges_RF);
end

% Create histogram using PD and PD-adjacent directions
occ_onset_distribution_PD_PDadj = histc(occ_cells_onset(:),binranges_onset);

% Cubic smoothing spline
if smoothSpline
    occ_onset_dist_smoothed = csaps(binranges_onset,...
        occ_onset_distribution_PD_PDadj,smoothParam,edges);
    figure
    hold on
    bar(binranges_onset,occ_onset_distribution_PD_PDadj,'histc')
    %plot(edges,occ_onset_dist_smoothed,'r','LineWidth',1.2)
    xlim([startTime endTime])
    ylim([0 20])
    xlabel('Spike response onset time (sec)')
    ylabel('Cell count')
    title('Occlusion PD and PD-adjacent directions')
end

% Fit beta distribution to histogram
if fitDist
    % Use PD and PD-adjacent directions
    occ_betaFit_onsetDist_PD_PDadj = fitdist(occ_cells_onset(:),'Weibull')
    figure
    occ_histfit = histfit(occ_cells_onset(:),27,'weibull');
    xlabel('Spike response onset time (sec)')
    xlim([startTime endTime])
    ylabel('Cell count')
    title('Occlusion PD and PD-adjacent directions')
    legend(occ_histfit(2),'Weibull')
end

% Fit Gaussian
% Use only PD
if fitGaussian
    occ_gaussianFit_onsetDist_PDonly = fitdist(occ_cells_onset(2,:)','normal')
    figure
    histfit(occ_cells_onset(2,:),27)
    xlabel('Spike response onset time (sec)')
    ylabel('Cell count')
    title('Occlusion Preferred direction (PD) only')
    % Combine PD and PD-adjacent directions
    occ_gaussianFit_onsetDist_PD_PDadj = fitdist(occ_cells_onset(:),'normal')
    figure
    histfit(occ_cells_onset(:),27)
    xlabel('Spike response onset time (sec)')
    ylabel('Cell count')
    title('Occlusion PD and PD-adjacent directions')
end

% Spike response duration
occ_resp_duration_alldir = occ_cells_offset_alldir - occ_cells_onset_alldir;
occ_resp_duration = occ_cells_offset - occ_cells_onset;
% Correlation between spike response onset and duration

% PD only
%{
R = corrcoef(occ_cells_onset(2,occ_cells_idx),occ_resp_duration(2,occ_cells_idx));
rho = R(2,1);   % Pearson correlation coefficient
occ_onset_duration_corr = polyfit(occ_cells_onset(2,occ_cells_idx),...  
    occ_resp_duration(2,occ_cells_idx),1);    % Find parameters of linear fit
if save_occ_pd
    save('occ_onset_duration_corr','occ_onset_duration_corr');  % Save parameters
end
%}

% PD and PD-adjacent
R = corrcoef(occ_cells_onset(~isnan(occ_cells_onset)),...
    occ_resp_duration(~isnan(occ_resp_duration)));
rho = R(2,1);   % Pearson correlation coefficient
occ_onset_duration_corr = polyfit(occ_cells_onset(~isnan(occ_cells_onset)),...
    occ_resp_duration(~isnan(occ_resp_duration)),1); % Parameters of linear fit
if save_occ_pd
    save('occ_onset_duration_corr','occ_onset_duration_corr')  % Save parameters
    save('occ_cells_onset','occ_cells_onset')                  % Save onset times
    save('occ_cells_offset','occ_cells_offset')                % Save offset times
end

% Visualization
if plot_occ_pd
    figure  % Plot spike onset time distribution
    bar(binranges_onset,occ_onset_distribution(:,2),'histc')
    title('Occlusion Spike Onset Time Distribution (PD)')
    xlabel('Time (sec)')
    ylabel('Cell Count')

    figure  % Plot receptive field size distribution
    bar(binranges_RF,occ_RF_distribution(:,2),'histc')
    title('Occlusion Receptive Field Radius Distribution (PD)')
    xlabel('Distance from soma (pixels)')
    ylabel('Cell Count')
    
    figure  % Plot correlation between PD onset time and response duration
    hold on
    % PD-1 (left neighbor)
    scatter(occ_cells_onset(1,occ_cells_idx),occ_resp_duration(1,occ_cells_idx),...
        'filled','MarkerFaceColor',[0.4660 0.6740 0.1880])
    % PD+1 (right neighbor)
    occ_PDadj = scatter(occ_cells_onset(3,occ_cells_idx),occ_resp_duration(3,occ_cells_idx),...
        'filled','MarkerFaceColor',[0.4660 0.6740 0.1880]);
    % PD
    occ_PD = scatter(occ_cells_onset(2,occ_cells_idx),occ_resp_duration(2,occ_cells_idx),'filled');
    xaxis = linspace(startTime,endTime/2);  % Create x-values of linear fit
    linFit = polyval(occ_onset_duration_corr,xaxis);    % Create y-values of linear fit
    plotLinFit = plot(xaxis,linFit,'k');    % Plot linear fit
    legend([occ_PDadj,occ_PD,plotLinFit],'PD-adjacent','PD','Fit')
    % Display Pearson correlation coefficient
    text(0.75,1,strcat('\rho =',{' '},num2str(rho)),'FontSize',11)
    xlabel('Onset time (sec)')
    ylabel('Duration (sec)')
    title('Spike response onset-duration correlation for occlusion')
end

%% Occlusion Null Response (ONR)
% Direction index of each cell's null direction
dir_idx_antiPD = mod(all_cells_prefdir(:,2)+180,360)/45+1; % Convert directions in degrees to indices

% Initialize
% All 8 directions
onr_cells_onset_alldir = NaN(num_dir,Ncells);   % Onset
onr_cells_offset_alldir = NaN(num_dir,Ncells);  % Offset
% anti-PD and anti-PD-adjacent directions
onr_cells_onset = NaN(3,Ncells);                % Onset
onr_cells_offset = NaN(3,Ncells);               % Offset
for icell = occ_cells_idx
    %---------------------------------------------------------------------%
    % PD and adjacent directions of i-th cell
    icell_dir_idx_antiPD = dir_idx_antiPD(icell);   % Direction index of i-th cell's PD
    antiPD_leftNeighbor = icell_dir_idx_antiPD-1;   % Neighbor left to PD
    antiPD_rightNeighbor =  icell_dir_idx_antiPD+1; % Neighbor right to PD
    if antiPD_leftNeighbor == 0
        % If PD = 0 degrees, set left neighbor to be 315 degrees
        antiPD_leftNeighbor = num_dir;
    end
    if antiPD_rightNeighbor > num_dir
        % If PD = 315 degrees, set right neighbor to be 0 degrees
        antiPD_rightNeighbor = 1;
    end
    % Array of direction indices for PD and adjacent directions
    antiPD_neighbors = [antiPD_leftNeighbor icell_dir_idx_antiPD antiPD_rightNeighbor];
    
    %---------------------------------------------------------------------%
    % Load data
    fn = strcat('occlusion_cell',num2str(icell)); % Filename
    filepath = strcat('C:\Users\alber\Dropbox\Occlusion Project\sampling_distributions\',fn);
    load(filepath) % Load data
    
    data = SortedDSData.Off;             % Use 'Off' tab for ONR
    ktrials = size(data.spikeCounts,2); %repetitions
    spikeTimes = data.spikeTimes;       % Spike times for all directions
    
    %---------------------------------------------------------------------%
    % Get baseline firing rate and threshold
    noStimFR = NaN(ktrials,num_dir);       % Initialize    
    for mdir = 1:num_dir
        % Loop through all directions
        dirSpikeTimes = spikeTimes{1,mdir}; % Get spike times for one direction
        for ktrial = 1:ktrials
            % Loop through all trials
            % Spike times for k-th trial
            ktrialSpikeTimes = dirSpikeTimes{ktrial,1};
            % Firing rate when NO stimulus is present
            noStimFR(ktrial,mdir) = max(histc(ktrialSpikeTimes,[noStimStartTime noStimEndTime]));
        end
    end
    % Threshold for change-point detection
    baseline = mean(noStimFR,'all'); % Baseline firing rate 
                                     % Averaged across all trials and all directions
    threshold = threshold_param*std(noStimFR,0,'all') + baseline;
    
    %---------------------------------------------------------------------%
    % Find onset and offset times
    onr_icell_onset_alldir = NaN(num_dir,1);        % Initialize onset array
    onr_icell_offset_alldir = NaN(num_dir,1);       % Initialize offset array
    for mdir = 1:num_dir
        % Loop through all directions
        dirSpikeTimes = spikeTimes{1,mdir};         % Get spike times for one direction
        dirSpikeCount = NaN(ktrials,length(edges)); % Initialize
        for ktrial = 1:ktrials
            % Loop through all trials
            % Spike times for k-th trial
            ktrialSpikeTimes = dirSpikeTimes{ktrial,1};
            % Bin spike times for k-th trial
            dirSpikeCount(ktrial,:) = histc(ktrialSpikeTimes,edges);
        end
        % Trial-averaged firing rate across bins for one direction
        % Onset
        dirFR = mean(dirSpikeCount)/binSize;  %spikes/sec
        for ibin = 1:Nbins
            % Find time bin at which response starts
            if dirFR(ibin)>threshold && dirFR(ibin+1)>threshold
                onr_icell_onset_alldir(mdir) = edges(ibin);
                break % Exit for-loop once onset has been found
            end
        end
        % Offset
        for ibin = Nbins:-1:3                   % Loop backwards
            % Offset defined as above-baseline firing in 2 out of 3
            % consecutive time bins
            if dirFR(ibin)>threshold            % Detect above-baseline firing in first time bin
                if dirFR(ibin-1)>threshold      % Detect above-baseline firing in second time bin
                    onr_icell_offset_alldir(mdir) = edges(ibin);
                    break       % Exit for-loop once offset has been found
                else % If second time bin has zero response
                    if dirFR(ibin-2)>threshold  % Detect above-baseline firing in third time bin
                        onr_icell_offset_alldir(mdir) = edges(ibin);
                        break   % Exit for-loop once offset has been found
                    end
                end
            end
        end
    end
    % Store spike onset times
    onr_cells_onset_alldir(:,icell) = onr_icell_onset_alldir;            % All 8 directions
    onr_cells_onset(:,icell) = onr_icell_onset_alldir(antiPD_neighbors); % Null and null-adjacent directions
    % Store spike offset times
    onr_cells_offset_alldir(:,icell) = onr_icell_offset_alldir;            % All 8 directions
    onr_cells_offset(:,icell) = onr_icell_offset_alldir(antiPD_neighbors); % Null and null-adjacent directions
end

% Convert onset times to receptive field radii
onr_cells_RF_alldir = onr_cells_onset_alldir*barspeed - movingField_radius; % All 8 directions
onr_cells_RF = onr_cells_onset*barspeed - movingField_radius;               % Null and null-adjacent directions

% Build distributions
% Initialize
onr_onset_distribution = NaN(length(binranges_onset),3); % Spike onset times
onr_RF_distribution = NaN(length(binranges_RF),3);       % Receptive field radii
for mdir = 1:3
    % Loop through null and null-adjacent directions
    onr_onset_distribution(:,mdir) = histc(onr_cells_onset(mdir,:),binranges_onset);
    onr_RF_distribution(:,mdir) = histc(onr_cells_RF(mdir,:),binranges_RF);
end

% Create histogram using PD and PD-adjacent directions
onr_onset_distribution_PD_PDadj = histc(onr_cells_onset(:),binranges_onset);

% Cubic smoothing spline
if smoothSpline
    onr_onset_dist_smoothed = csaps(binranges_onset,...
        onr_onset_distribution_PD_PDadj,smoothParam_onr,edges);
    figure
    hold on
    bar(binranges_onset,onr_onset_distribution_PD_PDadj,'histc')
    %plot(edges,onr_onset_dist_smoothed,'r','LineWidth',1.2)
    xlim([startTime endTime])
    ylim([0 50])
    xlabel('Spike response onset time (sec)')
    ylabel('Cell count')
    title('Occlusion null and null-adjacent directions')
end

% Fit gamma distribution to histogram
if fitDist
    % Use null and null-adjacent directions
    onr_gammaFit_onsetDist_PD_PDadj = fitdist(onr_cells_onset(:),'Gamma')
    figure
    onr_histfit = histfit(onr_cells_onset(:),15,'gamma');
    xlabel('Spike response onset time (sec)')
    xlim([startTime endTime])
    ylabel('Cell count')
    title('Occlusion null and null-adjacent directions')
    legend(onr_histfit(2),'Gamma')
end

% Fit Gaussian
% Use only PD
if fitGaussian
    onr_gaussianFit_onsetDist_PDonly = fitdist(onr_cells_onset(2,:)','normal')
    figure
    histfit(onr_cells_onset(2,:),15)
    xlabel('Spike response onset time (sec)')
    ylabel('Cell count')
    title('Occlusion null direction only')
    % Combine PD and PD-adjacent directions
    onr_gaussianFit_onsetDist_PD_PDadj = fitdist(onr_cells_onset(:),'normal')
    figure
    histfit(onr_cells_onset(:),15)
    xlabel('Spike response onset time (sec)')
    ylabel('Cell count')
    title('Occlusion null and null-adjacent directions')
end

% Spike response duration
onr_resp_duration_alldir = onr_cells_offset_alldir - onr_cells_onset_alldir;
onr_resp_duration = onr_cells_offset - onr_cells_onset;
% Correlation between spike response onset and duration

% ND only

R = corrcoef(onr_cells_onset(2,~isnan(onr_cells_onset(2,:))),...
    onr_resp_duration(2,~isnan(onr_cells_onset(2,:))));
rho = R(2,1)   % Pearson correlation coefficient
onr_onset_duration_corr_NDonly = polyfit(onr_cells_onset(2,~isnan(onr_cells_onset(2,:))),...  
    onr_resp_duration(2,~isnan(onr_cells_onset(2,:))),1);    % Find parameters of linear fit
if save_onr
    save('onr_onset_duration_corr_NDonly',...
        'onr_onset_duration_corr_NDonly');  % Save parameters
end

% PD and PD-adjacent
R = corrcoef(onr_cells_onset(~isnan(onr_cells_onset)),...
    onr_resp_duration(~isnan(onr_resp_duration)));
rho = R(2,1);   % Pearson correlation coefficient
onr_onset_duration_corr = polyfit(onr_cells_onset(~isnan(onr_cells_onset)),...
    onr_resp_duration(~isnan(onr_resp_duration)),1); % Parameters of linear fit
if save_onr
    save('onr_onset_duration_corr','onr_onset_duration_corr')  % Save parameters
    save('onr_cells_onset','onr_cells_onset')                  % Save onset times
    save('onr_cells_offset','onr_cells_offset')                % Save offset times
end

% Visualization
if plot_onr
    figure  % Plot spike onset time distribution
    bar(binranges_onset,onr_onset_distribution(:,2),'histc')
    title('Occlusion Spike Onset Time Distribution (Null)')
    xlabel('Time (sec)')
    ylabel('Cell Count')

    figure  % Plot receptive field size distribution
    bar(binranges_RF,onr_RF_distribution(:,2),'histc')
    title('Occlusion Receptive Field Radius Distribution (Null)')
    xlabel('Distance from soma (pixels)')
    ylabel('Cell Count')
    
    figure  % Plot correlation between null direction onset time and response duration
    hold on
    
    % PD-1 (left neighbor)
    scatter(onr_cells_onset(1,occ_cells_idx),onr_resp_duration(1,occ_cells_idx),...
        'filled','MarkerFaceColor','k')%[0.4660 0.6740 0.1880])
    % PD+1 (right neighbor)
    onr_NDadj = scatter(onr_cells_onset(3,occ_cells_idx),onr_resp_duration(3,occ_cells_idx),...
        'filled','MarkerFaceColor','k')%[0.4660 0.6740 0.1880]);
    % ND
    onr_ND = scatter(onr_cells_onset(2,occ_cells_idx),onr_resp_duration(2,occ_cells_idx),...
        'filled','MarkerFaceColor',[0 0.4470 0.7410]);
    xaxis = linspace(startTime/2,endTime);  % Create x-values of linear fit
    linFit = polyval(onr_onset_duration_corr,xaxis);            % Create y-values of linear fit
    linFit_ND = polyval(onr_onset_duration_corr_NDonly,xaxis);  % Create y-values of linear fit (ND only)
    plotLinFit = plot(xaxis,linFit,'k');                        % Plot linear fit
    plotLinFit_ND = plot(xaxis,linFit_ND,...
        'Color',[0 0.4470 0.7410]);                             % Plot linear fit (ND only)
    legend([onr_NDadj,onr_ND],'PD-adjacent','PD')
    % Display Pearson correlation coefficient
    text(1.65,0.4,strcat('\rho =',{' '},num2str(rho)),'FontSize',11)
    xlabel('Onset time (sec)')
    ylabel('Duration (sec)')
    xlim([1.2 1.8])
    ylim([-0.2 0.8])
    title('Spike response onset-duration correlation for occlusion null response')
end

