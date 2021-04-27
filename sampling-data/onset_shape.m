%{
(1) Plot traces of onset times for all cells with their preferred
    directions (PDs) aligned
(2) Same as above but align PD onset time
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
Ncells = max(max(ff_cells_idx),max(occ_cells_idx));

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

%{
Use data from bar speed 300 pixels/sec to get baseline firing rate
There is no firing (zero response) when no stimulus is present
%}
noStimStartTime = 5;        %sec
noStimEndTime = 6;          %sec

manualAdjust = 1;   % Manually adjust? (Usually yes)

% Plotting
theta_axis = -180:45:135;   % Horizontal (theta) axis
siglevel = 0.01;            % Significan level

% Save files?
save_ff_onset_aligned = 0;      % Save full-field onset times aligned by PD?
save_occ_onset_aligned = 0;     % Save occlusion PD onset times aligned by PD?
save_onr_onset_aligned = 0;     % Save occlusion null response onset times aligned by anti-PD?

%% Load Data
% Direction index of each cell's PD
cd('C:\Users\alber\Dropbox\Occlusion Project\sampling_distributions') % Specify directory
load('all_cells_prefdir')                                             % Load appropriate file
% First column for full-field
% Second column for occlusion

%% Full-field only
dir_idx_PD = all_cells_prefdir(:,1)/45+1; % Convert directions in degrees to indices 

% Initialize
ff_cells_onset_alldir = NaN(num_dir,Ncells);            % All 8 directions
ff_cells_onset_aligned = NaN(num_dir+num_dir/2,Ncells); % All 8 directions aligned
for icell = ff_cells_idx
    %---------------------------------------------------------------------%
    % PD of i-th cell
    icell_dir_idx_PD = dir_idx_PD(icell);   % Direction index of i-th cell's PD
        
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
    noStimFR = NaN(ktrials,num_dir);     % Initialize    
    for mdir = 1:num_dir
        % Loop through all directions
        dirSpikeTimes = spikeTimes{1,mdir};    % Get spike times for one direction
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
    threshold = threshold_param*std(noStimFR,0,'all') + baseline; % Adjust threshold_param
    
    %---------------------------------------------------------------------%
    % Find onset times
    ff_icell_onset_alldir = NaN(num_dir,1);   % Initialize
    for mdir = 1:num_dir
        % Loop through all directions
        dirSpikeTimes = spikeTimes{1,mdir};   % Get spike times for one direction
        dirSpikeCount = NaN(ktrials,Nbins);   % Initialize
        for ktrial = 1:ktrials
            % Loop through all trials
            % Spike times for k-th trial
            ktrialSpikeTimes = dirSpikeTimes{ktrial,1};
            % Bin spikes for k-th trial
            dirSpikeCount(ktrial,:) = histc(ktrialSpikeTimes,edges);
        end
        % Trial-averaged firing rate across bins for one direction
        dirFR = mean(dirSpikeCount)/binSize;  %spikes/sec
        for ibin = 1:Nbins
            % Find time bin at which response starts
            if dirFR(ibin)>threshold && dirFR(ibin+1)>threshold
                ff_icell_onset_alldir(mdir) = edges(ibin);
                break % Exit for-loop once onset has been found
            end
        end
    end
    % Store spike onset times
    ff_cells_onset_alldir(:,icell) = ff_icell_onset_alldir; % All 8 directions    
    % Align PD
    ff_cells_onset_aligned(num_dir/2+1:end,icell) = ...     
        [ff_icell_onset_alldir(icell_dir_idx_PD:end); ...   % Right of PD
        ff_icell_onset_alldir(1:icell_dir_idx_PD-1)];       % Left of PD
end
ff_cells_onset_aligned(1:num_dir/2,:) = ff_cells_onset_aligned(num_dir+1:end,:);
ff_cells_onset_aligned = ff_cells_onset_aligned(1:num_dir,:);
% Manual adjustment
if manualAdjust
    ff_cells_onset_aligned(num_dir/2+1,25) = 0.8;     % Cell 25 PD
    ff_cells_onset_aligned(num_dir/2+2,63) = 0.775;   % Cell 63 PD+1 
    ff_cells_onset_aligned(num_dir/2,89) = 0.625;     % Cell 89 PD-1
end
ff_cells_meanOnset = nanmean(ff_cells_onset_aligned,2);     % Average across cells

% (2) Align PD onset time
ff_cells_onset_aligned_timeFixed = ff_cells_onset_aligned - ...
    repmat(ff_cells_onset_aligned(num_dir/2+1,:),num_dir,1);                % Subtract PD onset time from all entries
ff_cells_meanOnset_timeFixed = nanmean(ff_cells_onset_aligned_timeFixed,2); % Average across cells

% (1) Compute standard deviation and error bars (SEM)
ff_cells_onset_std = nanstd(ff_cells_onset_aligned,0,2);    % Sigma of onset time for each direction
ff_cells_onset_sem = ff_cells_onset_std/sqrt(ff_Ncells);    % SEM for each direction

% (2) Determine whether PD-adjacent delays are significant via Student's t-test
[significant(1) p_value(1)] = ttest(ff_cells_onset_aligned_timeFixed(num_dir/2,:));   % PD left-adjacent
[significant(2) p_value(2)] = ttest(ff_cells_onset_aligned_timeFixed(num_dir/2+2,:)); % PD right-adjacent
% (2) Compute error bars (SEM)
ff_cells_onset_timeFixed_sem = nanstd(ff_cells_onset_aligned_timeFixed,0,2)/sqrt(ff_Ncells);

% Visualization
figure  % (1) Plot traces of onset times for all cells with their preferred
        %     directions (PDs) aligned along with error bars (SEM)
hold on % Plot multiple traces of individual cells
for icell = ff_cells_idx
    plot(theta_axis,ff_cells_onset_aligned(:,icell))
    %plot(theta_axis,ff_cells_onset_aligned(:,icell),...
    %    'Color',[0.3010 0.7450 0.9330])
end
errorbar(theta_axis,ff_cells_meanOnset,ff_cells_onset_sem,...
    'Color',[0 0.4470 0.7410],'LineWidth',1.6)
text(-210,0,'\sigma =','FontSize',12)
for mdir = 1:num_dir
    % Display standard deviation (sigma) for each direction
    text(theta_axis(mdir),0,num2str(round(ff_cells_onset_std(mdir),2)))
end
xlim([-225 180])
xlabel('Deviation from PD')
ylim([-0.1 1.8])
ylabel('Spike Onset Time (sec)')
title('Full-field only')

figure  % (2) Same as above but align PD onset time
hold on % Plot multiple traces of individual cells
for icell = ff_cells_idx
    plot(theta_axis,ff_cells_onset_aligned_timeFixed(:,icell))
    %plot(theta_axis,ff_cells_onset_aligned_timeFixed(:,icell),...
    %    'Color',[0.3010 0.7450 0.9330])
end
errorbar(theta_axis,ff_cells_meanOnset_timeFixed,ff_cells_onset_timeFixed_sem,...
    'Color',[0 0.4470 0.7410],'LineWidth',1.8)
if p_value(1) < siglevel   % PD left-adjacent
    % If significant, draw significance asterisk
    text(theta_axis(num_dir/2),0.5,'*','FontSize',16)
end
if p_value(2) < siglevel   % PD right-adjacent
    % If significant, draw significance asterisk
    text(theta_axis(num_dir/2+2),0.5,'*','FontSize',16)
end
text(theta_axis(end)-5,1,'*p < 0.01','FontSize',10)
xlim([-225 180])
xticks([-225:45:180])
xlabel('Deviation from PD')
%ylim([0 1])
ylabel('Spike Onset Time Relative to PD Spike Onset Time (sec)')
title('Full-field only')

if save_ff_onset_aligned
    save('ff_cells_onset_aligned','ff_cells_onset_aligned')
end

%% Occlusion only
dir_idx_PD = all_cells_prefdir(:,2)/45+1; % Convert directions in degrees to indices 

% Initialize
occ_cells_onset_alldir = NaN(num_dir,Ncells);            % All 8 directions
occ_cells_onset_aligned = NaN(num_dir+num_dir/2,Ncells); % All 8 directions aligned
for icell = occ_cells_idx
    %---------------------------------------------------------------------%
    % PD of i-th cell
    icell_dir_idx_PD = dir_idx_PD(icell);   % Direction index of i-th cell's PD
        
    %---------------------------------------------------------------------%
    % Load data
    fn = strcat('occlusion_cell',num2str(icell)); % Filename
    filepath = strcat('C:\Users\alber\Dropbox\Occlusion Project\sampling_distributions\',fn);
    load(filepath) % Load data
    
    %---------------------------------------------------------------------%
    % Get baseline firing rate and threshold
    data = SortedDSData.total;               % Get data structure
                                             % Use 'total' tab for threshold !!!
    ktrials = size(data.spikeCounts,2);      %repetitions
    spikeTimes = data.spikeTimes;            % Spike times for all directions
    

    noStimFR = NaN(ktrials,num_dir);     % Initialize    
    for mdir = 1:num_dir
        % Loop through all directions
        dirSpikeTimes = spikeTimes{1,mdir};    % Get spike times for one direction
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
    threshold = threshold_param*std(noStimFR,0,'all') + baseline; % Adjust threshold_param

    %---------------------------------------------------------------------%
    % Find onset times
    data = SortedDSData.On;                  % Get data structure
                                             % Use 'On' tab for occlusion PD response !!!
    ktrials = size(data.spikeCounts,2);      %repetitions
    spikeTimes = data.spikeTimes;            % Spike times for all directions
    
    occ_icell_onset_alldir = NaN(num_dir,1);   % Initialize
    for mdir = 1:num_dir
        % Loop through all directions
        dirSpikeTimes = spikeTimes{1,mdir};   % Get spike times for one direction
        dirSpikeCount = NaN(ktrials,Nbins);   % Initialize
        for ktrial = 1:ktrials
            % Loop through all trials
            % Spike times for k-th trial
            ktrialSpikeTimes = dirSpikeTimes{ktrial,1};
            % Bin spikes for k-th trial
            dirSpikeCount(ktrial,:) = histc(ktrialSpikeTimes,edges);
        end
        % Trial-averaged firing rate across bins for one direction
        dirFR = mean(dirSpikeCount)/binSize;  %spikes/sec
        for ibin = 1:Nbins
            % Find time bin at which response starts
            if dirFR(ibin)>threshold && dirFR(ibin+1)>threshold
                occ_icell_onset_alldir(mdir) = edges(ibin);
                break % Exit for-loop once onset has been found
            end
        end
    end
    % Store spike onset times
    occ_cells_onset_alldir(:,icell) = occ_icell_onset_alldir; % All 8 directions
    % Align PD
    occ_cells_onset_aligned(num_dir/2+1:end,icell) = ...     
        [occ_icell_onset_alldir(icell_dir_idx_PD:end); ...   % Right of PD
        occ_icell_onset_alldir(1:icell_dir_idx_PD-1)];       % Left of PD
end
occ_cells_onset_aligned(1:num_dir/2,:) = occ_cells_onset_aligned(num_dir+1:end,:);
occ_cells_onset_aligned = occ_cells_onset_aligned(1:num_dir,:);
% Manual adjustment
if manualAdjust
    occ_cells_onset_aligned(num_dir/2+1,12) = 0.425;      % Cell 12 PD
    occ_cells_onset_aligned(num_dir/2+1,28) = 0.475;      % Cell 28 PD
    occ_cells_onset_aligned(num_dir/2,30) = 0.45;         % Cell 30 PD-1
    occ_cells_onset_aligned(num_dir/2+1,39) = 0.6;        % Cell 39 PD
end
occ_cells_meanOnset = nanmean(occ_cells_onset_aligned,2);    % Average across cells

% (2) Align PD onset time
occ_cells_onset_aligned_timeFixed = occ_cells_onset_aligned - ...
    repmat(occ_cells_onset_aligned(num_dir/2+1,:),num_dir,1);                 % Subtract PD onset time from all entries
occ_cells_meanOnset_timeFixed = nanmean(occ_cells_onset_aligned_timeFixed,2); % Average across cells

% (1) Compute standard deviation and error bar (SEM)
occ_cells_onset_std = nanstd(occ_cells_onset_aligned,0,2);    % Sigma of onset time for each direction
occ_cells_onset_sem = occ_cells_onset_std/sqrt(occ_Ncells);   % SEM for each direction

% (2) Determine whether PD-adjacent delays are significant via Student's t-test
[significant(1) p_value(1)] = ttest(occ_cells_onset_aligned_timeFixed(num_dir/2,:));   % PD left-adjacent
[significant(2) p_value(2)] = ttest(occ_cells_onset_aligned_timeFixed(num_dir/2+2,:)); % PD right-adjacent
% (2) Compute error bars (SEM)
occ_cells_onset_timeFixed_sem = nanstd(occ_cells_onset_aligned_timeFixed,0,2)/sqrt(occ_Ncells);

% Visualization
figure
hold on % Plot multiple traces of individual cells
for icell = occ_cells_idx %ff_cells_idx
    plot(theta_axis,occ_cells_onset_aligned(:,icell))
    %plot(theta_axis,occ_cells_onset_aligned(:,icell),...
    %    'Color',[0.3010 0.7450 0.9330])
end
errorbar(theta_axis,occ_cells_meanOnset,occ_cells_onset_sem,...
    'Color',[0 0.4470 0.7410],'LineWidth',1.6)
text(-210,0,'\sigma =','FontSize',12)
for mdir = 1:num_dir
    % Display standard deviation (sigma) for each direction
    text(theta_axis(mdir),0,num2str(round(occ_cells_onset_std(mdir),2)))
end
xlim([-225 180])
xlabel('Deviation from PD')
ylim([-0.1 1])
ylabel('Spike Onset Time (sec)')
title('Entering occlusion only')

figure  % (2) Same as above but align PD onset time
hold on % Plot multiple traces of individual cells
for icell = occ_cells_idx %ff_cells_idx
    plot(theta_axis,occ_cells_onset_aligned_timeFixed(:,icell))
    %plot(theta_axis,occ_cells_onset_aligned_timeFixed(:,icell),...
    %    'Color',[0.3010 0.7450 0.9330])
end
errorbar(theta_axis,occ_cells_meanOnset_timeFixed,occ_cells_onset_timeFixed_sem,...
    'Color',[0 0.4470 0.7410],'LineWidth',1.8)
if p_value(1) < siglevel   % PD left-adjacent
    % If significant, draw significance asterisk 
    text(theta_axis(num_dir/2),0.5,'*','FontSize',16)
end
if p_value(2) < siglevel   % PD right-adjacent
    % If significant, draw significance asterisk
    text(theta_axis(num_dir/2+2),0.5,'*','FontSize',16)
end
text(theta_axis(end)-5,0.9,'*p < 0.01','FontSize',10) 
xlim([-225 180])
xlabel('Deviation from PD')
%ylim([0 1])
ylabel('Spike Onset Time Relative to PD Spike Onset Time (sec)')
title('Occlusion only')

if save_occ_onset_aligned
    save('occ_cells_onset_aligned','occ_cells_onset_aligned')
end

%% Occlusion null response    
all_cells_nulldir = mod(all_cells_prefdir+180,360); % Find null directions
dir_idx_antiPD = all_cells_nulldir(:,2)/45+1;       % Convert directions in degrees to indices 

% Initialize
onr_cells_onset_alldir = NaN(num_dir,Ncells);            % All 8 directions
onr_cells_onset_aligned = NaN(num_dir+num_dir/2,Ncells); % All 8 directions aligned
for icell = occ_cells_idx
    %---------------------------------------------------------------------%
    % PD of i-th cell
    icell_dir_idx_antiPD = dir_idx_antiPD(icell);   % Direction index of i-th cell's PD
        
    %---------------------------------------------------------------------%
    % Load data
    fn = strcat('occlusion_cell',num2str(icell)); % Filename
    filepath = strcat('C:\Users\alber\Dropbox\Occlusion Project\sampling_distributions\',fn);
    load(filepath) % Load data
    
    data = SortedDSData.Off;                 % Get data structure
    ktrials = size(data.spikeCounts,2);      %repetitions
    spikeTimes = data.spikeTimes;            % Spike times for all directions
    
    %---------------------------------------------------------------------%
    % Get baseline firing rate and threshold
    noStimFR = NaN(ktrials,num_dir);         % Initialize    
    for mdir = 1:num_dir
        % Loop through all directions
        dirSpikeTimes = spikeTimes{1,mdir};  % Get spike times for one direction
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
    threshold = threshold_param*std(noStimFR,0,'all') + baseline; % Adjust threshold_param
    
    %---------------------------------------------------------------------%
    % Find onset times
    onr_icell_onset_alldir = NaN(num_dir,1); % Initialize
    for mdir = 1:num_dir
        % Loop through all directions
        dirSpikeTimes = spikeTimes{1,mdir};  % Get spike times for one direction
        dirSpikeCount = NaN(ktrials,Nbins);  % Initialize
        for ktrial = 1:ktrials
            % Loop through all trials
            % Spike times for k-th trial
            ktrialSpikeTimes = dirSpikeTimes{ktrial,1};
            % Bin spikes for k-th trial
            dirSpikeCount(ktrial,:) = histc(ktrialSpikeTimes,edges);
        end
        % Trial-averaged firing rate across bins for one direction
        dirFR = mean(dirSpikeCount)/binSize;  %spikes/sec
        for ibin = 1:Nbins
            % Find time bin at which response starts
            if dirFR(ibin)>threshold && dirFR(ibin+1)>threshold
                onr_icell_onset_alldir(mdir) = edges(ibin);
                break % Exit for-loop once onset has been found
            end
        end
    end
    % Store spike onset times
    onr_cells_onset_alldir(:,icell) = onr_icell_onset_alldir; % All 8 directions
    % Align anti-PD
    onr_cells_onset_aligned(num_dir/2+1:end,icell) = ...     
        [onr_icell_onset_alldir(icell_dir_idx_antiPD:end); ...   % Right of PD
        onr_icell_onset_alldir(1:icell_dir_idx_antiPD-1)];       % Left of PD
end
onr_cells_onset_aligned(1:num_dir/2,:) = onr_cells_onset_aligned(num_dir+1:end,:);
onr_cells_onset_aligned = onr_cells_onset_aligned(1:num_dir,:);
onr_cells_meanOnset = nanmean(onr_cells_onset_aligned,2);    % Average across cells

% (2) Align PD onset time
onr_cells_onset_aligned_timeFixed = onr_cells_onset_aligned - ...
    repmat(onr_cells_onset_aligned(num_dir/2+1,:),num_dir,1);                 % Subtract PD onset time from all entries
onr_cells_meanOnset_timeFixed = nanmean(onr_cells_onset_aligned_timeFixed,2); % Average across cells

% (1) Compute standard deviation and error bar (SEM)
onr_cells_onset_std = nanstd(onr_cells_onset_aligned,0,2);    % Sigma of onset time for each direction
onr_cells_onset_sem = onr_cells_onset_std/sqrt(occ_Ncells);   % SEM for each direction

% (2) Determine whether PD-adjacent delays are significant via Student's t-test
[significant(1) p_value(1)] = ttest(onr_cells_onset_aligned_timeFixed(num_dir/2,:));   % PD left-adjacent
[significant(2) p_value(2)] = ttest(onr_cells_onset_aligned_timeFixed(num_dir/2+2,:)); % PD right-adjacent
% (2) Compute error bars (SEM)
onr_cells_onset_timeFixed_sem = nanstd(onr_cells_onset_aligned_timeFixed,0,2)/sqrt(occ_Ncells);

% Visualization
figure
hold on % Plot multiple traces of individual cells
for icell = occ_cells_idx %ff_cells_idx
    plot(theta_axis,onr_cells_onset_aligned(:,icell))
    %plot(theta_axis,onr_cells_onset_aligned(:,icell),...
    %    'Color',[0.3010 0.7450 0.9330])
end
errorbar(theta_axis,onr_cells_meanOnset,onr_cells_onset_sem,...
    'Color',[0 0.4470 0.7410],'LineWidth',1.6)
text(-210,1.1,'\sigma =','FontSize',12)
for mdir = 1:num_dir
    % Display standard deviation (sigma) for each direction
    text(theta_axis(mdir),1.1,num2str(round(onr_cells_onset_std(mdir),2)))
end
xlim([-225 180])
xlabel('Deviation from null direction')
ylim([1 2])
ylabel('Spike Onset Time (sec)')
title('Occlusion null response')

figure  % (2) Same as above but align PD onset time
hold on % Plot multiple traces of individual cells
for icell = occ_cells_idx %ff_cells_idx
    plot(theta_axis,onr_cells_onset_aligned_timeFixed(:,icell))
    %plot(theta_axis,onr_cells_onset_aligned_timeFixed(:,icell),...
    %    'Color',[0.3010 0.7450 0.9330])
end
errorbar(theta_axis,onr_cells_meanOnset_timeFixed,onr_cells_onset_timeFixed_sem,...
    'Color',[0 0.4470 0.7410],'LineWidth',1.8)
if p_value(1) < siglevel   % PD left-adjacent
    % If significant, draw significance asterisk 
    text(theta_axis(num_dir/2),1.5,'*','FontSize',16)
end
if p_value(2) < siglevel   % PD right-adjacent
    % If significant, draw significance asterisk
    text(theta_axis(num_dir/2+2),1.5,'*','FontSize',16)
end
text(theta_axis(end)-5,0.6,'*p < 0.01','FontSize',10) 
xlim([-225 180])
xticks([-225:45:180])
xlabel('Deviation from null direction')
ylim([-0.3 0.7])
ylabel('Spike Onset Time Relative to anti-PD Spike Onset Time (sec)')
title('Occlusion null response')

if save_onr_onset_aligned
    save('onr_cells_onset_aligned','onr_cells_onset_aligned')
end


