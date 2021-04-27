%% Parameters
% List of full-field cell indices
ff_cells_idx = [1:25 32:47 49:52 54:60 63:64 66:79 82:83 87:90];   
ff_Ncells = length(ff_cells_idx);                   % Number of full-field cells
% List of occlusion cell indices
occ_cells_idx = [1:16 18 21 23 26:48 53 55 58:64 66 68 71:73 76 79:91];
occ_Ncells = length(occ_cells_idx);                 % Number of occlusion cells
Ncells = max(max(ff_cells_idx),max(occ_cells_idx)); % Total number of cells

num_dir = 8;           % Number of directions

binSize = 0.025;       %sec
startTime = 0;         %sec
endTime = 2;           %sec
edges = startTime:binSize:endTime; % Create time bins
Nbins = length(edges); % Number of bins

histBin = 12;          % Bin size of histogram

% Moving mean
win = 0.100;           %sec

% Fitting
fitDist = 0;           % Fit distribution to histogram?

save_peakFR = 0;       % Save peak firing rates?

% Plotting
plot_ff = 0;           % Plot full-field?
plot_occ = 0;          % Plot occlusion PD?
plot_onr = 0;          % Plot occlusion null?

% Onset-peak correlation
onset_peak_corr = 1;   % Analyze onset-peak correlation?

% Load direction indices
cd('C:\Users\alber\Dropbox\Occlusion Project\sampling_distributions') % Specify directory
load('all_cells_prefdir')                                             % Load appropriate file
% 1st column is for full-field
% 2nd column is for occlusion


%% Full-field
% Direction index of each cell's PD
dir_idx_PD = all_cells_prefdir(:,1)/45+1; % Convert directions in degrees to indices 

ff_cells_peakFR = NaN(3,Ncells);            % Initialize
for icell = ff_cells_idx
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
    
    % Load data
    fn = strcat('fullfield_cell',num2str(icell)); % Filename
    filepath = strcat('C:\Users\alber\Dropbox\Occlusion Project\sampling_distributions\',fn);
    load(filepath) % Load data
    
    data = SortedDSData.total;                % Get data structure
    ktrials = size(data.spikeCounts,2);       %repetitions
    spikeTimes = data.spikeTimes;             % Spike times for all directions
    
    heatmap = NaN(num_dir,length(edges));     % Initialize heat map
    for mdir = 1:num_dir
        % Loop through all directions
        dirSpikeTimes = spikeTimes{1,mdir};   % Get spike times for one direction
        dirSpikeCount = NaN(ktrials,length(edges)); % Initialize
        for ktrial = 1:ktrials
            % Loop through all trials
            % Bin spike times for k-th trial
            ktrialSpikeTimes = dirSpikeTimes{ktrial,1};
            dirSpikeCount(ktrial,:) = histc(ktrialSpikeTimes,edges);
            
            % Moving mean
            %temp = histc(ktrialSpikeTimes,edges);
            %dirSpikeCount(ktrial,:) = movmean((temp),win/binSize);
            
        end
        heatmap(mdir,:) = mean(dirSpikeCount);      % Average across trials
    end
    heatmap = heatmap/binSize;                % Convert spike count to firing rate
    icell_peakFR_alldir = max(heatmap,[],2);  % Peak firing rate of i-th cell in all directions
    % Peak firing rate of i-th cell in PD and PD-adjacent directions
    icell_peakFR = icell_peakFR_alldir(PD_neighbors);
    ff_cells_peakFR(:,icell) = icell_peakFR;  % Store peak firing rate
end
ff_peakFR_PD = ff_cells_peakFR(2,:);                % Get PD peak firing rates
ff_peakFR_PD = ff_peakFR_PD(~isnan(ff_peakFR_PD));  % Remove NaN values
% Save PD peak firing rates
if save_peakFR
    save('ff_peakFR_PD','ff_peakFR_PD')
end

binranges = 0:histBin:max(ff_cells_peakFR,[],'all')+5;
ff_peakFR_distribution = NaN(length(binranges),3); % Initialize
for mdir = 1:3
    % Loop through PD and PD-adjacent directions
    ff_peakFR_distribution(:,mdir) = histc(ff_cells_peakFR(mdir,:),binranges);
end

% Fit gamma distribution to histogram
if fitDist
    ff_gammaFit_peakFRdist_PDonly = fitdist(ff_cells_peakFR(2,:)','gamma')
    figure
    histfit(ff_cells_peakFR(2,:),12,'Gamma')
    xlabel('Peak firing rate (Hz)')
    ylabel('Cell count')
    title('Full-field Preferred direction (PD) only')
end

% Combine PD left- and right-adjacent neighbors
ff_peakFR_neighbors = [ff_cells_peakFR(1,:) ff_cells_peakFR(end,:)];
ff_peakFR_distribution_neighbors = histc(ff_peakFR_neighbors,binranges);

if plot_ff
    figure % PD-1
    bar(binranges,ff_peakFR_distribution(:,1),'histc')
    title(strcat('Full-field Preferred Direction Left Adjacent',...
        ' (bin size =',{' '},num2str(binSize*1000),' ms)'))
    xlabel('Peak Firing Rate (Hz)')
    ylabel('Cell Count')
    figure % PD
    bar(binranges,ff_peakFR_distribution(:,2),'histc')
    title(strcat('Full-field Preferred Direction',...
        ' (bin size =',{' '},num2str(binSize*1000),' ms)'))
    xlabel('Peak Firing Rate (Hz)')
    ylabel('Cell Count')
    figure % PD+1
    bar(binranges,ff_peakFR_distribution(:,3),'histc')
    title(strcat('Full-field Preferred Direction Right Adjacent',...
        ' (bin size =',{' '},num2str(binSize*1000),' ms)'))
    xlabel('Peak Firing Rate (Hz)')
    ylabel('Cell Count')
    figure % PD neighbors
    bar(binranges,ff_peakFR_distribution_neighbors,'histc')
    title(strcat('Full-field Preferred Direction Adjacent Neighbors',...
        ' (bin size =',{' '},num2str(binSize*1000),' ms)'))
    xlabel('Peak Firing Rate (Hz)')
    ylabel('Cell Count')
end


%% Occlusion PD response
% Direction index of each cell's PD
dir_idx_PD = all_cells_prefdir(:,2)/45+1; % Convert directions in degrees to indices 

occ_cells_peakFR = NaN(3,Ncells);           % Initialize
for icell = occ_cells_idx
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
    
    % Load data
    fn = strcat('occlusion_cell',num2str(icell)); % Filename
    filepath = strcat('C:\Users\alber\Dropbox\Occlusion Project\sampling_distributions\',fn);
    load(filepath) % Load data
    
    data = SortedDSData.On;                   % Use 'On' tab for occlusion PD response
    ktrials = size(data.spikeCounts,2);       %repetitions
    spikeTimes = data.spikeTimes;             % Spike times for all directions
    
    heatmap = NaN(num_dir,length(edges));     % Initialize heat map
    for mdir = 1:num_dir
        % Loop through all directions
        dirSpikeTimes = spikeTimes{1,mdir};   % Get spike times for one direction
        dirSpikeCount = NaN(ktrials,length(edges)); % Initialize
        for ktrial = 1:ktrials
            % Loop through all trials
            % Bin spike times for k-th trial
            ktrialSpikeTimes = dirSpikeTimes{ktrial,1};
            dirSpikeCount(ktrial,:) = histc(ktrialSpikeTimes,edges);
            
            % Moving mean
            %dirSpikeCount(ktrial,:) = movmean(dirSpikeCount(ktrial,:),3);
            
        end
        heatmap(mdir,:) = mean(dirSpikeCount);      % Average across trials
    end
    heatmap = heatmap/binSize;                % Convert spike count to firing rate
    icell_peakFR_alldir = max(heatmap,[],2);  % Peak firing rate of i-th cell in all directions
    % Peak firing rate of i-th cell in PD and PD-adjacent directions
    icell_peakFR = icell_peakFR_alldir(PD_neighbors);
    occ_cells_peakFR(:,icell) = icell_peakFR;  % Store peak firing rate
end
occ_peakFR_PD = occ_cells_peakFR(2,:);                % Get PD peak firing rates
occ_peakFR_PD = occ_peakFR_PD(~isnan(occ_peakFR_PD)); % Remove NaN values
% Save PD peak firing rates
if save_peakFR
    save('occ_peakFR_PD','occ_peakFR_PD')
end

binranges = 0:histBin:max(occ_cells_peakFR,[],'all')+5;
occ_peakFR_distribution = NaN(length(binranges),3); % Initialize
for mdir = 1:3
    % Loop through PD and PD-adjacent directions
    occ_peakFR_distribution(:,mdir) = histc(occ_cells_peakFR(mdir,:),binranges);
end

% Fit gamma distribution to histogram
if fitDist
    occ_gammaFit_peakFRdist_PDonly = fitdist(occ_cells_peakFR(2,:)','gamma')
    figure
    histfit(occ_cells_peakFR(2,:),9,'Gamma')
    xlabel('Peak firing rate (Hz)')
    ylabel('Cell count')
    title('Occlusion Preferred direction (PD) only')
end

% Combine PD left- and right-adjacent neighbors
occ_peakFR_neighbors = [occ_cells_peakFR(1,:) occ_cells_peakFR(end,:)];
occ_peakFR_distribution_neighbors = histc(occ_peakFR_neighbors,binranges);

if plot_occ
    figure
    bar(binranges,occ_peakFR_distribution(:,1),'histc')
    title(strcat('Occlusion Preferred Direction Left Adjacent',...
        ' (bin size =',{' '},num2str(binSize*1000),' ms)'))
    xlabel('Peak Firing Rate (Hz)')
    ylabel('Cell Count')
    figure
    bar(binranges,occ_peakFR_distribution(:,2),'histc')
    title(strcat('Occlusion Preferred Direction',...
        ' (bin size =',{' '},num2str(binSize*1000),' ms)'))
    xlabel('Peak Firing Rate (Hz)')
    ylabel('Cell Count')
    figure
    bar(binranges,occ_peakFR_distribution(:,3),'histc')
    title(strcat('Occlusion Preferred Direction Right Adjacent',...
        ' (bin size =',{' '},num2str(binSize*1000),' ms)'))
    xlabel('Peak Firing Rate (Hz)')
    ylabel('Cell Count')
    figure % PD neighbors
    bar(binranges,occ_peakFR_distribution_neighbors,'histc')
    title(strcat('Occlusion Preferred Direction Adjacent Neighbors',...
        ' (bin size =',{' '},num2str(binSize*1000),' ms)'))
    xlabel('Peak Firing Rate (Hz)')
    ylabel('Cell Count')
end


%% Occlusion null response (ONR)
% Direction index of each cell's null direction
dir_idx_antiPD = mod(all_cells_prefdir(:,2)+180,360)/45+1; % Convert directions in degrees to indices

onr_cells_peakFR = NaN(3,Ncells);                   % Initialize
for icell = occ_cells_idx
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
    
    % Load data
    fn = strcat('occlusion_cell',num2str(icell)); % Filename
    filepath = strcat('C:\Users\alber\Dropbox\Occlusion Project\sampling_distributions\',fn);
    load(filepath) % Load data
    
    data = SortedDSData.Off;                  % Use 'Off' tab for null direction
    ktrials = size(data.spikeCounts,2);       %repetitions
    spikeTimes = data.spikeTimes;             % Spike times for all directions
    
    heatmap = NaN(num_dir,length(edges));     % Initialize heat map
    for mdir = 1:num_dir
        % Loop through all directions
        dirSpikeTimes = spikeTimes{1,mdir};    % Get spike times for one direction
        dirSpikeCount = NaN(ktrials,length(edges)); % Initialize
        for ktrial = 1:ktrials
            % Loop through all trials
            % Bin spike times for k-th trial
            ktrialSpikeTimes = dirSpikeTimes{ktrial,1};
            dirSpikeCount(ktrial,:) = histc(ktrialSpikeTimes,edges);
            
            % Moving mean
            %dirSpikeCount(ktrial,:) = movmean(dirSpikeCount(ktrial,:),3);
            
        end
        heatmap(mdir,:) = mean(dirSpikeCount); % Average across trials
    end
    heatmap = heatmap/binSize;                % Convert spike count to firing rate
    icell_peakFR_alldir = max(heatmap,[],2);  % Peak firing rate of i-th cell in all directions
    % Peak firing rate of i-th cell in anti-PD and anti-PD-adjacent directions
    icell_peakFR = icell_peakFR_alldir(antiPD_neighbors);
    onr_cells_peakFR(:,icell) = icell_peakFR;  % Store peak firing rate
end
onr_peakFR_ND = onr_cells_peakFR(2,:);                % Get PD peak firing rates
onr_peakFR_ND = onr_peakFR_ND(~isnan(onr_peakFR_ND)); % Remove NaN values
% Save PD peak firing rates
if save_peakFR
    save('onr_peakFR_ND','onr_peakFR_ND')
end

binranges = 0:histBin:max(onr_cells_peakFR,[],'all')+5;
onr_peakFR_distribution = NaN(length(binranges),3); % Initialize
for mdir = 1:3
    % Loop through PD and PD-adjacent directions
    onr_peakFR_distribution(:,mdir) = histc(onr_cells_peakFR(mdir,:),binranges);
end

% Fit gamma distribution to histogram
if fitDist
    onr_gammaFit_peakFRdist_PDonly = fitdist(onr_cells_peakFR(2,:)','normal')
    figure
    histfit(onr_cells_peakFR(2,:),9,'Normal')
    xlabel('Peak firing rate (Hz)')
    ylabel('Cell count')
    title('Occlusion null direction only')
end

% Combine null direction left- and right-adjacent neighbors
onr_peakFR_neighbors = [onr_cells_peakFR(1,:) onr_cells_peakFR(end,:)];
onr_peakFR_distribution_neighbors = histc(onr_peakFR_neighbors,binranges);

if plot_onr
    figure
    bar(binranges,onr_peakFR_distribution(:,1),'histc')
    title(strcat('Occlusion Null Direction Left Adjacent',...
        ' (bin size =',{' '},num2str(binSize*1000),' ms)'))
    xlabel('Peak Firing Rate (Hz)')
    ylabel('Cell Count')
    figure
    bar(binranges,onr_peakFR_distribution(:,2),'histc')
    title(strcat('Occlusion Null Direction',...
        ' (bin size =',{' '},num2str(binSize*1000),' ms)'))
    xlabel('Peak Firing Rate (Hz)')
    ylabel('Cell Count')
    figure
    bar(binranges,onr_peakFR_distribution(:,3),'histc')
    title(strcat('Occlusion Null Direction Right Adjacent',...
        ' (bin size =',{' '},num2str(binSize*1000),' ms)'))
    xlabel('Peak Firing Rate (Hz)')
    ylabel('Cell Count')
    figure % Null direction neighbors
    bar(binranges,onr_peakFR_distribution_neighbors,'histc')
    title(strcat('Occlusion Null Direction Adjacent Neighbors',...
        ' (bin size =',{' '},num2str(binSize*1000),' ms)'))
    xlabel('Peak Firing Rate (Hz)')
    ylabel('Cell Count')
end


%% Correlation between peak firing rate and onset time
if onset_peak_corr
    % 08/27/2019
    
    % Full-field
    % Load data
    load('ff_peakFR_PD')    % Load peak firing rates
    load('ff_cells_onset')  % Load onset times
    ff_cells_idx = [1:18 20:25 32:47 49:52 54:60 63:64 66:79 82:83 87:90];  
    ff_peakFR_PD(19) = [];  % Exclude Cell 19
    ff_onset_PD = ff_cells_onset(2,:);              % Get onset times for PD only
    ff_onset_PD = ff_onset_PD(~isnan(ff_onset_PD)); % Remove NaN values
    R = corrcoef(ff_onset_PD,ff_peakFR_PD); 
    rho = R(2,1);   % Pearson correlation coefficient
    ff_onset_peakFR_corr = polyfit(ff_onset_PD,ff_peakFR_PD,1); % Parameters of linear fit
    % Visualization
    figure
    hold on
    scatter(ff_onset_PD,ff_peakFR_PD,'filled')
    text(1.5,70,strcat('\rho = ',{' '},num2str(rho)))
    xlim([startTime endTime])
    %ylim([0 130])
    xlabel('Spike response onset time (sec)')
    ylabel('Peak firing rate (Hz)')
    title('Preferred direction (PD) only')

    R = corrcoef(ff_cells_onset(:,ff_cells_idx),ff_cells_peakFR(:,ff_cells_idx)); 
    rho = R(2,1);   % Pearson correlation coefficient
    figure
    hold on
    scatter(ff_cells_onset(1,ff_cells_idx),ff_cells_peakFR(1,ff_cells_idx),'filled',...
        'MarkerFaceColor','k')%[0.4660 0.6740 0.1880])
    scatter(ff_cells_onset(3,ff_cells_idx),ff_cells_peakFR(3,ff_cells_idx),'filled',...
        'MarkerFaceColor','k')%[0.4660 0.6740 0.1880])
    scatter(ff_cells_onset(2,ff_cells_idx),ff_cells_peakFR(2,ff_cells_idx),'filled','k')
    text(1,100,strcat('\rho = ',{' '},num2str(rho)))
    xlim([0 1.2])
    ylim([30 140])
    xlabel('Spike response onset time (sec)')
    ylabel('Peak firing rate (Hz)')
    title('PD and PD-adjacent directions')

    % Occlusion null response
    % Load data
    load('onr_cells_onset')  % Load onset times
    occ_cells_idx = [1:16 18 21 23 26:45 47:48 53 58:64 66 68 71:73 76 79:91];  
    onr_peakFR_ND = onr_cells_peakFR(2,:);  
    onr_peakFR_ND = onr_peakFR_ND(occ_cells_idx); % Exclude Cell 46 and 55
    onr_onset_ND = onr_cells_onset(2,:); % Get onset times for ND only
    onr_onset_ND = onr_onset_ND(~isnan(onr_onset_ND)); % Remove NaN values
    R = corrcoef(onr_onset_ND,onr_peakFR_ND); 
    rho = R(2,1);   % Pearson correlation coefficient
    onr_onset_peakFR_corr = polyfit(ff_onset_PD,ff_peakFR_PD,1); % Parameters of linear fit
    % Visualization
    figure
    hold on
    scatter(onr_onset_ND,onr_peakFR_ND,'filled')
    text(1.6,70,strcat('\rho = ',{' '},num2str(rho)))
    xlim([startTime endTime])
    %ylim([0 130])
    xlabel('Spike response onset time (sec)')
    ylabel('Peak firing rate (Hz)')
    title('Null direction only')

    tempOnset = onr_cells_onset(:,occ_cells_idx);
    tempPeak = onr_cells_peakFR(:,occ_cells_idx);
    R = corrcoef(tempOnset(~isnan(tempOnset)),tempPeak(~isnan(tempOnset))); 
    rho = R(2,1);   % Pearson correlation coefficient
    figure
    hold on
    scatter(onr_cells_onset(1,occ_cells_idx),onr_cells_peakFR(1,occ_cells_idx),'filled',...
        'MarkerFaceColor','k')%[0.4660 0.6740 0.1880])
    scatter(onr_cells_onset(3,occ_cells_idx),onr_cells_peakFR(3,occ_cells_idx),'filled',...
        'MarkerFaceColor','k')%[0.4660 0.6740 0.1880])
    scatter(onr_cells_onset(2,occ_cells_idx),onr_cells_peakFR(2,occ_cells_idx),'filled','k')
    text(1.6,70,strcat('\rho = ',{' '},num2str(rho)))
    xlim([1.2 1.8])
    ylim([0 130])
    xlabel('Spike response onset time (sec)')
    ylabel('Peak firing rate (Hz)')
    title('Null and null-adjacent directions')
end
