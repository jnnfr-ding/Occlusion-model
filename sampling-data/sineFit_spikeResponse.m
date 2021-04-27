%% Parameters
% List of full-field cell indices
%ff_cells_idx = [1:18 20:25 32:47 49:52 54:60 63:64 66:79 82:83 87:90];
ff_cells_idx = [25 63 89];              % Indices of manually adjusted cells
                                        % Exclude Cell 19!!!
ff_Ncells = length(ff_cells_idx);       % Number of full-field cells
% List of occlusion cell indices
occ_cells_idx = 1; %[1:16 18 21 23 26:45 47:48 53 55 58:64 66 68 71:73 76 79:91];
%occ_cells_idx = [12 28 30 39 65];       % Indices of manually adjusted cells
                                        % Exclude Cell 46!!!
occ_Ncells = length(occ_cells_idx);     % Number of occlusion cells
Ncells = max(max(ff_cells_idx),max(occ_cells_idx));

num_dir = 8;           % Number of directions

binSize = 0.025;       %sec
startTime = 0;         %sec
endTime = 2;           %sec
edges = startTime:binSize:endTime; % Create time bins
Nbins = length(edges); % Number of bins

% Figures
% Full-field
plot_ff_psth = 1;      % Plot PSTH?
save_ff_figures = 0;   % Save figures?
% Occlusion
plot_occ_psth = 1;
save_occ_figures = 0;

%% Load preferred directions
cd('C:\Users\alber\Dropbox\Occlusion Project\sampling_distributions') % Specify directory
% Direction index of each cell's PD
load('all_cells_prefdir')   % Load file containing preferred directions
                            % First column for full-field
                            % Second column for occlusion
                            
%% Full-field
dir_idx_PD = all_cells_prefdir(:,1)/45+1; % Convert directions in degrees to indices 
% Spike response onset time
load('ff_cells_onset_aligned')  % Row:      8 directions
                                % Column:   Ncells
                                % PD is 5-th row
% Load spike response onset-duration correlation
% PD only
%load('onset_duration_corr')     % 1st variable: slope
                                % 2nd variable: vertical offset
% PD and PD-adjacent
load('ff_onset_duration_corr')  % 1st variable: slope
                                % 2nd variable: vertical offset
% Load offset times
load('ff_cells_offset')         % PD and PD-adjacent
                                
for icell = ff_cells_idx
    icell_dir_idx_PD = dir_idx_PD(icell);     % PD of i-th cell
    icell_onset_time = ff_cells_onset_aligned(:,icell); % Onset time for i-th cell
    
    % Load data
    fn = strcat('fullfield_cell',num2str(icell)); % Filename
    filepath = strcat('C:\Users\alber\Dropbox\Occlusion Project\sampling_distributions\',fn);
    load(filepath) % Load data
    
    data = SortedDSData.total;                % Get data structure
    ktrials = size(data.spikeCounts,2);       %repetitions
    spikeTimes = data.spikeTimes;             % Spike times for all directions

    heatmap_aligned = NaN(1.5*num_dir,Nbins); % Initialize heat map (12 x Nbins)
    heatmap = NaN(num_dir,Nbins);             % Initialize heat map (8 x Nbins)
    for mdir = 1:num_dir
        % Loop through all directions
        dirSpikeTimes = spikeTimes{1,mdir};   % Get spike times for one direction
        dirSpikeCount = NaN(ktrials,Nbins);   % Initialize PSTH for one direction
        for ktrial = 1:ktrials
            % Loop through all trials
            % Bin spike times for k-th trial
            ktrialSpikeTimes = dirSpikeTimes{ktrial,1};
            dirSpikeCount(ktrial,:) = histc(ktrialSpikeTimes,edges);
        end
        heatmap(mdir,:) = mean(dirSpikeCount); % Average across trials
    end
    heatmap = heatmap/binSize;                % Convert spike count to firing rate
    
    % Align PD
    heatmap_aligned(num_dir/2+1:end,:) = ...
        [heatmap(icell_dir_idx_PD:end,:); ... % Right of PD
        heatmap(1:icell_dir_idx_PD-1,:)];     % Left of PD
    heatmap_aligned(1:num_dir/2,:) = heatmap_aligned(num_dir+1:end,:);
    heatmap_aligned = heatmap_aligned(1:num_dir,:); % Heat map aligned by PD (5-th row)
    
    % Sine fit
    sinFit_param = NaN(num_dir,3);  % Initialize array sine fit parameters
                                    % 1st column: Amplitude (peak firing)
                                    % 2nd column: Period (two times spike response duration)
                                    % 3rd column: Phase (spike response onset time)
    sinFit_param(:,1) = max(heatmap_aligned,[],2);  % Get and store peak firing
    sinFit_param(:,3) = icell_onset_time;           % Get and store onset times
    sinFit_param(:,2) = 2.*(ff_onset_duration_corr(2)+...  % Convert onset time to spike response duration
        icell_onset_time.*ff_onset_duration_corr(1));      % Convert spike response duration to period
        
    % Plotting
    if plot_ff_psth
        
        figure
        plot(edges,heatmap_aligned(5,:))          % PD
        hold on
        s5 = sineFit(sinFit_param(5,:));          % 5th row
        s5_edges = sinFit_param(5,3):binSize:...  % Time bins for sine fit
            sinFit_param(5,3)+sinFit_param(5,2)/2+binSize;
        plot(s5_edges,s5(s5_edges))               % Plot sine fit
        text(0.2,125,num2str(ff_cells_onset_aligned(5,icell))) % Display onset time
        text(1.2,125,num2str(ff_cells_offset(2,icell)))        % Display offset time
        xlabel('Time (sec)')
        ylabel('Firing Rate (Hz)')
        ylim([0 150])
        
        figure
        plot(edges,heatmap_aligned(4,:))          % PD-1
        hold on     % Add sine fit
        s4 = sineFit(sinFit_param(4,:));          % 4th row
        s4_occ_edges = sinFit_param(4,3):binSize:...  % Time bins for sine fit
            sinFit_param(4,3)+sinFit_param(4,2)/2+binSize;
        plot(s4_occ_edges,s4(s4_occ_edges))       % Plot sine fit
        text(0.2,125,num2str(ff_cells_onset_aligned(4,icell))) % Display onset time
        text(1.2,125,num2str(ff_cells_offset(1,icell)))        % Display offset time
        xlabel('Time (sec)')
        ylabel('Firing Rate (Hz)')
        ylim([0 150])
        
        figure
        plot(edges,heatmap_aligned(6,:))          % PD+1
        hold on
        s6 = sineFit(sinFit_param(6,:));          % 6th row
        s6_edges = sinFit_param(6,3):binSize:...  % Time bins for sine fit
            sinFit_param(6,3)+sinFit_param(6,2)/2+binSize;
        plot(s6_edges,s6(s6_edges))               % Plot sine fit
        text(0.2,125,num2str(ff_cells_onset_aligned(6,icell))) % Display onset time
        text(1.2,125,num2str(ff_cells_offset(3,icell)))        % Display offset time
        xlabel('Time (sec)')
        ylabel('Firing Rate (Hz)')
        ylim([0 150])
        
        
        % PD is always pointed along 90 degrees
        figure                                    % PSTH and sine fit
        hold on
        
        subplot(3,3,1)                            % 135 degrees
        plot(edges,heatmap_aligned(4,:))          % PD-1
        hold on     % Add sine fit
        s4 = sineFit(sinFit_param(4,:));          % 4th row
        s4_occ_edges = sinFit_param(4,3):binSize:...  % Time bins for sine fit
            sinFit_param(4,3)+sinFit_param(4,2)/2+binSize;
        plot(s4_occ_edges,s4(s4_occ_edges))       % Plot sine fit
        text(0.2,125,num2str(ff_cells_onset_aligned(4,icell))) % Display onset time
        text(1.2,125,num2str(ff_cells_offset(1,icell)))        % Display offset time
        xlabel('Time (sec)')
        ylabel('Firing Rate (Hz)')
        ylim([0 150])
                     
        subplot(3,3,2)                            % 90 degrees
        plot(edges,heatmap_aligned(5,:))          % PD
        hold on
        s5 = sineFit(sinFit_param(5,:));          % 5th row
        s5_edges = sinFit_param(5,3):binSize:...  % Time bins for sine fit
            sinFit_param(5,3)+sinFit_param(5,2)/2+binSize;
        plot(s5_edges,s5(s5_edges))               % Plot sine fit
        text(0.2,125,num2str(ff_cells_onset_aligned(5,icell))) % Display onset time
        text(1.2,125,num2str(ff_cells_offset(2,icell)))        % Display offset time
        xlabel('Time (sec)')
        ylabel('Firing Rate (Hz)')
        ylim([0 150])
        
        subplot(3,3,3)                            % 45 degrees
        plot(edges,heatmap_aligned(6,:))          % PD+1
        hold on
        s6 = sineFit(sinFit_param(6,:));          % 6th row
        s6_edges = sinFit_param(6,3):binSize:...  % Time bins for sine fit
            sinFit_param(6,3)+sinFit_param(6,2)/2+binSize;
        plot(s6_edges,s6(s6_edges))               % Plot sine fit
        text(0.2,125,num2str(ff_cells_onset_aligned(6,icell))) % Display onset time
        text(1.2,125,num2str(ff_cells_offset(3,icell)))        % Display offset time
        xlabel('Time (sec)')
        ylabel('Firing Rate (Hz)')
        ylim([0 150])
        
        subplot(3,3,4)                            % 180 degrees
        plot(edges,heatmap_aligned(3,:))          % PD-2
        hold on
        s3 = sineFit(sinFit_param(3,:));          % 3rd row
        s3_edges = sinFit_param(3,3):binSize:...  % Time bins for sine fit
            sinFit_param(3,3)+sinFit_param(3,2)/2+binSize;
        plot(s3_edges,s3(s3_edges))               % Plot sine fit
        xlabel('Time (sec)')
        ylabel('Firing Rate (Hz)')
        ylim([0 150])

        subplot(3,3,6)                            % 0 degrees
        plot(edges,heatmap_aligned(7,:))          % PD+2
        hold on
        s7 = sineFit(sinFit_param(7,:));          % 7th row
        s7_edges = sinFit_param(7,3):binSize:...  % Time bins for sine fit
            sinFit_param(7,3)+sinFit_param(7,2)/2+binSize;
        plot(s7_edges,s7(s7_edges))               % Plot sine fit
        xlabel('Time (sec)')
        ylabel('Firing Rate (Hz)')
        ylim([0 150])
        
        subplot(3,3,7)                            % 225 degrees
        plot(edges,heatmap_aligned(2,:))          % PD-3
        hold on
        s2 = sineFit(sinFit_param(2,:));          % 2nd row
        s2_edges = sinFit_param(2,3):binSize:...  % Time bins for sine fit
            sinFit_param(2,3)+sinFit_param(2,2)/2+binSize;
        plot(s2_edges,s2(s2_edges))               % Plot sine fit
        xlabel('Time (sec)')
        ylabel('Firing Rate (Hz)')
        ylim([0 150])
        
        subplot(3,3,8)                            % 270 degrees
        plot(edges,heatmap_aligned(1,:))          % PD-4
        hold on
        s1 = sineFit(sinFit_param(1,:));          % 1st row
        s1_edges = sinFit_param(1,3):binSize:...  % Time bins for sine fit
            sinFit_param(1,3)+sinFit_param(1,2)/2+binSize;
        plot(s1_edges,s1(s1_edges))               % Plot sine fit
        xlabel('Time (sec)')
        ylabel('Firing Rate (Hz)')
        ylim([0 150])
        
        subplot(3,3,9)                            % 315 degrees
        plot(edges,heatmap_aligned(8,:))          % PD+3
        hold on
        s8 = sineFit(sinFit_param(8,:));          % 8th row
        s8_edges = sinFit_param(8,3):binSize:...  % Time bins for sine fit
            sinFit_param(8,3)+sinFit_param(8,2)/2+binSize;
        plot(s8_edges,s8(s8_edges))               % Plot sine fit
        xlabel('Time (sec)')
        ylabel('Firing Rate (Hz)')
        ylim([0 150])
        
        sgtitle(strcat('Full-field',' Cell',num2str(icell),...
            ' (bin size =',{' '},num2str(binSize*1000),' ms)'))
        if save_ff_figures
            saveas(gcf,strcat('ff_PSTH_alldir_cell',num2str(icell)),'jpg')
        end
    end
end

%% Occlusion
% Convert directions in degrees to indices
dir_idx_PD = all_cells_prefdir(:,2)/45+1;                  % PD
dir_idx_antiPD = mod(all_cells_prefdir(:,2)+180,360)/45+1; % Null direction

% Spike response onset time
load('occ_cells_onset_aligned') % Row:      8 directions
                                % Column:   Ncells
                                % PD is 5th row
load('onr_cells_onset_aligned') % NULL DIRECTION IS 5TH ROW!!!!!

% Load spike response onset-duration correlation
load('occ_onset_duration_corr') % 1st variable: slope
                                % 2nd variable: vertical offset
load('onr_onset_duration_corr')

% Spike response offset time
load('occ_cells_offset','occ_cells_offset') % Occlusion PD
                                            % PD is 5th row
load('onr_cells_offset','onr_cells_offset') % Occlusion null response (ONR)
                                            % NULL DIRECTION IS 5TH ROW!!!
for icell = occ_cells_idx    
    icell_dir_idx_PD = dir_idx_PD(icell);           % PD of i-th cell
    icell_dir_idx_antiPD = dir_idx_antiPD(icell);   % Null direction of i-th cell
    icell_occ_onset_time = occ_cells_onset_aligned(:,icell); % Onset time for occlusion PD response
    icell_onr_onset_time = onr_cells_onset_aligned(:,icell); % Onset time for occlusion null reponse
    
    % Load data
    fn = strcat('occlusion_cell',num2str(icell)); % Filename
    filepath = strcat('C:\Users\alber\Dropbox\Occlusion Project\sampling_distributions\',fn);
    load(filepath) % Load data
    
    %---------------------------------------------------------------------%
    % Occlusion PD
    data = SortedDSData.On;                   % Get data structure
                                              % Use 'On' tab!!!
    ktrials = size(data.spikeCounts,2);       %repetitions
    spikeTimes = data.spikeTimes;             % Spike times for all directions

    heatmap_aligned = NaN(1.5*num_dir,Nbins); % Initialize heat map (12 x Nbins)
    heatmap = NaN(num_dir,Nbins);             % Initialize heat map (8 x Nbins)
    for mdir = 1:num_dir
        % Loop through all directions
        dirSpikeTimes = spikeTimes{1,mdir};   % Get spike times for one direction
        dirSpikeCount = NaN(ktrials,Nbins);   % Initialize PSTH for one direction
        for ktrial = 1:ktrials
            % Loop through all trials
            % Bin spike times for k-th trial
            ktrialSpikeTimes = dirSpikeTimes{ktrial,1};
            dirSpikeCount(ktrial,:) = histc(ktrialSpikeTimes,edges);
        end
        heatmap(mdir,:) = mean(dirSpikeCount); % Average across trials
    end
    heatmap = heatmap/binSize;                % Convert spike count to firing rate
    
    % Align PD
    heatmap_aligned(num_dir/2+1:end,:) = ...
        [heatmap(icell_dir_idx_PD:end,:); ... % Right of PD
        heatmap(1:icell_dir_idx_PD-1,:)];     % Left of PD
    heatmap_aligned(1:num_dir/2,:) = heatmap_aligned(num_dir+1:end,:);
    heatmap_aligned = heatmap_aligned(1:num_dir,:); % Heat map aligned by PD (5-th row)
    
    % Sine fit
    occ_sinFit_param = NaN(num_dir,3);  % Initialize array sine fit parameters
                                        % 1st column: Amplitude (peak firing)
                                        % 2nd column: Period (two times spike response duration)
                                        % 3rd column: Phase (spike response onset time)
    occ_sinFit_param(:,1) = max(heatmap_aligned,[],2);  % Get and store peak firing
    occ_sinFit_param(:,3) = icell_occ_onset_time;       % Get and store onset times
    occ_sinFit_param(:,2) = 2.*(occ_onset_duration_corr(2)+...  % Convert onset time to spike response duration
        icell_occ_onset_time.*occ_onset_duration_corr(1));      % Convert spike response duration to period
        
    %---------------------------------------------------------------------%
    % Occlusion null response (ONR)
    % Null direction is 5th row!!!
    data = SortedDSData.Off;                  % Get data structure
                                              % Use 'Off' tab!!!
    ktrials = size(data.spikeCounts,2);       %repetitions
    spikeTimes = data.spikeTimes;             % Spike times for all directions

    heatmap_aligned = NaN(1.5*num_dir,Nbins); % Initialize heat map (12 x Nbins)
    heatmap = NaN(num_dir,Nbins);             % Initialize heat map (8 x Nbins)
    for mdir = 1:num_dir
        % Loop through all directions
        dirSpikeTimes = spikeTimes{1,mdir};   % Get spike times for one direction
        dirSpikeCount = NaN(ktrials,Nbins);   % Initialize PSTH for one direction
        for ktrial = 1:ktrials
            % Loop through all trials
            % Bin spike times for k-th trial
            ktrialSpikeTimes = dirSpikeTimes{ktrial,1};
            dirSpikeCount(ktrial,:) = histc(ktrialSpikeTimes,edges);
        end
        heatmap(mdir,:) = mean(dirSpikeCount); % Average across trials
    end
    heatmap = heatmap/binSize;                % Convert spike count to firing rate
    
    % Align PD
    heatmap_aligned(num_dir/2+1:end,:) = ...
        [heatmap(icell_dir_idx_antiPD:end,:); ... % Right of anti-PD
        heatmap(1:icell_dir_idx_antiPD-1,:)];     % Left of anti-PD
    heatmap_aligned(1:num_dir/2,:) = heatmap_aligned(num_dir+1:end,:);
    heatmap_aligned = heatmap_aligned(1:num_dir,:); % Heat map aligned by anti-PD 
                                                    % Anti-PD is 5th row!!!
    
    % Sine fit
    onr_sinFit_param = NaN(num_dir,3);  % Initialize array sine fit parameters
                                        % 1st column: Amplitude (peak firing)
                                        % 2nd column: Period (two times spike response duration)
                                        % 3rd column: Phase (spike response onset time)
    onr_sinFit_param(:,1) = max(heatmap_aligned,[],2);  % Get and store peak firing
    onr_sinFit_param(:,3) = icell_onr_onset_time;       % Get and store onset times
    onr_sinFit_param(:,2) = 2.*(onr_onset_duration_corr(2)+...  % Convert onset time to spike response duration
        icell_onr_onset_time.*onr_onset_duration_corr(1));      % Convert spike response duration to period
        
    %---------------------------------------------------------------------%
    % PSTH
    data = SortedDSData.total;                % Get data structure
                                              % Use 'total' tab!!!
    ktrials = size(data.spikeCounts,2);       %repetitions
    spikeTimes = data.spikeTimes;             % Spike times for all directions

    heatmap_aligned = NaN(1.5*num_dir,Nbins); % Initialize heat map (12 x Nbins)
    heatmap = NaN(num_dir,Nbins);             % Initialize heat map (8 x Nbins)
    for mdir = 1:num_dir
        % Loop through all directions
        dirSpikeTimes = spikeTimes{1,mdir};   % Get spike times for one direction
        dirSpikeCount = NaN(ktrials,Nbins);   % Initialize PSTH for one direction
        for ktrial = 1:ktrials
            % Loop through all trials
            % Bin spike times for k-th trial
            ktrialSpikeTimes = dirSpikeTimes{ktrial,1};
            dirSpikeCount(ktrial,:) = histc(ktrialSpikeTimes,edges);
        end
        heatmap(mdir,:) = mean(dirSpikeCount); % Average across trials
    end
    heatmap = heatmap/binSize;                % Convert spike count to firing rate
    
    % Align PD
    heatmap_aligned(num_dir/2+1:end,:) = ...
        [heatmap(icell_dir_idx_PD:end,:); ... % Right of PD
        heatmap(1:icell_dir_idx_PD-1,:)];     % Left of PD
    heatmap_aligned(1:num_dir/2,:) = heatmap_aligned(num_dir+1:end,:);
    heatmap_aligned = heatmap_aligned(1:num_dir,:); % Heat map aligned by PD (5-th row)
    
    %---------------------------------------------------------------------%
    
    
    % Plotting
    if plot_occ_psth
        
        figure
        hold on                                   % PD-4
        plot(edges,heatmap_aligned(1,:))          % (1)
        s1_occ = sineFit(occ_sinFit_param(1,:));          % 1st row
        s1_occ_edges = occ_sinFit_param(1,3):binSize:...  % Time bins for sine fit
            occ_sinFit_param(1,3)+occ_sinFit_param(1,2)/2+binSize;
        plot(s1_occ_edges,s1_occ(s1_occ_edges))   % (2)
        s5_onr = sineFit(onr_sinFit_param(5,:));          % 5th row
        s5_onr_edges = onr_sinFit_param(5,3):binSize:...
            onr_sinFit_param(5,3)+onr_sinFit_param(5,2)/2+binSize;
        plot(s5_onr_edges,s5_onr(s5_onr_edges),...
            'Color',[0.4660, 0.6740, 0.1880])     % (3)
        text(0.2,125,num2str(onr_cells_onset_aligned(5,icell))) % Display onset time for ONR antiPD
        text(1.2,125,num2str(onr_cells_offset(2,icell)))        % Display offset time for ONR antiPD
        xlabel('Time (sec)')
        ylabel('Firing Rate (Hz)')
        ylim([0 150])
        
        
        % PD is always pointed along 90 degrees
        figure                                    
        hold on
        
        % Each subplot has 
        % (1) PSTH 
        % (2) sine fit to occlusion PD response
        % (3) sine fit to occlusion null response
        
        subplot(3,3,1)                            % 135 degrees
        hold on                                   % PD-1
        plot(edges,heatmap_aligned(4,:))          % (1)
        s4_occ = sineFit(occ_sinFit_param(4,:));          % 4th row
        s4_occ_edges = occ_sinFit_param(4,3):binSize:...  % Time bins for sine fit
            occ_sinFit_param(4,3)+occ_sinFit_param(4,2)/2+binSize;
        plot(s4_occ_edges,s4_occ(s4_occ_edges))   % (2)
        s8_onr = sineFit(onr_sinFit_param(8,:));          % 8th row
        s8_onr_edges = onr_sinFit_param(8,3):binSize:...
            onr_sinFit_param(8,3)+onr_sinFit_param(8,2)/2+binSize;
        plot(s8_onr_edges,s8_onr(s8_onr_edges),...
            'Color',[0.4660, 0.6740, 0.1880])     % (3)
        text(0.2,125,num2str(occ_cells_onset_aligned(4,icell))) % Display onset time for occlusion PD-1
        text(1.2,125,num2str(occ_cells_offset(1,icell))) % Display offset time for occlusion PD-1
        xlabel('Time (sec)')
        ylabel('Firing Rate (Hz)')
        ylim([0 150])
                     
        subplot(3,3,2)                            % 90 degrees
        hold on                                   % PD
        plot(edges,heatmap_aligned(5,:))          % (1)      
        s5_occ = sineFit(occ_sinFit_param(5,:));          % 5th row
        s5_occ_edges = occ_sinFit_param(5,3):binSize:...  % Time bins for sine fit
            occ_sinFit_param(5,3)+occ_sinFit_param(5,2)/2+binSize;
        plot(s5_occ_edges,s5_occ(s5_occ_edges))   % (2)
        s1_onr = sineFit(onr_sinFit_param(1,:));          % 1st row
        s1_onr_edges = onr_sinFit_param(1,3):binSize:...
            onr_sinFit_param(1,3)+onr_sinFit_param(1,2)/2+binSize;
        plot(s1_onr_edges,s1_onr(s1_onr_edges),...
            'Color',[0.4660, 0.6740, 0.1880])    % (3)
        text(0.2,125,num2str(occ_cells_onset_aligned(5,icell))) % Display onset time for occlusion PD
        text(1.2,125,num2str(occ_cells_offset(2,icell)))        % Display offset time for occlusion PD
        xlabel('Time (sec)')
        ylabel('Firing Rate (Hz)')
        ylim([0 150])
        
        subplot(3,3,3)                            % 45 degrees
        hold on                                   % PD+1
        plot(edges,heatmap_aligned(6,:))          % (1)
        s6_occ = sineFit(occ_sinFit_param(6,:));          % 6th row
        s6_occ_edges = occ_sinFit_param(6,3):binSize:...  % Time bins for sine fit
            occ_sinFit_param(6,3)+occ_sinFit_param(6,2)/2+binSize;
        plot(s6_occ_edges,s6_occ(s6_occ_edges))   % (2)
        s2_onr = sineFit(onr_sinFit_param(2,:));          % 2nd row
        s2_onr_edges = onr_sinFit_param(2,3):binSize:...
            onr_sinFit_param(2,3)+onr_sinFit_param(2,2)/2+binSize;
        plot(s2_onr_edges,s2_onr(s2_onr_edges),...
            'Color',[0.4660, 0.6740, 0.1880])     % (3)
        text(0.2,125,num2str(occ_cells_onset_aligned(6,icell))) % Display onset time for occlusion PD+1
        text(1.2,125,num2str(occ_cells_offset(3,icell)))        % Display offset time for occlusion PD+1
        xlabel('Time (sec)')
        ylabel('Firing Rate (Hz)')
        ylim([0 150])
        
        subplot(3,3,4)                            % 180 degrees
        hold on                                   % PD-2
        plot(edges,heatmap_aligned(3,:))          % (1)
        s3_occ = sineFit(occ_sinFit_param(3,:));          % 3rd row
        s3_occ_edges = occ_sinFit_param(3,3):binSize:...  % Time bins for sine fit
            occ_sinFit_param(3,3)+occ_sinFit_param(3,2)/2+binSize;
        plot(s3_occ_edges,s3_occ(s3_occ_edges))   % (2)
        s7_onr = sineFit(onr_sinFit_param(7,:));          % 7th row
        s7_onr_edges = onr_sinFit_param(7,3):binSize:...
            onr_sinFit_param(7,3)+onr_sinFit_param(7,2)/2+binSize;
        plot(s7_onr_edges,s7_onr(s7_onr_edges),...
            'Color',[0.4660, 0.6740, 0.1880])     % (3)
        xlabel('Time (sec)')
        ylabel('Firing Rate (Hz)')
        ylim([0 150])

        subplot(3,3,6)                            % 0 degrees
        hold on                                   % PD+2
        plot(edges,heatmap_aligned(7,:))          % (1)
        s7_occ = sineFit(occ_sinFit_param(7,:));          % 7th row
        s7_occ_edges = occ_sinFit_param(7,3):binSize:...  % Time bins for sine fit
            occ_sinFit_param(7,3)+occ_sinFit_param(7,2)/2+binSize;
        plot(s7_occ_edges,s7_occ(s7_occ_edges))   % (2)
        s3_onr = sineFit(onr_sinFit_param(3,:));          % 3rd row
        s3_onr_edges = onr_sinFit_param(3,3):binSize:...
            onr_sinFit_param(3,3)+onr_sinFit_param(3,2)/2+binSize;
        plot(s3_onr_edges,s3_onr(s3_onr_edges),...
            'Color',[0.4660, 0.6740, 0.1880])     % (3)
        xlabel('Time (sec)')
        ylabel('Firing Rate (Hz)')
        ylim([0 150])
        
        subplot(3,3,7)                            % 225 degrees
        hold on                                   % PD-3
        plot(edges,heatmap_aligned(2,:))          % (1)
        s2_occ = sineFit(occ_sinFit_param(2,:));          % 2nd row
        s2_occ_edges = occ_sinFit_param(2,3):binSize:...  % Time bins for sine fit
            occ_sinFit_param(2,3)+occ_sinFit_param(2,2)/2+binSize;
        plot(s2_occ_edges,s2_occ(s2_occ_edges))   % (2)
        s6_onr = sineFit(onr_sinFit_param(6,:));          % 6th row
        s6_onr_edges = onr_sinFit_param(6,3):binSize:...
            onr_sinFit_param(6,3)+onr_sinFit_param(6,2)/2+binSize;
        plot(s6_onr_edges,s6_onr(s6_onr_edges),...
            'Color',[0.4660, 0.6740, 0.1880])     % (3)
        text(0.2,125,num2str(onr_cells_onset_aligned(6,icell))) % Display onset time for ONR antiPD-1
        text(1.2,125,num2str(onr_cells_offset(3,icell)))        % Display offset time for ONR antiPD-1
        xlabel('Time (sec)')
        ylabel('Firing Rate (Hz)')
        ylim([0 150])
        
        subplot(3,3,8)                            % 270 degrees
        hold on                                   % PD-4
        plot(edges,heatmap_aligned(1,:))          % (1)
        s1_occ = sineFit(occ_sinFit_param(1,:));          % 1st row
        s1_occ_edges = occ_sinFit_param(1,3):binSize:...  % Time bins for sine fit
            occ_sinFit_param(1,3)+occ_sinFit_param(1,2)/2+binSize;
        plot(s1_occ_edges,s1_occ(s1_occ_edges))   % (2)
        s5_onr = sineFit(onr_sinFit_param(5,:));          % 5th row
        s5_onr_edges = onr_sinFit_param(5,3):binSize:...
            onr_sinFit_param(5,3)+onr_sinFit_param(5,2)/2+binSize;
        plot(s5_onr_edges,s5_onr(s5_onr_edges),...
            'Color',[0.4660, 0.6740, 0.1880])     % (3)
        text(0.2,125,num2str(onr_cells_onset_aligned(5,icell))) % Display onset time for ONR antiPD
        text(1.2,125,num2str(onr_cells_offset(2,icell)))        % Display offset time for ONR antiPD
        xlabel('Time (sec)')
        ylabel('Firing Rate (Hz)')
        ylim([0 150])
        
        subplot(3,3,9)                            % 315 degrees
        hold on                                   % PD+3
        plot(edges,heatmap_aligned(8,:))          % (1)
        s8_occ = sineFit(occ_sinFit_param(8,:));          % 8th row
        s8_occ_edges = occ_sinFit_param(8,3):binSize:...  % Time bins for sine fit
            occ_sinFit_param(8,3)+occ_sinFit_param(8,2)/2+binSize;
        plot(s8_occ_edges,s8_occ(s8_occ_edges))   % (2)
        s4_onr = sineFit(onr_sinFit_param(4,:));          % 4th row
        s4_onr_edges = onr_sinFit_param(4,3):binSize:...
            onr_sinFit_param(4,3)+onr_sinFit_param(4,2)/2+binSize;
        plot(s4_onr_edges,s4_onr(s4_onr_edges),...
            'Color',[0.4660, 0.6740, 0.1880])     % (3)
        text(0.2,125,num2str(onr_cells_onset_aligned(4,icell))) % Display onset time for ONR antiPD+1
        text(1.2,125,num2str(onr_cells_offset(1,icell)))        % Display offset time for ONR antiPD+1
        xlabel('Time (sec)')
        ylabel('Firing Rate (Hz)')
        ylim([0 150])
        
        sgtitle(strcat('Occlusion',' Cell',num2str(icell),...
            ' (bin size =',{' '},num2str(binSize*1000),' ms)'))
        if save_occ_figures
            saveas(gcf,strcat('occ_PSTH_alldir_cell',num2str(icell)),'jpg')
        end
    end
end



