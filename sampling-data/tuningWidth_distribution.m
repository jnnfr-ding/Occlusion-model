%{
Use CircStat toolbox by Philipp Berens 
P. Berens, CircStat: A Matlab Toolbox for Circular Statistics, Journal of Statistical Software, Volume 31, Issue 10, 2009 
%}

%% Parameters
% List of full-field cell indices
ff_cells_idx = [1:18 20:25 32:47 49:52 54:60 63:64 66:79 82:83 87:90];
ff_Ncells = length(ff_cells_idx);   % Number of full-field cells
% List of occlusion cell indices
occ_cells_idx = [1:16 18 21 23 26:45 47:48 53 55 58:64 66 68 71:73 76 79:91];
occ_Ncells = length(occ_cells_idx); % Number of occlusion cells
Ncells = max(max(ff_cells_idx),max(occ_cells_idx));

num_dir = 8;            % Number of directions

binranges = 0:0.05:1;    % Edges for binning tuning widths                

%% Full-field
% Direction index of each cell's PD
cd('C:\Users\alber\Dropbox\Occlusion Project\sampling_distributions') % Specify directory
load('all_cells_prefdir')                                             % Load appropriate file
% First column for full-field
% Second column for occlusion
dir_idx_PD = all_cells_prefdir(:,1)/45+1; % Convert directions in degrees to indices 

% Initialize
ff_tuningWidths = NaN(Ncells,1);   % Initialize array of circular variances
ff_vonMises_params = NaN(Ncells,2);    % Initialize array of von Mises parameters
                                          % First column:   mu
                                          % Second column:  kappa
for icell = ff_cells_idx
    % Load data
    fn = strcat('fullfield_cell',num2str(icell));   % Filename
    filepath = strcat('C:\Users\alber\Dropbox\Occlusion Project\sampling_distributions\',fn);
    load(filepath) % Load data
    
    data = SortedDSData.total;                      % Get data structure
    directions = circ_ang2rad(data.directions);     % Directions in radians
    ktrials = size(data.spikeCounts,2);             % Number of repetitions
    spikeCounts = data.spikeCounts;                 % Spike counts for all directions
    % Trial-averaged spike count for each direction
    meanSpikeCount = mean(spikeCounts,2);
    
    icell_distribution = [];    % Initialize data array
    for mdir = 1:num_dir
        % Weight direction by mean spike count
        % Repeat direction (in radians) by mean spike count (rounded to nearest integer)
        dirWeighted = repmat(directions(mdir),round(meanSpikeCount(mdir)),1);
        % Append to data array
        icell_distribution = [icell_distribution; dirWeighted];
    end
    ff_tuningWidths(icell) = circ_var(icell_distribution);
    ff_vonMises_params(icell,:) = circ_vmpar(icell_distribution);
end

% Build distribution
ff_tuningWidth_distribution = histc(ff_tuningWidths,binranges); % Histogram
% Visualize
figure
bar(binranges,ff_tuningWidth_distribution,'histc')
xlabel('Circular Variance')
ylabel('Cell Count')
title('Distribution of Tuning Curve Widths')

% Gaussian fit to distribution
ff_gaussianFit_tuningWidth_dist = fitdist(ff_tuningWidths,'normal')
tuningWidth_gaussFit_param = ff_gaussianFit_tuningWidth_dist.ParameterValues;
save('tuningWidth_gaussFit_param','tuningWidth_gaussFit_param')
figure
histfit(ff_tuningWidths,14,'normal')
xlabel('Circular variance')
ylabel('Cell count')
title('Gaussian fit to tuning width distribution')

