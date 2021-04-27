%{
Plot tuning curves for all 73 full-field cells
Rescale tuning curve height by peak firing rate
Rescale tuning curve width by circular variance

Plot peak of tuning curve vs. square-root of circular variance

Use CircStat toolbox by Philipp Berens 
P. Berens, CircStat: A Matlab Toolbox for Circular Statistics, Journal of Statistical Software, Volume 31, Issue 10, 2009 
%}

%% Parameters
% List of full-field cell indices
ff_cells_idx = [1:18 20:25 32:47 49:52 54:60 63:64 66:79 82:83 87:90];
ff_Ncells = length(ff_cells_idx);       % Number of full-field cells
Ncells = max(ff_cells_idx);             % Total number of cells

num_dir = 8;                % Number of directions

theta_axis = -180:45:135;   % Horizontal (theta) axis

% Plotting
plot_collapsedTuningCurves = 0;     % Plot collapsed tuning curves?
plot_gaussianFit = 1;               % Plot Gaussian fit to tuning curves?
plot_corrPeakWidth = 0;             % Plot peak-width correlation?

%% Full-field
% Direction index of each cell's PD
cd('C:\Users\alber\Dropbox\Occlusion Project\sampling_distributions') % Specify directory
load('all_cells_prefdir')                                             % Load file
% First column for full-field
% Second column for occlusion
dir_idx_PD = all_cells_prefdir(:,1)/45+1; % Convert directions in degrees to indices 

% Initialize
ff_cells_meanSpikeCount_aligned = NaN(num_dir+num_dir/2,Ncells); % All 8 directions aligned
                                                                 % 12 by Ncells
ff_tuningWidths = NaN(Ncells,1);          % Initialize array of circular variances
for icell = ff_cells_idx
    icell_dir_idx_PD = dir_idx_PD(icell); % Direction index of i-th cell's PD
    
    % Load data
    fn = strcat('fullfield_cell',num2str(icell));   % Filename
    filepath = strcat('C:\Users\alber\Dropbox\Occlusion Project\sampling_distributions\',fn);
    load(filepath) % Load data
    
    data = SortedDSData.total;                      % Get data structure
    directions = circ_ang2rad(data.directions);     % Directions in radians
    ktrials = size(data.spikeCounts,2);             % Number of repetitions
    spikeCounts = data.spikeCounts;                 % Spike counts for all directions
    
    meanSpikeCount = mean(spikeCounts,2);           % Trial-averaged spike count for each direction

    % Align PD and store mean spike count of i-th cell
    % 5-th entry is PD
    ff_cells_meanSpikeCount_aligned(num_dir/2+1:end,icell) = ...
        [meanSpikeCount(icell_dir_idx_PD:end); ...          % Right of PD
        meanSpikeCount(1:icell_dir_idx_PD-1)];              % Left of PD
    
    icell_distribution = [];    % Initialize data array
    % Want to build a distribution of directions where the height is given
    % by mean spike count
    for mdir = 1:num_dir
        % Weight direction by mean spike count
        % Repeat direction (in radians) by mean spike count (rounded to nearest integer)
        dirWeighted = repmat(directions(mdir),round(meanSpikeCount(mdir)),1);
        % Append to data array
        icell_distribution = [icell_distribution; dirWeighted];
    end
    ff_tuningWidths(icell) = circ_var(icell_distribution);  % Store circular variance for i-th cell
end
% Align PD
ff_cells_meanSpikeCount_aligned(1:num_dir/2,:) = ...        % Chop off anything beyond 8-th entry and
    ff_cells_meanSpikeCount_aligned(num_dir+1:end,:);       % paste into 1st-4th entries
ff_cells_meanSpikeCount_aligned = ...
    ff_cells_meanSpikeCount_aligned(1:num_dir,:);           % 8 by Ncells
maxMeanSpikeCount = max(ff_cells_meanSpikeCount_aligned);   % Max mean spike count of each cell
ff_cells_meanSpikeCount_aligned = ...                       % Rescale (divide) by PD raw spike count
    ff_cells_meanSpikeCount_aligned./repmat(maxMeanSpikeCount,num_dir,1);

% Fit Gaussian to tuning curves
ydataGaussianFit = reshape(ff_cells_meanSpikeCount_aligned,...  % Reshape spike counts for each cell
    1,[prod(size(ff_cells_meanSpikeCount_aligned))]);           % into a single column
ydataGaussianFit = ydataGaussianFit(~isnan(ydataGaussianFit));  % Remove NaN
xdataGaussianFit = [];      % Initialize x data array for Gaussian fit
for icell = ff_cells_idx
    xdataGaussianFit = [xdataGaussianFit ...        % Concatenate theta (x) axis
        theta_axis/sqrt(2*ff_tuningWidths(icell))]; % rescaled by angular deviation
end
GaussianFit = fit(xdataGaussianFit',ydataGaussianFit','gauss1') % Display fitted parameters
                                                                % in command window

%% Visualization
if plot_collapsedTuningCurves
    figure
    hold on % Plot multiple tuning curves on top of one another
    for icell = ff_cells_idx
        plot(theta_axis/sqrt(2*ff_tuningWidths(icell)),...      % Line plot
            ff_cells_meanSpikeCount_aligned(:,icell))
        %scatter(theta_axis/sqrt(2*ff_tuningWidths(icell)),...  % Scatter plot
        %    ff_cells_meanSpikeCount_aligned(:,icell))
    end
    xlabel('Deviation from PD (Angular deviation)')
    ylabel('Spike Count (Units of max spike count)')
    ylim([-0.1 1.2])
    title(strcat('Rescaled Tuning Curves of All Cells (N=',num2str(ff_Ncells),')'))
end
% Plot Gaussian fit to tuning curves
if plot_gaussianFit
    figure
    hold on % Plot individual traces first
    for icell = ff_cells_idx
        plot(theta_axis/sqrt(2*ff_tuningWidths(icell)),...      % Line plot
            ff_cells_meanSpikeCount_aligned(:,icell))
    end
    % Plot Gaussian fit
    h = plot(GaussianFit,xdataGaussianFit,ydataGaussianFit);
    set(h,'LineWidth',2.5)
    legend(gca,'off')   % Turn off legend
    xlabel('Deviation from PD (Angular deviation)')
    ylabel('Spike Count (Units of max spike count)')
    ylim([-0.1 1.2])
    title(strcat('Gaussian Fit to Tuning Curves of All Cells (N=',num2str(ff_Ncells),')'))
    
end

% Find parameters of linear fit
p = polyfit(maxMeanSpikeCount(ff_cells_idx),ff_tuningWidths(ff_cells_idx)',1);
% Find correlation coefficient
R = corrcoef(maxMeanSpikeCount(ff_cells_idx),ff_tuningWidths(ff_cells_idx)');
rho = R(2,1);   % Pearson correlation coefficient
if plot_corrPeakWidth
    % Plot correlation between tuning curve peak and width
    xaxis = linspace(0,max(maxMeanSpikeCount)+10);  % Create x-values of linear fit
    linFit = polyval(p,xaxis);                      % Create y-values of linear fit
    
    figure
    hold on
    % Plot scatter plot
    scatter(maxMeanSpikeCount(ff_cells_idx),ff_tuningWidths(ff_cells_idx),'filled')
    % Plot linear fit
    plot(xaxis,linFit,'k')
    % Display Pearson correlation coefficient
    text(70,0.35,strcat('\rho =',{' '},num2str(rho)),'FontSize',11)
    xlabel('Max Spike Count')
    ylabel('Tuning Width (Square-root of Circular Variance)')
    title('Tuning Curve Peak-Width Correlation')
end


