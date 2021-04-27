%{
PD = preferred direction
ND = null direction
ONR = occlusion null response
%}

function runModel_coincidenceDetection(varargin)

%% Check inputs and set up default parameters.
p = inputParser;
v = @validateattributes;

addParameter(p, 'sampling_dir', 'sampling_distributions/', @(x) v(x,{'char'},{'nonempty'},mfilename,'sampling_dir')); % name of directory containing sample data .mat files
addParameter(p, 'occ_xpos',           1800, @(x) v(x,{'numeric'},{'scalar','positive'},mfilename,'occ_xpos')); % Fix this magic number (microns)!!! OR specify after setting up soma positions
addParameter(p, 'baselineFR',            0, @(x) v(x,{'numeric'},{'scalar','>=',0},mfilename,'baselineFR')); % spikes per second
addParameter(p, 'barSpeed',            330, @(x) v(x,{'numeric'},{'scalar','positive'},mfilename,'barSpeed')); % Moving bar stimulus speed (microns/sec)
addParameter(p, 'plot_ROCcurve',         1, @(x) v(x,{'numeric','logical'},{'scalar','binary'},mfilename,'plot_ROCcurve')) % Plot ROC curve?
addParameter(p, 'save_ROCcurve',         1, @(x) v(x,{'numeric','logical'},{'scalar','binary'},mfilename,'save_ROCcurve')) % Save ROC curve?
addParameter(p, 'save_rates',            1, @(x) v(x,{'numeric','logical'},{'scalar','binary'},mfilename,'save_rates')) % Save hit and false alarm (FA) rates?
addParameter(p, 'out_dir',  'Outputs/ROC/', @(x) v(x,{'char'},{'nonempty'},mfilename,'out_directory')); % Directory to save rates if save_rates == 1
addParameter(p, 'fn_alpha',             [], @(x) v(x,{'char'},{'nonempty'},mfilename,'fn_alpha')); % Filename for FPR
addParameter(p, 'fn_beta',              [], @(x) v(x,{'char'},{'nonempty'},mfilename,'fn_beta')); % Filename for TPR
addParameter(p, 'theta_thresholds',     [], @(x) v(x,{'numeric'},{'vector'},mfilename,'theta_thresholds'));
addParameter(p, 'Rruns',               100, @(x) v(x,{'numeric'},{'scalar','integer','positive'},mfilename,'Rruns')); % R different soma arrangements, i.e. retinas
addParameter(p, 'Ktrials',               3, @(x) v(x,{'numeric'},{'scalar','integer','positive'},mfilename,'Ktrials')); % K trials or repetitions for each soma arrangement
addParameter(p, 'nrows',                10, @(x) v(x,{'numeric'},{'scalar','integer','positive'},mfilename,'nrows')); % n 1D arrays
addParameter(p, 'mcells',              100, @(x) v(x,{'numeric'},{'scalar','integer','positive'},mfilename,'mcells')); % m cells per array for each subtype
addParameter(p, 'dt',                0.025, @(x) v(x,{'numeric'},{'scalar','positive'},mfilename,'dt')); % seconds
addParameter(p, 'prefSide_width',      110, @(x) v(x,{'numeric'},{'scalar','positive'},mfilename,'prefSide_width')); % Width of region in which cell is sensitive to emerging edge moving in ND (microns)
addParameter(p, 'somaDistance_mean_PD', 39, @(x) v(x,{'numeric'},{'scalar','positive'},mfilename,'somaDistance_mean_PD')); % Inter-soma distance mean for PD-preferring (microns)
addParameter(p, 'somaDistance_mean_ND', 39, @(x) v(x,{'numeric'},{'scalar','positive'},mfilename,'somaDistance_mean_ND')); % Inter-soma distance mean for ND-preferring (microns)
addParameter(p, 'somaDistance_sd',      16, @(x) v(x,{'numeric'},{'scalar','positive'},mfilename,'somaDistance_sd')); % Inter-soma distance standard deviation (microns)
addParameter(p, 'radius_movingField',  330, @(x) v(x,{'numeric'},{'scalar','positive'},mfilename,'radius_movingField')); % Moving field radius (microns)
addParameter(p, 'alignPD',               0, @(x) v(x,{'numeric','logical'},{'scalar','binary'},mfilename,'alignPD')); % Align preferred direction (PD) axes of all DSGCs?
addParameter(p, 'alignMargin',        14.1, @(x) v(x,{'numeric'},{'scalar','positive'},mfilename,'alignMargin')); % PD +/- ___ degrees. Taken from uniform distribution
addParameter(p, 'spikeRespVar',        0.4, @(x) v(x,{'numeric'},{'scalar','positive'},mfilename,'spikeRespVar')); % Sub-Poisson noise model variance
addParameter(p, 'filter_center_A',     120, @(x) v(x,{'numeric'},{'scalar'},mfilename,'filter_center_A')); % Height of center Gaussian for temporal filter
addParameter(p, 'filter_surround_A',    50, @(x) v(x,{'numeric'},{'scalar'},mfilename,'filter_center_A')); % Height of surround Gaussian for temporal filter
addParameter(p, 'filter_center_sd',    100, @(x) v(x,{'numeric'},{'scalar','>=',0},mfilename,'filter_center_sd')); % Width of center Gaussian for temporal filter (milliseconds)
addParameter(p, 'filter_surround_sd',  250, @(x) v(x,{'numeric'},{'scalar','>=',0},mfilename,'filter_surround_sd')); % Width of surround Gaussian for temporal filter (milliseconds)
addParameter(p, 'filter_t0',           210, @(x) v(x,{'numeric'},{'scalar','>=',0},mfilename,'filter_t0')); % Mean of Gaussians (milliseconds)
addParameter(p, 'filter_tLength',     1000, @(x) v(x,{'numeric'},{'scalar','positive'},mfilename,'filter_tLength')); % Temporal length of filter (milliseconds)
addParameter(p, 'occ_window',         -2:2, @(x) v(x,{'numeric'},{'vector','integer'},mfilename,'occ_window')); % Time window around occ_tbin defining occlusion zone
addParameter(p, 'noise_only',            0, @(x) v(x,{'numeric','logical'},{'scalar','binary'},mfilename,'noise_only')); % Simulate noise only?

parse(p, varargin{:});

% Write results from input parser.
sampling_dir = p.Results.sampling_dir;
occ_xpos = p.Results.occ_xpos;
baselineFR = p.Results.baselineFR;
barSpeed = p.Results.barSpeed;
plot_ROCcurve = p.Results.plot_ROCcurve;
save_ROCcurve = p.Results.save_ROCcurve;
save_rates = p.Results.save_rates;
out_dir = p.Results.out_dir;
fn_alpha = p.Results.fn_alpha;
fn_beta = p.Results.fn_beta;
theta_thresholds = p.Results.theta_thresholds;
Rruns = p.Results.Rruns;
Ktrials = p.Results.Ktrials;
nrows = p.Results.nrows;
mcells = p.Results.mcells;
dt = p.Results.dt;
prefSide_width = p.Results.prefSide_width;
somaDistance_mean_PD = p.Results.somaDistance_mean_PD;
somaDistance_mean_ND = p.Results.somaDistance_mean_ND;
somaDistance_sd = p.Results.somaDistance_sd;
radius_movingField = p.Results.radius_movingField;
alignPD = p.Results.alignPD;
alignMargin = p.Results.alignMargin;
spikeRespVar = p.Results.spikeRespVar;
filter_center_A = p.Results.filter_center_A;
filter_surround_A = p.Results.filter_surround_A;
filter_center_sd = p.Results.filter_center_sd;
filter_surround_sd = p.Results.filter_surround_sd;
filter_t0 = p.Results.filter_t0;
filter_tLength = p.Results.filter_tLength;
occ_window = p.Results.occ_window;
noise_only = p.Results.noise_only;

clearvars p varargin v

if ~plot_ROCcurve
    save_ROCcurve = 0;
end

%% Create dependent parameters.
if save_rates
    if isempty(fn_alpha)
        fn_alpha = strcat('FPR_bkgdFR',num2str(baselineFR*1000,'%.6i'),'_windowpm',num2str(occ_window(end)),'_noiseonly',num2str(noise_only));
    end
    if isempty(fn_beta)
        fn_beta = strcat('TPR_bkgdFR',num2str(baselineFR*1000,'%.6i'),'_windowpm',num2str(occ_window(end)),'_noiseonly',num2str(noise_only));
    end
end
    
if isempty(theta_thresholds)
    if baselineFR <= 1
        theta_thresholds = 5e5:2e4:2e6;
    elseif baselineFR >= 100
        theta_thresholds = 1e7:1e5:1.5e7;
    else
        theta_thresholds = 5e5:1e4:17e5;
    end
end

occ_tbin = ceil(occ_xpos/barSpeed/dt);

spikeResp_sd = sqrt(spikeRespVar); % Sub-Poisson noise model standard deviation

%% Create temporal filter.
[filter, filterLength] = makeTempFilter_diffGaussians(filter_center_A, filter_center_sd, filter_surround_A, filter_surround_sd, filter_t0, dt*1000, filter_tLength);
clearvars filter_center_A filter_center_sd filter_surround_A filter_surround_sd filter_t0 filter_tLength

%% Load sampling data.
runModel_dir = pwd;
cd(sampling_dir)
% Load spike response onset times
load('ff_cells_onset.mat','ff_cells_onset')
% load('occ_cells_onset.mat','occ_cells_onset')
load('onr_cells_onset.mat','onr_cells_onset')
% Load spike response onset-duration correlation
load('ff_onset_duration_corr.mat','ff_onset_duration_corr')
% load('occ_onset_duration_corr.mat','occ_onset_duration_corr')
load('onr_onset_duration_corr.mat','onr_onset_duration_corr')
% Load spike response peak firing rates
load('ff_peakFR_PD.mat','ff_peakFR_PD')
% load('occ_peakFR_PD.mat','occ_peakFR_PD')
load('onr_peakFR_ND.mat','onr_peakFR_ND')
% Load direction tuning curve parameters
load('tuningWidth_gaussFit_param.mat','tuningWidth_gaussFit_param')
load('gaussFit_collapsedTuningCurve_SD.mat','gaussFit_collapsedTuningCurve_SD')
cd(runModel_dir)
clearvars runModel_dir sampling_dir

%% Run model
falsePositiveRate = NaN(1,length(theta_thresholds)); % Initialize false alarm rates
truePositiveRate = NaN(1,length(theta_thresholds)); % Initialize hit rates

ithreshold_truePositiveRate = NaN(length(theta_thresholds),Rruns,Ktrials);
ithreshold_falsePositiveRate = NaN(length(theta_thresholds),Rruns,Ktrials);

for ithreshold = 1:length(theta_thresholds)
    threshold = theta_thresholds(ithreshold);
    irun_TPR = zeros(Rruns,Ktrials);
    irun_FPR = zeros(Rruns,Ktrials);
    
    for irun = 1:Rruns
        [TPR, FPR] = model_coincidenceDetection(threshold, filter, filterLength, nrows, mcells, somaDistance_mean_PD, somaDistance_mean_ND, somaDistance_sd, radius_movingField, dt, barSpeed, Ktrials, baselineFR, alignPD, alignMargin, spikeResp_sd, occ_xpos, occ_tbin, prefSide_width, occ_window, noise_only, ff_cells_onset, onr_cells_onset, ff_onset_duration_corr, onr_onset_duration_corr, ff_peakFR_PD, onr_peakFR_ND, tuningWidth_gaussFit_param, gaussFit_collapsedTuningCurve_SD);
        irun_TPR(irun,:) = TPR;
        irun_FPR(irun,:) = FPR;
        disp(['Threshold ' num2str(ithreshold) '/' num2str(length(theta_thresholds)) ' Run ' num2str(irun) '/' num2str(Rruns)])
    end
    
    ithreshold_truePositiveRate(ithreshold,:,:) = irun_TPR;
    ithreshold_falsePositiveRate(ithreshold,:,:) = irun_FPR;
    
    falsePositiveRate(ithreshold) = mean2(ithreshold_falsePositiveRate(ithreshold,:,:));
    truePositiveRate(ithreshold) = mean2(ithreshold_truePositiveRate(ithreshold,:,:));
    
end

%% Visualize ROC curve and save hit and FA rates
if plot_ROCcurve
    figure
    scatter(falsePositiveRate,truePositiveRate,10,'filled')
    xlabel('False positive rate'); ylabel('Hit rate')
    xlim([-0.1 1.1]); ylim([-0.1 1.1])
    if save_ROCcurve
        savefig(strcat(out_dir,'ROC_bkgdFR',num2str(baselineFR*1000,'%.6i'),'_windowpm',num2str(occ_window(end)),'_noiseonly',num2str(noise_only)));
    end
end

if save_rates
    save(strcat(out_dir,fn_alpha),'falsePositiveRate')
    save(strcat(out_dir,fn_beta),'truePositiveRate')
end

end


