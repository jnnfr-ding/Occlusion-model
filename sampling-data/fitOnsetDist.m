% Run onset_distribution.m first
onset_distribution

%% Occlusion PD and PD-adjacent directions
numData = 100000;    % Simulate 100k data points
muTemp = nanmean(ff_cells_onset,'all');
sigmaTemp = nanstd(ff_cells_onset,0,'all');
skewTemp = skewness(ff_cells_onset,0,'all');
kurtTemp = kurtosis(ff_cells_onset,0,'all');

%%
num = 1;
type = NaN(num,1);
for i = 1:num
    [dataTemp type(i)]= pearsrnd(muTemp,sigmaTemp,skewTemp,kurtTemp,numData,1);
    y = histc(dataTemp,edges);
    y = y/sum(y);   % Normalize

    figure
    hold on
    %histogram(occ_cells_onset(:),edges,'FaceColor',[0 0.4470 0.7410])
    plot(edges,ff_onset_distribution_PD_PDadj./sum(ff_onset_distribution_PD_PDadj),...
        'Color',[0 0.4470 0.7410],'LineWidth',1.2)
    %histogram(dataTemp,edges,'FaceColor',[0.9290 0.6940 0.1250])
    plot(edges,y,'Color',[0.9290 0.6940 0.1250],'LineWidth',1.2) % Plot PDF
    xlabel('Spike response onset time (sec)')
    ylabel('Probability')
    title('Beta distribution fit to occlusion null onset times')
    legend('PDF of data','Fit','Location','best')
end

