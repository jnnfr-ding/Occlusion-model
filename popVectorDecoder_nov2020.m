function [estimates,labels_noisy] = popVectorDecoder_nov2020(...
    DSGCpopResp,labels,RFwidths)
%{
Compute (position) estimate using a population vector, which is the
firing-rate weighted average of the neurons' preferred values or assigned
labels. N.B. The population vector is the optimal max likelihood (ML)
decoder for independent neurons with Poisson firing and Gaussian tuning
curves with uniform widths.

INPUTS:
    DSGCpopResp = 2NxTxK    Array of firing-rate responses
                            N cells
                            T time bins
                            K trials
    labels      = 2Nx1      Labels
    RFwidths    = 2Nx1      Receptive field (RF) widths
OUTPUTS:
    estimates       =  KxT  Estimates computed by decoder
    labels_noisy    = 2Nx1  Noisy labels used by decoder    
%}

[Ncells,Tbins,Ktrials] = size(DSGCpopResp); % Get dimensions of response array
Ncells = Ncells/2;
estimates = zeros(Ktrials,Tbins);           % Initialize array of estimates
labels_noisy = normrnd(labels,RFwidths);    % Add RF noise to (position) labels
for tbin = 1:Tbins                          % Loop through T time bins
    % Estimate at t-th time bin
    % Firing-rate weighted average of labels
    estimates(:,tbin) = (squeeze(DSGCpopResp(:,tbin,:))'*(labels_noisy./(RFwidths.^2)))...
        ./sum(squeeze(DSGCpopResp(:,tbin,:))'./repmat((RFwidths.^2)',[Ktrials,1]),2);
        %*% Numerator
        % Weight by firing rate (K x 2N)
        % Soma positions with noise (2N x 1)
        % Divide by RF widths squared (2N x 1)
        %*% Denominator
        % Divide response by RF width squared for all K trials (K x 2N)
        % Responses summed over 2N cells (K x 1)
end


end