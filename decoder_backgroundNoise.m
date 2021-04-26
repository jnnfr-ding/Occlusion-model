function [estimates,labels_noisy] = decoder_backgroundNoise(...
    DSGCpopResp,labels,RFwidths,barPos,decodeRange)
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
    barPos_x
    decodeRange = 
OUTPUTS:
    estimates       =  KxT  Estimates computed by decoder
    labels_noisy    = 2Nx1  Noisy labels used by decoder    
%}

[Ncells,Tbins,Ktrials] = size(DSGCpopResp); % Get dimensions of response array
Ncells = Ncells/2;
estimates = zeros(Ktrials,Tbins);           % Initialize array of estimates
labels_noisy = normrnd(labels,RFwidths);    % Add RF noise to (position) labels
for tbin = 1:Tbins                          % Loop through T time bins
    % Find cells within 400 microns of bar leading edge
    cell_idx = find(barPos(tbin)-decodeRange < labels &...
        labels < barPos(tbin)+decodeRange);
    
    % Estimate at t-th time bin
    % Firing-rate weighted average of labels
    if length(cell_idx) == 1
        estimates(:,tbin) = (squeeze(DSGCpopResp(cell_idx,tbin,:))'*(labels_noisy(cell_idx)./(RFwidths(cell_idx).^2)))'...
            ./sum(squeeze(DSGCpopResp(cell_idx,tbin,:))'./repmat((RFwidths(cell_idx).^2)',[Ktrials,1]),2);
        % Weight by firing rate (K x 2N)
        % Soma positions with noise (2N x 1)
        % Divide by responses summed over 2N cells (K x 1)
    else
        estimates(:,tbin) = (squeeze(DSGCpopResp(cell_idx,tbin,:))'*(labels_noisy(cell_idx)./(RFwidths(cell_idx).^2)))...
            ./sum(squeeze(DSGCpopResp(cell_idx,tbin,:))'./repmat((RFwidths(cell_idx).^2)',[Ktrials,1]),2);
        % Weight by firing rate (K x 2N)
        % Soma positions with noise (2N x 1)
        % Divide by responses summed over 2N cells (K x 1)
    end
    
end


end