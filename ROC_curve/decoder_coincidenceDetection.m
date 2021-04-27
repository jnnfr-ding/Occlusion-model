function [threshold_cross, conv_outputs] = decoder_coincidenceDetection(...
    resp,filter,filterLength,Tbins,threshold)
%{
%}
% INPUTS:
%   resp            = 1xT Firing-rate response (spikes/sec)
%   filter          = 1xL Temporal filter
%   filterLength    = 1x1 Length of temporal filter, L, (bins)
%   Tbins           = 1x1 Length of experimental trial, T, (bins)
%   threshold       = 1x1 Threshold for detection
% OUTPUTS:
%   threshold_cross = 1x(T-L) Above-threshold time bins
%   conv_outputs    = 1x(T-L) Output of convolving response with filter

% Initialize
threshold_cross = zeros(1,Tbins-filterLength);  % Above-threshold bins 
conv_outputs = NaN(1,Tbins-filterLength);       % Convolution outputs

for tbin = 1:Tbins-filterLength         % Loop through T time bins
    % Extract time bins that fall in temporal filter window
    % Convolve response with temporal filter   
    conv_output = sum(filter.*resp(tbin:tbin+filterLength-1));
    if conv_output > threshold
        % If convolution output is above threshold
        threshold_cross(tbin) = 1;
    end
    conv_outputs(tbin) = conv_output;   % Store convolution output
end

end



