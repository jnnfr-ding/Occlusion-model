function spatialOnset = speed_adjustOnset(onr,spikeRespOnset,...
    barSpeed,barSpeed_baseline,slope,radius_movingField)
%{
INPUTS:
    spikeRespOnset      = Nx1   Onset times for bar speed 330 micron/sec
    barspeed            = 1x1   Current bar speed (micron/sec)
    barspeed_baseline   = 1x1   Baseline bar speed (i.e. 330 micron/sec)
    slope               = 1x1   Slope (sec) of spatial onset vs. bar speed
                                0.095 sec for full-field
                                0.083 sec for occlusion null response                     
OUTPUTS:
    spatialOnset        = Nx1   Spatial onset accounting for overall (both
                                speed-dependent and retinal processing, i.e.
                                70-ms) lag
%}
% Spatial onset at baseline bar speed(330 micron/s)
if onr
    % Occlusion null response
    spatialOnset_baseline = spikeRespOnset*barSpeed_baseline-radius_movingField;
else
    % Full-field
    spatialOnset_baseline = spikeRespOnset*barSpeed_baseline;
end

% Spatial onset at current bar speed
spatialOnset = spatialOnset_baseline +...
    slope*(barSpeed-barSpeed_baseline);

end