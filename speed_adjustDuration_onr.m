function duration_speedAdjusted = speed_adjustDuration_onr(...
    spikeRespDuration,barSpeed,barSpeed_baseline,powerExp)
%{
INPUTS:
    spikeRespDuration   = Nx1
    barSpeed            = 1x1
    barSpeed_baseline   = 1x1
    powerExp            = 1x1 ... -0.561
OUTPUTS:
    duration_speedAdjusted  = Nx1
%}

r = barSpeed/barSpeed_baseline;
duration_speedAdjusted = spikeRespDuration.*power(r,powerExp);

end