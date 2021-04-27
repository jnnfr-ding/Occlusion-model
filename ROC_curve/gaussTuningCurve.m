function f = gaussTuningCurve(x,a,sigma)
%{
Find modulated/attenuated peak firing rates due to jittered PD alignment
from von Mises direction tuning curve approximated as Gaussian.
%}
% INPUTS:
%   x     = 1xN Jitters from preferred direction (deg)
%   a     = 1xN Peak firing rates, which are the amplitudes of the tuning
%               curves (Hz)
%   sigma = 1xN Standard deviations of the tuning curves (deg)
% OUTPUT:
%   f     = 1xN Firing rates at jittered directions (Hz)

f = a.*exp(-0.5*(x./sigma).^2);

end