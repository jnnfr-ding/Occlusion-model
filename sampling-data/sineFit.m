function  y = sineFit(x)
%{
Compute sine fit to instantaneous firing rate given by PSTH. Amplitude is
given by peak firing rate; period is given by twice the spike response
duration; phase is given by spike response onset time.
%}
% INPUTS
%   x(1) = amplitude
%   x(2) = period
%   x(3) = phase
% OUTPUT
%   y    = function handle for sine fit

y = @(t) x(1).*sin(2.*pi./x(2).*(t - x(3)));

end