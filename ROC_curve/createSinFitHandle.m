function y = createSinFitHandle(p)
%{
Create function handle for sine fit to spike response. Amplitude is given 
by peak firing rate; period is given by twice the spike response duration; 
phase is given by spike response onset time.
%}
% INPUTS:
%   p(1) = amplitude
%   p(2) = period
%   p(3) = phase
% OUTPUT:
%   y    = function handle for sine fit

y = @(t) p(1).*sin(2.*pi./p(2).*(t - p(3)));

end