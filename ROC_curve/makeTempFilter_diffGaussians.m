function [filter, N] = makeTempFilter_diffGaussians(...
    B_c,sd_c,B_s,sd_s,t_0,dt,tLength)
%{
Make... temporal filter... difference of Gaussians
%}
% INPUTS:
%   B_c     = 1x1 Height of center Gaussian (msec)
%   sd_c    = 1x1 Width of center Gaussian (msec)
%   B_s     = 1x1 Height of surround Gaussian (msec)
%   sd_s    = 1x1 Width of surround Gaussian (msec)
%   t_0     = 1x1 Mean of Gaussians (msec)
%   dt      = 1x1 Size of time step or bin size (msec)
%   tLength = 1x1 Temporal length of filter (msec)
% OUTPUTS:
%   filter  = 1xN Difference of Gaussians temporal filter
%   N       = 1x1 Number of time bins in temporal filter

tbins = 0:dt:tLength-dt;
N = length(tbins);

filter = B_c*exp(-(tbins-t_0).^2./(2*sd_c^2))...
    - B_s*exp(-(tbins-t_0).^2./(2*sd_s^2));

end





