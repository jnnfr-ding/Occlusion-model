nrows = 10;             % n 1D arrays
mcells = 100;           % m cells per array for each subtype
Ncells = nrows*mcells;

somaDist_mean_PD = 39;  % Inter-soma distance mean for PD-preferring (microns)
somaDist_mean_ND = 39;  % Inter-soma distance mean for ND-preferring (microns)
somaDist_SD = 16;       % Inter-soma distance standard deviation (microns)
radius_soma = 10;          % DSGC soma radius (microns)
radius_RF = 110;            % DSGC RF radius (microns)
radius_movingField = 330;   % Moving field radius (microns)

[somaPos_x_rightPop,somaPos_y_rightPop,...
    somaPos_x_leftPop,somaPos_y_leftPop] = ...
    plotSomaPositions(nrows,mcells,somaDist_mean_PD,somaDist_mean_ND,...
    somaDist_SD,radius_movingField);

% Visualization
figure; hold on
% RF
viscircles([somaPos_x_leftPop,somaPos_y_leftPop],...
    repmat(radius_RF,[Ncells,1]),'Color',[0.75 0.75 0.75],'LineWidth',0.1,...
    'LineStyle','--')
viscircles([somaPos_x_rightPop,somaPos_y_rightPop],...
    repmat(radius_RF,[Ncells,1]),'Color',[1 0.75 0.75],'LineWidth',0.1,...
    'LineStyle','--')
% Soma
viscircles([somaPos_x_leftPop,somaPos_y_leftPop],...
    repmat(radius_soma,[Ncells,1]),'Color','k','LineWidth',0.8)
viscircles([somaPos_x_rightPop,somaPos_y_rightPop],...
    repmat(radius_soma,[Ncells,1]),'Color','r','LineWidth',0.8)
% scatter(somaPos_x_leftPop,somaPos_y_leftPop,10,'k','filled')
% scatter(somaPos_x_rightPop,somaPos_y_rightPop,10,'r','filled')
xlim([1000 2000])
% ylim([600 900])

function [somaPos_x_rightPop,somaPos_y_rightPop,...
    somaPos_x_leftPop,somaPos_y_leftPop] = ...
    plotSomaPositions(nrows,mcells,somaDist_mean_PD,somaDist_mean_ND,...
    somaDist_SD,radius_movingField)
%{
Generate soma positions with Gaussian inter-soma (nearest neighbor)
distance distribution

INPUTS:
    nrows       = n rows
    mcells      = m cells per row
    meanDist    = Mean inter-soma distance
    sigmaDist   = SD of inter-soma distance
    bufferDist  = Buffer distance from periphery of viewing field
OUTPUTS:
    pos         = n x m    Soma positions
%}
%*% SOMA POSITIONS
% Generate soma positions
somaPos_x_rightPop = somaPositions(...  % PD-preferring pop x-coordinates
    nrows,mcells,somaDist_mean_PD,somaDist_SD,2*radius_movingField);
somaPos_y_rightPop = somaPositions(...  % PD-preferring pop y-coordinates
    mcells,nrows,somaDist_mean_PD,somaDist_SD,2*radius_movingField)';
somaPos_x_leftPop = somaPositions(...   % ND-preferring pop x-coordinates
    nrows,mcells,somaDist_mean_ND,somaDist_SD,2*radius_movingField);
somaPos_y_leftPop = somaPositions(...   % ND-preferring pop y-coordinates
    mcells,nrows,somaDist_mean_ND,somaDist_SD,2*radius_movingField)';
% Reshape into single column vector (Nx1)
somaPos_x_rightPop = somaPos_x_rightPop';   
somaPos_x_rightPop = somaPos_x_rightPop(:); % PD-preferring pop x-coordinates
somaPos_y_rightPop = somaPos_y_rightPop';
somaPos_y_rightPop = somaPos_y_rightPop(:); % PD-preferring pop y-coordinates
somaPos_x_leftPop = somaPos_x_leftPop';
somaPos_x_leftPop = somaPos_x_leftPop(:);   % ND-preferring pop x-coordinates
somaPos_y_leftPop = somaPos_y_leftPop';
somaPos_y_leftPop = somaPos_y_leftPop(:);   % ND-preferring pop y-coordinates

end