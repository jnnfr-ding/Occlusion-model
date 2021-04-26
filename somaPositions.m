function pos = somaPositions(...
    nrows,mcells,meanDist,sigmaDist,bufferDist)
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

% nxm Soma positions with Gaussian distributed inter-soma distances
pos = sigmaDist/sqrt(2).*randn(nrows,mcells) + bufferDist + ...
    repmat([0:meanDist:meanDist*(mcells-1)],nrows,1);

% nxm Arrange soma positions in ascending order
pos = sort(pos,2);  % Sort elements of each row


end