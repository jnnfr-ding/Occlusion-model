function RFwidths = receptiveFieldWidth(Ncells,...
    radius_dendriticField_mean,radius_dendriticField_SD,dendriticField2RF)
%{
INPUTS:
    Ncells                      = 
    radius_dendriticField_mean  =
    radius_dendriticField_SD    =
    dendriticField2RF           =
OUTPUTS:
    RFwidths    = Nx1
%}

radius_dendriticField = normrnd(radius_dendriticField_mean,...
    radius_dendriticField_SD,[Ncells,1]);   % SD for Gaussian RF
radius_dendriticField(radius_dendriticField<0) = ...
    radius_dendriticField_mean;             % Make sure SD nonnegative

RFwidths = dendriticField2RF*radius_dendriticField;

end