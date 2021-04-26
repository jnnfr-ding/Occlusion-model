%% Parameters
dt = 0.010;         % Time bin size (sec)
t_axis = 0:dt:50;   % Time axis (sec)
numOccPos = 4;
barSpeeds = [330,660,1320,1980,2640];

%% Load data
% Change directory to where .mat data files are stored
directory = 'C:\Users\jnnfr\Dropbox\Occlusion Project\occlusion_model_jan2021\20210114';
cd(directory)
% Get a list of all .mat data files
filePattern = fullfile(directory,'*.mat');  % File pattern ending in .mat
fileList = dir(filePattern);                % List of files in folder
for file = 1:length(fileList)
    fn = fileList(file).name;               % File name
    load(fn)                                % Load data file
end

% %% Baseline FR = 0 spikes/s; occluder position 1100 microns
% left_x = 1000;  % Left spatial boundary (microns)
% right_x = 4000; % Right spatial boundary (microns)
% % Concatenate full-field data for NO baseline firing
% ff_allspeeds_bkgdFR0000 = [ff_speed0330_bkgdFR0000; ff_speed0660_bkgdFR0000;...
%     ff_speed1320_bkgdFR0000; ff_speed1980_bkgdFR0000; ff_speed2640_bkgdFR0000];
% % Concatenate occlusion data for occluder position 1100 microns and NO baseline firing
% occ_allspeeds_bkgdFR0000_occPos1100 = [occ_speed0330_bkgdFR0000_occPos1100;...
%     occ_speed0660_bkgdFR0000_occPos1100; occ_speed1320_bkgdFR0000_occPos1100;...
%     occ_speed1980_bkgdFR0000_occPos1100; occ_speed2640_bkgdFR0000_occPos1100];
% 
% RMSEdecrease = NaN(2,length(barSpeeds));    % Initialize
%                                             % 1st row is % decrease
%                                             % 2nd row is raw RMSE decrease
% for jspeed = 1:length(barSpeeds)    % Loop through bar speeds
%     barSpeed = barSpeeds(jspeed);   % j-th bar speed
%     % Use only RMSE between spatial positions left_x and right_x
%     temp_ff = ff_allspeeds_bkgdFR0000(jspeed,:);
%     temp_ff = temp_ff(ceil(left_x/barSpeed/dt):ceil(right_x/barSpeed/dt));
%     temp_occ = occ_allspeeds_bkgdFR0000_occPos1100(jspeed,:);
%     temp_occ = temp_occ(ceil(left_x/barSpeed/dt):ceil(right_x/barSpeed/dt));
%     
%     RMSEdecrease(1,jspeed) = max((temp_ff-temp_occ)./temp_ff)*100;
%     [ampDecrease,idx] = max(temp_ff-temp_occ);
%     RMSEdecrease(2,jspeed) = ampDecrease;
%     checkOccPos = (idx + ceil(left_x/barSpeed/dt) - 1)*barSpeed*dt;
%     disp(checkOccPos)
% end
% 
% figure  % Percent decrease
% plot(barSpeeds,RMSEdecrease(1,:),'LineWidth',1.2)
% xlabel('Bar speed (micron/s)'); ylabel('% decrease')
% ylim([0 100])
% title('% decrease in position estimate RMSE at 0 spike/s baseline FR, occluder position 1100 microns')
% 
% figure  % Absolute RMSE decrease
% plot(barSpeeds,RMSEdecrease(2,:),'LineWidth',1.2)
% xlabel('Bar speed (micron/s)'); ylabel('RMSE decrease (deg)')
% ylim([0 5])
% title('Absolute amplitude decrease in position estimate RMSE at 0 spike/s baseline FR, occluder position 1100 microns')

%% Baseline FR = 0 spikes/s; occluder position 1800 microns
left_x = 1000;  % Left spatial boundary (microns)
right_x = 4000; % Right spatial boundary (microns)
% Concatenate full-field data for NO baseline firing
ff_allspeeds_bkgdFR0000 = [ff_speed0330_bkgdFR0000; ff_speed0660_bkgdFR0000;...
    ff_speed1320_bkgdFR0000; ff_speed1980_bkgdFR0000; ff_speed2640_bkgdFR0000];
% Concatenate occlusion data for occluder position 1100 microns and NO baseline firing
occ_allspeeds_bkgdFR0000_occPos1800 = [occ_speed0330_bkgdFR0000_occPos1800;...
    occ_speed0660_bkgdFR0000_occPos1800; occ_speed1320_bkgdFR0000_occPos1800;...
    occ_speed1980_bkgdFR0000_occPos1800; occ_speed2640_bkgdFR0000_occPos1800];

RMSEdecrease_baselineFR0000 = NaN(2,length(barSpeeds));    % Initialize
                                            % 1st row is % decrease
                                            % 2nd row is raw RMSE decrease
for jspeed = 1:length(barSpeeds)    % Loop through bar speeds
    barSpeed = barSpeeds(jspeed);   % j-th bar speed
    % Use only RMSE between spatial positions left_x and right_x
    temp_ff = ff_allspeeds_bkgdFR0000(jspeed,:);
    temp_ff = temp_ff(ceil(left_x/barSpeed/dt):ceil(right_x/barSpeed/dt));
    temp_occ = occ_allspeeds_bkgdFR0000_occPos1800(jspeed,:);
    temp_occ = temp_occ(ceil(left_x/barSpeed/dt):ceil(right_x/barSpeed/dt));
    
    RMSEdecrease_baselineFR0000(1,jspeed) = max((temp_ff-temp_occ)./temp_ff)*100;
    [ampDecrease,idx] = max(temp_ff-temp_occ);
    RMSEdecrease_baselineFR0000(2,jspeed) = ampDecrease;
    checkOccPos = (idx + ceil(left_x/barSpeed/dt) - 1)*barSpeed*dt;
    disp(checkOccPos)
end

figure  % Percent decrease
plot(barSpeeds,RMSEdecrease_baselineFR0000(1,:),'LineWidth',1.2)
xlabel('Bar speed (micron/s)'); ylabel('% decrease')
ylim([0 100])
title('% decrease in position estimate RMSE at 0 spike/s baseline FR, occluder position 1800 microns')

figure  % Absolute RMSE decrease
plot(barSpeeds,RMSEdecrease_baselineFR0000(2,:),'LineWidth',1.2)
xlabel('Bar speed (micron/s)'); ylabel('RMSE decrease (deg)')
ylim([0 5])
title('Absolute amplitude decrease in position estimate RMSE at 0 spike/s baseline FR, occluder position 1800 microns')

%% Baseline FR = 0.025 spikes/s; occluder position 1800 microns
left_x = 1000;  % Left spatial boundary (microns)
right_x = 4000; % Right spatial boundary (microns)
% Concatenate full-field data for NO baseline firing
ff_allspeeds_bkgdFR0025 = [ff_speed0330_bkgdFR0025; ff_speed0660_bkgdFR0025;...
    ff_speed1320_bkgdFR0025; ff_speed1980_bkgdFR0025; ff_speed2640_bkgdFR0025];
% Concatenate occlusion data for occluder position 1100 microns and NO baseline firing
occ_allspeeds_bkgdFR0025_occPos1800 = [occ_speed0330_bkgdFR0025_occPos1800;...
    occ_speed0660_bkgdFR0025_occPos1800; occ_speed1320_bkgdFR0025_occPos1800;...
    occ_speed1980_bkgdFR0025_occPos1800; occ_speed2640_bkgdFR0025_occPos1800];

RMSEdecrease_baselineFR0025 = NaN(2,length(barSpeeds));    % Initialize
                                            % 1st row is % decrease
                                            % 2nd row is raw RMSE decrease
for jspeed = 1:length(barSpeeds)    % Loop through bar speeds
    barSpeed = barSpeeds(jspeed);   % j-th bar speed
    % Use only RMSE between spatial positions left_x and right_x
    temp_ff = ff_allspeeds_bkgdFR0025(jspeed,:);
    temp_ff = temp_ff(ceil(left_x/barSpeed/dt):ceil(right_x/barSpeed/dt));
    temp_occ = occ_allspeeds_bkgdFR0025_occPos1800(jspeed,:);
    temp_occ = temp_occ(ceil(left_x/barSpeed/dt):ceil(right_x/barSpeed/dt));
    
    RMSEdecrease_baselineFR0025(1,jspeed) = max((temp_ff-temp_occ)./temp_ff)*100;
    [ampDecrease,idx] = max(temp_ff-temp_occ);
    RMSEdecrease_baselineFR0025(2,jspeed) = ampDecrease;
    checkOccPos = (idx + ceil(left_x/barSpeed/dt) - 1)*barSpeed*dt;
    disp(checkOccPos)
end

figure  % Percent decrease
plot(barSpeeds,RMSEdecrease_baselineFR0025(1,:),'LineWidth',1.2)
xlabel('Bar speed (micron/s)'); ylabel('% decrease')
ylim([0 100])
title('% decrease in position estimate RMSE at 0.025 spike/s baseline FR, occluder position 1800 microns')

figure  % Absolute RMSE decrease
plot(barSpeeds,RMSEdecrease_baselineFR0025(2,:),'LineWidth',1.2)
xlabel('Bar speed (micron/s)'); ylabel('RMSE decrease (deg)')
ylim([0 5])
title('Absolute amplitude decrease in position estimate RMSE at 0.025 spike/s baseline FR, occluder position 1800 microns')

%% Baseline FR = 0.1 spikes/s; occluder position 1800 microns
left_x = 1000;  % Left spatial boundary (microns)
right_x = 4000; % Right spatial boundary (microns)
% Concatenate full-field data for NO baseline firing
ff_allspeeds_bkgdFR0100 = [ff_speed0330_bkgdFR0100; ff_speed0660_bkgdFR0100;...
    ff_speed1320_bkgdFR0100; ff_speed1980_bkgdFR0100; ff_speed2640_bkgdFR0100];
% Concatenate occlusion data for occluder position 1100 microns and NO baseline firing
occ_allspeeds_bkgdFR0100_occPos1800 = [occ_speed0330_bkgdFR0100_occPos1800;...
    occ_speed0660_bkgdFR0100_occPos1800; occ_speed1320_bkgdFR0100_occPos1800;...
    occ_speed1980_bkgdFR0100_occPos1800; occ_speed2640_bkgdFR0100_occPos1800];

RMSEdecrease_baselineFR0100 = NaN(2,length(barSpeeds));    % Initialize
                                            % 1st row is % decrease
                                            % 2nd row is raw RMSE decrease
for jspeed = 1:length(barSpeeds)    % Loop through bar speeds
    barSpeed = barSpeeds(jspeed);   % j-th bar speed
    % Use only RMSE between spatial positions left_x and right_x
    temp_ff = ff_allspeeds_bkgdFR0100(jspeed,:);
    temp_ff = temp_ff(ceil(left_x/barSpeed/dt):ceil(right_x/barSpeed/dt));
    temp_occ = occ_allspeeds_bkgdFR0100_occPos1800(jspeed,:);
    temp_occ = temp_occ(ceil(left_x/barSpeed/dt):ceil(right_x/barSpeed/dt));
    
    RMSEdecrease_baselineFR0100(1,jspeed) = max((temp_ff-temp_occ)./temp_ff)*100;
    [ampDecrease,idx] = max(temp_ff-temp_occ);
    RMSEdecrease_baselineFR0100(2,jspeed) = ampDecrease;
    checkOccPos = (idx + ceil(left_x/barSpeed/dt) - 1)*barSpeed*dt;
    disp(checkOccPos)
end

figure  % Percent decrease
plot(barSpeeds,RMSEdecrease_baselineFR0100(1,:),'LineWidth',1.2)
xlabel('Bar speed (micron/s)'); ylabel('% decrease')
ylim([0 100])
title('% decrease in position estimate RMSE at 0.1 spike/s baseline FR, occluder position 1800 microns')

figure  % Absolute RMSE decrease
plot(barSpeeds,RMSEdecrease_baselineFR0100(2,:),'LineWidth',1.2)
xlabel('Bar speed (micron/s)'); ylabel('RMSE decrease (deg)')
ylim([0 5])
title('Absolute amplitude decrease in position estimate RMSE at 0.1 spike/s baseline FR, occluder position 1800 microns')

%%
figure; hold on     % Percent decrease
plot(barSpeeds,RMSEdecrease_baselineFR0000(1,:),'LineWidth',1.2,...
    'DisplayName','Baseline FR = 0 spike/s')
plot(barSpeeds,RMSEdecrease_baselineFR0025(1,:),'LineWidth',1.2,...
    'DisplayName','Baseline FR = 0.025 spike/s')
plot(barSpeeds,RMSEdecrease_baselineFR0100(1,:),'LineWidth',1.2,...
    'DisplayName','Baseline FR = 0.1 spike/s')
xlabel('Bar speed'); ylabel('% decrease')
ylim([0 100])
legend

figure; hold on     % Absolute RMSE decrease
plot(barSpeeds,RMSEdecrease_baselineFR0000(2,:),'LineWidth',1.2,...
    'DisplayName','Baseline FR = 0 spike/s')
plot(barSpeeds,RMSEdecrease_baselineFR0025(2,:),'LineWidth',1.2,...
    'DisplayName','Baseline FR = 0.025 spike/s')
plot(barSpeeds,RMSEdecrease_baselineFR0100(2,:),'LineWidth',1.2,...
    'DisplayName','Baseline FR = 0.1 spike/s')
xlabel('Bar speed'); ylabel('Absolute RMSE decrease (deg)')
ylim([0 5])
legend

%% Baseline FR = 0 spikes/s; occluder position 1800 microns
left_x = 1000;  % Left spatial boundary (microns)
right_x = 4000; % Right spatial boundary (microns)
% Concatenate full-field data for NO baseline firing
ff_allspeeds_bkgdFR0000 = [ff_speed0330_bkgdFR0000; ff_speed0660_bkgdFR0000;...
    ff_speed1320_bkgdFR0000; ff_speed1980_bkgdFR0000; ff_speed2640_bkgdFR0000];
% Concatenate occlusion data for occluder position 1100 microns and NO baseline firing
occ_allspeeds_bkgdFR0000_occPos1800 = [occ_speed0330_bkgdFR0000_occPos1800;...
    occ_speed0660_bkgdFR0000_occPos1800; occ_speed1320_bkgdFR0000_occPos1800;...
    occ_speed1980_bkgdFR0000_occPos1800; occ_speed2640_bkgdFR0000_occPos1800];


RMSEdecrease_baselineFR0000_ff_occ = NaN(2,length(barSpeeds));    % Initialize
                                            % 1st row is full-field
                                            % 2nd row is occlusion
for jspeed = 1:length(barSpeeds)    % Loop through bar speeds
    barSpeed = barSpeeds(jspeed);   % j-th bar speed
    % Use only RMSE between spatial positions left_x and right_x
    temp_ff = ff_allspeeds_bkgdFR0000(jspeed,:);
    temp_ff = temp_ff(ceil(left_x/barSpeed/dt):ceil(right_x/barSpeed/dt));
    temp_occ = occ_allspeeds_bkgdFR0000_occPos1800(jspeed,:);
    temp_occ = temp_occ(ceil(left_x/barSpeed/dt):ceil(right_x/barSpeed/dt));
    
    [ampDecrease,idx] = max(temp_ff-temp_occ);
    % Full-field
    RMSEdecrease_baselineFR0000_ff_occ(1,jspeed) = temp_ff(idx);
    % Occlusion
    RMSEdecrease_baselineFR0000_ff_occ(2,jspeed) = temp_occ(idx);
    checkOccPos = (idx + ceil(left_x/barSpeed/dt) - 1)*barSpeed*dt;
    disp(checkOccPos)    
end

figure; hold on
plot(barSpeeds,RMSEdecrease_baselineFR0000_ff_occ(1,:),...
    'DisplayName','Full-field','LineWidth',1.2)
plot(barSpeeds,RMSEdecrease_baselineFR0000_ff_occ(2,:),...
    'DisplayName','Occlusion','LineWidth',1.2)
xlabel('Bar speed (micron/s)'); ylabel('RMSE at occlusion (deg)')
legend


%% Baseline FR = 0.025 spikes/s; occluder position 1800 microns
left_x = 1000;  % Left spatial boundary (microns)
right_x = 4000; % Right spatial boundary (microns)
% Concatenate full-field data for NO baseline firing
ff_allspeeds_bkgdFR0025 = [ff_speed0330_bkgdFR0025; ff_speed0660_bkgdFR0025;...
    ff_speed1320_bkgdFR0025; ff_speed1980_bkgdFR0025; ff_speed2640_bkgdFR0025];
% Concatenate occlusion data for occluder position 1100 microns and NO baseline firing
occ_allspeeds_bkgdFR0025_occPos1800 = [occ_speed0330_bkgdFR0025_occPos1800;...
    occ_speed0660_bkgdFR0025_occPos1800; occ_speed1320_bkgdFR0025_occPos1800;...
    occ_speed1980_bkgdFR0025_occPos1800; occ_speed2640_bkgdFR0025_occPos1800];


RMSEdecrease_baselineFR0025_ff_occ = NaN(2,length(barSpeeds));    % Initialize
                                            % 1st row is full-field
                                            % 2nd row is occlusion
for jspeed = 1:length(barSpeeds)    % Loop through bar speeds
    barSpeed = barSpeeds(jspeed);   % j-th bar speed
    % Use only RMSE between spatial positions left_x and right_x
    temp_ff = ff_allspeeds_bkgdFR0025(jspeed,:);
    temp_ff = temp_ff(ceil(left_x/barSpeed/dt):ceil(right_x/barSpeed/dt));
    temp_occ = occ_allspeeds_bkgdFR0025_occPos1800(jspeed,:);
    temp_occ = temp_occ(ceil(left_x/barSpeed/dt):ceil(right_x/barSpeed/dt));
    
    [ampDecrease,idx] = max(temp_ff-temp_occ);
    % Full-field
    RMSEdecrease_baselineFR0025_ff_occ(1,jspeed) = temp_ff(idx);
    % Occlusion
    RMSEdecrease_baselineFR0025_ff_occ(2,jspeed) = temp_occ(idx);
    checkOccPos = (idx + ceil(left_x/barSpeed/dt) - 1)*barSpeed*dt;
    disp(checkOccPos)    
end

figure; hold on
plot(barSpeeds,RMSEdecrease_baselineFR0025_ff_occ(1,:),...
    'DisplayName','Full-field','LineWidth',1.2)
plot(barSpeeds,RMSEdecrease_baselineFR0025_ff_occ(2,:),...
    'DisplayName','Occlusion','LineWidth',1.2)
xlabel('Bar speed (micron/s)'); ylabel('RMSE at occlusion (deg)')
title('Baseline FR = 0.025 spikes/s')
legend