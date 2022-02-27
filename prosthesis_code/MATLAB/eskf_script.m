clear
clc
close all
%% Load Prosthesis_data 
% ----------------------------------------------------------------------------------------------------------------------------------------------------------
% -- Download from : https://purdue0-my.sharepoint.com/:f:/g/personal/solimana_purdue_edu/Elv6C4eptuhOr7Q72VWuZWoBpcJ1ic5i27-tr9B9FnpJ1w?e=OdUiqj
load('/Users/ahmedsoliman/Documents/GitHub/os-eskf/Data/Prosthesis_data/cable_driven_prosthesis_020622.mat')


%% Preproc vars - sec 1 - after load

tr = trial(3);

rf1 = eskf_cdprosthesis.rf1;
rs2 = eskf_cdprosthesis.rs2;
mag_north = eskf_cdprosthesis.MAG_NORTH' / sqrt(sum(eskf_cdprosthesis.MAG_NORTH.^2));
Ta1 = eskf_cdprosthesis.Ta1; 
Tw1 = eskf_cdprosthesis.Tw1;
Ta2 = eskf_cdprosthesis.Ta2;
Tw2 = eskf_cdprosthesis.Tw2;
Am2 = eskf_cdprosthesis.Am2 / eskf_cdprosthesis.expmfs2;
Am1 = eskf_cdprosthesis.Am1 / eskf_cdprosthesis.expmfs1;
T1 = quat2rotm(rotm2quat(Ta1 + Tw1));
T2 = quat2rotm(rotm2quat(Ta2 + Tw2));
abias1 = eskf_cdprosthesis.abias1;
abias2 = eskf_cdprosthesis.abias2;
wbias1 = eskf_cdprosthesis.wbias1;
wbias2 = eskf_cdprosthesis.wbias2;
mbias1 = eskf_cdprosthesis.mbias1 + [0.573917728061859,-0.444902858259223,17.0118052633663]*1;
mbias2 = eskf_cdprosthesis.mbias2 + [6.55176218522236,-0.824504064840478,-1.76627658024604]*1;
Fs = 400;

t = (20: 1/Fs : tr.mtime(end))'; %% skip to start of expereiemnt 

% loads synchronized time - t_raw
load('/Users/ahmedsoliman/Documents/GitHub/os-eskf/Data/Prosthesis_data/t_raw.mat')

% ensures 400 Hz sampling rate , turns components into foot and shank
% frames (from OMC calibration)
% 1 - foot
% 2 - shank
am2 = interp1(t_raw, tr.packet.A1 - abias2, t) * inv(Ta2)';
wm2 = interp1(t_raw, tr.packet.W1 - wbias2, t) * inv(Tw2)';
mm2 = interp1(t_raw, tr.packet.M1 - mbias2, t) * Am2 * T2; 
am1 = interp1(t_raw, tr.packet.A2 - abias1, t) * inv(Ta1)';
wm1 = interp1(t_raw, tr.packet.W2 - wbias1, t) * inv(Tw1)';
mm1 = interp1(t_raw, tr.packet.M2 - mbias1, t) * Am1 * T1;
voltage = interp1(t_raw, tr.packet.Sg, t);

% finds midstance 
is_stance = sum(abs(voltage - nanmedian(voltage)), 2) > 0.02;

p1 = fillmissing(tr.ftrans + quatrotate(quatinv(tr.fquat), rf1'), 'linear', 'EndValues', 'nearest');
p2 = fillmissing(tr.strans + quatrotate(quatinv(tr.squat), rs2'), 'linear', 'EndValues', 'nearest');

% load ground truth measurements of foot and shank quaternion
load('/Users/ahmedsoliman/Documents/GitHub/os-eskf/Data/Prosthesis_data/fquat.mat')
load('/Users/ahmedsoliman/Documents/GitHub/os-eskf/Data/Prosthesis_data/squat.mat')

v1 = deriv_sgolay(p1, 183, [3, 5]);
v2 = deriv_sgolay(p2, 183, [3, 5]);

awmbias = zeros(length(tr.mtime), 3*3);
block_yaw = repmat(eskf_cdprosthesis.BLOCK_YAW, [length(tr.mtime), 1]);
north = repmat(mag_north, [length(tr.mtime), 1]);

NominalState = [interp1(tr.mtime, [p1, v1, fquat, awmbias], t), ...
                interp1(tr.mtime, [p2, v2, squat, awmbias, ...
                block_yaw, north], t)];
            
aangle = interp1(tr.mtime, log_quat(quatmultiply(quatinv(squat), fquat)), t);
            meas = [am2,wm2,mm2,am1,wm1,mm1,voltage];
            
%% ESKF - section 2

            
dt = 1 / Fs;
dlen = length(t);

ft = eskf_cdprosthesis();

x0 = project_state_to_constrain(ft, NominalState(1,1:3)', am1(1,:)', mm1(1,:)', am2(1,:)', mm2(1,:)'); % rt initialize

ft.set_nominal_state(x0); %

Measurement = zeros(dlen, ft.ny);
MeasurementCorrected = zeros(dlen, ft.ny);
State = zeros(dlen, ft.nx + 2);
StateVar = zeros(dlen, ft.nx);

for i = 1 : dlen
    State(i,:) = ft.nominalState; %  rt 
    StateVar(i,:) = diag(ft.errorCov_); 
    um = [am1(i,:), wm1(i,:), am2(i,:), wm2(i,:)]'; %
    ym = [mm1(i,:), mm2(i,:)]; %
    
    ft.correct_measurement(mm1(i,:)', mm2(i,:)', is_stance(i), um); %
    
    
    Measurement(i,1:length(ym)) = ym;
    MeasurementCorrected(i,:) = ft.estimate_measurement(um);

    ft.predict(um); %
end

StateVar = [StateVar(:,1:6), zeros(dlen,1), StateVar(:,7:18), ...
            StateVar(:,19:24), zeros(dlen,1), StateVar(:,25:end)];

aangle = log_quat(quatmultiply(quatinv(NominalState(:,26:29)), NominalState(:,7:10)));
aangleh = log_quat(quatmultiply(quatinv(State(:,26:29)), State(:,7:10)));


figure, set(gcf, 'position', [432.2000  151.4000  463.2000  524.0000])
h1 = subplot(211); ax = plot(t, aangleh*180/pi);
hold on, set(gca, 'ColorOrderIndex', 1)
plot(t, aangle*180/pi, '--'), grid on, title('ankle angle [deg]')
h2 = subplot(212); plot(t, (aangleh - aangle)*180/pi), grid on, title('error [deg]')
linkaxes([h1, h2], 'x'), sgtitle('ESKF'), 
legend(ax, 'IE', 'EI', 'DP')




            