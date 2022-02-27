clear
clc
close all

%% Load Human_data 
% ----------------------------------------------------------------------------------------------------------------------------------------------------------
% -- Download from : https://purdue0-my.sharepoint.com/:f:/g/personal/solimana_purdue_edu/EqfMwu4GRrFFm5kUIA9fDiIBGnNZ-kNACIa6cl1AKxoxsg?e=brtQnd

f = load('/Users/ahmedsoliman/Documents/GitHub/os-eskf/Data/Human_data/calib_imu1_foot_a_w.mat');
s =load('/Users/ahmedsoliman/Documents/GitHub/os-eskf/Data/Human_data/calib_imu2_shank_a_w.mat');
f_m = load('/Users/ahmedsoliman/Documents/GitHub/os-eskf/Data/Human_data/calib_imu1_foot_m.mat');
s_m =load('/Users/ahmedsoliman/Documents/GitHub/os-eskf/Data/Human_data/calib_imu2_shank_m.mat');

r1 = f.r1;
Tw1 = f.Tw1;
Ta1 = f.Ta1;
Tm1 = f_m.Tm1;
abias1 = mean(f.cabias1);
wbias1 = mean(f.cwbias1);
mbias1 = (f_m.mbias1);

r2 = s.r1;
Tw2 = s.Tw1;
Ta2 = s.Ta1;
Tm2 = s_m.Tm2;
abias2 = mean(s.cabias1);
wbias2 = mean(s.cwbias1);
mbias2 = (s_m.mbias2);

rfa = f.rfa;
rsa = f.rsa;

mag_north = eskf_ankle.MAG_NORTH;
Fs = 400;

load('/Users/ahmedsoliman/Documents/GitHub/os-eskf/Data/Human_data/outdoor walk 40min.mat')


t = trial;

%% 
time = (t.itime(1) : 1/Fs : 1800)';

dlen = length(time);

wf = interp1(t.itime, (t.w1 - wbias1) * inv(Tw1)', time);
af = interp1(t.itime, (t.a1 - abias1) * inv(Ta1)', time);
mf = interp1(t.itime, (t.m1 - mbias1) * Tm1, time);

ws = interp1(t.itime, (t.w2 - wbias2) * inv(Tw2)', time);
as = interp1(t.itime, (t.a2 - abias2) * inv(Ta2)', time);
ms = interp1(t.itime, (t.m2 - mbias2) * Tm2, time);

% Label stance
wftol = .55;
aftol = 1;
pitchtol = 30 * pi/180;

yf = af ./ sqrt(sum(af.^2, 2));

is_stance = (sqrt(sum(wf.^2, 2)) < wftol) & ...
    (abs(sqrt(sum(af.^2, 2)) - 9.81) < aftol) & ...
    (abs(yf(:,2)) > cos(pitchtol));



NominalState = zeros(dlen, 41);
NominalState(:,end-2:end) = repmat(mag_north', dlen, 1);


%% ESKF

warning('This code takes up to 5 minutes to run (depending on machine)');

ft = eskf_ankle();

% x0 = NominalState(1,:)';
s10 = NominalState(1,1:3)';
x0 = ft.project_state_to_constrain(s10, af(1,:)', mf(1,:)', as(1,:)', ms(1,:)');

ft.set_nominal_state(x0);

Measurement = zeros(dlen, ft.ny);
MeasurementCorrected = zeros(dlen, ft.ny);
State = zeros(dlen, ft.nx + 2);
StateVar = zeros(dlen, ft.nx);

for i = 1 : dlen
    State(i,:) = ft.nominalState;
    StateVar(i,:) = diag(ft.errorCov_);
    um = [af(i,:), wf(i,:), as(i,:), ws(i,:)]';
    ym = [mf(i,:), ms(i,:)];
    
    ft.correct_measurement(mf(i,:)', ms(i,:)', is_stance(i), um);
    
    Measurement(i,1:length(ym)) = ym;
    MeasurementCorrected(i,:) = ft.estimate_measurement(um);

    ft.predict(um);
    if mod(i, round(dlen/10)) == 0
        fprintf('Done : %d%%, ', round(i*100/dlen))
        fprintf('\n')
    end
end
fprintf('\n')

StateVar = [StateVar(:,1:6), zeros(dlen,1), StateVar(:,7:18), ...
            StateVar(:,19:24), zeros(dlen,1), StateVar(:,25:end)];

aangle = log_quat(quatmultiply(quatinv(NominalState(:,26:29)), NominalState(:,7:10)));
aangleh = log_quat(quatmultiply(quatinv(State(:,26:29)), State(:,7:10)));

show_stance = max(max(aangleh*180/pi)) / 10;

figure, set(gcf, 'position', [432.2000  151.4000  463.2000  524.0000])
plot(time, aangleh*180/pi);
grid on, title('ankle angle [deg]')
legend('IE', 'EI', 'DP')
xlim([0 2000])

