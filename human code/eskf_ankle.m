classdef eskf_ankle < handle
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    % Sola, Joan. "Quaternion kinematics for the error-state Kalman filter." arXiv preprint arXiv:1711.02508 (2017).
    % https://arxiv.org/abs/1711.02508
    
    properties
        
        nominalState
        errorCov_
        
        processNoise_
        measurementCov_
    end
    
    properties (Constant)
        % update your magnetic north when changing cities:
        %   https://www.ngdc.noaa.gov/geomag/calculators/magcalc.shtml#igrfwmm
%         MAG_NORTH = [1, 0, 0; 0, 0, 1; 0, -1, 0] * [-1.5, 19.9, -48.6]' % QTM, West Lafayette
        MAG_NORTH = Ry(10*pi/180) * [0;-39.2642881932957;-20.7460923566855]
        GRAV = [0; -9.80124; 0]
        
        
%                 Human subject parameters
        rfa = [-0.00400468951415383;-0.00696457661573292;0.00981857670826987];
        rsa = [0.000353893495041279;-0.212277798170128;0.00120069587260607];
        rf1 =[0.0913052623003641;-0.00601773804459769;0.00576289430437782];
        rs2 = [-0.0120421461130547;-0.0209334536188618;-0.0406870370938923];

        % IMU intrinsic calib (5h calib), 1 slope
        Sa = ones(1,3) * 2.4241e-02
        Sw = ones(1,3) * 7.3556e-03
        Sm = ones(1,3) * 1.5205e+00 * 5  % multiply to compensate mag disturbance
        Sab = ones(1,3) * 5.7844e-04
        Swb = ones(1,3) * 7.5423e-05
        Smb = ones(1,3) * 1.1920e-01
        
        DT = 1 / 400
        R1A = -eskf_ankle.rf1 + eskf_ankle.rfa
        R2A = -eskf_ankle.rs2 + eskf_ankle.rsa
        
        RCOP = [0;0;0]  % ankle to cop. Disable when walking is not realistic
        
        nx = 18 * 2 + 3
        nw = 15 * 2
        
        ny = 3 * 2 + 3 + 3 + 3
    end
    
    methods
        function obj = eskf_ankle()
            
            obj.processNoise_ = diag([obj.Sa, obj.Sw, obj.Sab, obj.Swb, obj.Smb, ...
                                      obj.Sa, obj.Sw, obj.Sab, obj.Swb, obj.Smb].^2);

            obj.measurementCov_ = diag([obj.Sm * 10, obj.Sm, ...
                ones(1,3)*1e-2, ones(1,3)*2e-2, ones(1,3)*10e-1].^2);

            StdQuat0 = [1, 1, 1];  % 57 deg std error
            StdTrans0 = [0.5, 0.5, 0.5];  % shank trans error irt foot
            StdV0 = [2,2,2];  % use peak vel during walk
            StdWb = obj.Swb * 3;  % should be proportional to bias walk noise
            StdAb = obj.Sab * 3;
            StdMb = obj.Smb * 3;
            StdMagNorth = [0, 10, 10] * 0;
            obj.errorCov_ = diag([StdTrans0*0, StdV0, StdQuat0, StdAb, StdWb, StdMb, ...
                                  StdTrans0,   StdV0, StdQuat0, StdAb, StdWb, StdMb, ...
                                  StdMagNorth].^2);  % foot trans is unobservable, so make it null
        end
        
        function xproj = project_state_to_constrain(obj, s1, am1, mm1, am2, mm2)

            yaxis = am1 / sqrt(am1' * am1);
            mag_north_body = -mm1 / sqrt(mm1' * mm1);  % towards (-Y-Z), X~=0
            xaxis = cross(yaxis, mag_north_body);
            xaxis = xaxis / sqrt(xaxis' * xaxis);
            zaxis = cross(xaxis, yaxis);
            Rf = [xaxis, yaxis, zaxis]';
            
            yaxis = am2 / sqrt(am2' * am2);
            mag_north_body = -mm2 / sqrt(mm2' * mm2);  % towards (-Y-Z), X~=0
            xaxis = cross(yaxis, mag_north_body);
            xaxis = xaxis / sqrt(xaxis' * xaxis);
            zaxis = cross(xaxis, yaxis);
            Rs = [xaxis, yaxis, zaxis]';
            
            fquat = rotm2quat(Rf);
            squat = rotm2quat(Rs);
            s2 = s1 + Rf * obj.R1A - Rs * obj.R2A;
                               
            xproj = [s1; zeros(3,1); fquat'; zeros(9,1)
                     s2; zeros(3,1); squat'; zeros(9,1)
                     obj.MAG_NORTH];
        end
        
        function set_nominal_state(obj, varargin)
            assert(nargin == 2, 'wrong input args')
            
            x = varargin{1};
            obj.nominalState = x;
        end
        
        function x = get_nominal_state(obj)
            x = obj.nominalState;
        end
        
        
        function [yh, H] = estimate_measurement(obj, U)
            x_nominal = get_nominal_state(obj);
            % X, U, grav, mag_north, dt, r1a, r2a, rcop
            [yh, H] = eskf_ankle_measurement(x_nominal, U, obj.GRAV, obj.DT, obj.R1A, obj.R2A, obj.RCOP);
        end
        
        function obj = predict(obj, U)
            
            x = get_nominal_state(obj);
            processNoise = obj.processNoise_;
            
            % error state ================================================
            [xk, A, G] = eskf_ankle_prediction(x, U, obj.GRAV, obj.DT);  % before nominal?
            set_nominal_state(obj, xk);

            errorCov = get_error_state(obj);
            errorCov = A * errorCov * A' + G * processNoise * G';
            set_error_state(obj, errorCov);
        end


        function obj = correct_measurement(obj, mm1, mm2, is_stance, U)
            [yh, H] = estimate_measurement(obj, U);
            ym = [mm1; mm2; zeros(3,1); zeros(3,1); zeros(3,1)];
            % y: [m1; m2; atrans; vfstance, wfstance]
            
            rows = [1:6, 7:9, 10:12, 13:15];
            
            % ignore mag if disturbance is too large
            mag_tol = 20;
            mag_norm = sqrt(obj.MAG_NORTH' * obj.MAG_NORTH);
            
            mm1_abs = sqrt(mm1' * mm1);
            if abs(mm1_abs - mag_norm) > mag_tol
                rows = setdiff(rows, 1:3);
            end
            
            mm2_abs = sqrt(mm2' * mm2);
            if abs(mm2_abs - mag_norm) > mag_tol
                rows = setdiff(rows, 4:6);
            end
            
            if ~is_stance
                rows = setdiff(rows, 10:15);
            end
            
            measurementCov = obj.measurementCov_(rows,rows);
            H = H(rows,:);
            res = ym(rows) - yh(rows);
            
            % KF correct equations ==================
            errorCov = get_error_state(obj);
            
            I = eye(obj.nx);
            K = errorCov * H' / (H * errorCov * H' + measurementCov);
            dx = K * res;
            errorCov = (I - K*H) * errorCov * (I - K*H)' + K * measurementCov * K';
            
            % inject error into nominal states ===========================
            x_nominal = get_nominal_state(obj);
            [x_post, G] = eskf_ankle_reset(x_nominal, dx); % X, dX 
            set_nominal_state(obj, x_post);
            
            % reset error ================================================            
            errorCov = G * errorCov * G.';
            set_error_state(obj, errorCov);
        end

        
        function set_error_state(obj, errorCov)
            obj.errorCov_ = errorCov;
        end
        
        function [errorCov] = get_error_state(obj)
           errorCov = obj.errorCov_;
        end
        
    end
    
    methods (Static)
        
        function generate_functions()
        
            % constant params
            mag_north = sym('mag_north', [3, 1]);
            dmag_north = sym('dmag_north', [3, 1]);
            mag_north_true = mag_north + dmag_north;
            
            grav = sym('grav', [3, 1]);
            syms dt
            r1a = sym('r1a', [3, 1]);
            r2a = sym('r2a', [3, 1]);
            rcop = sym('rcop', [3, 1]);
            
            [X_f, dX_f, Xt_f, dXk_f, Xk_f, y_f, U_f, Wn_f, reset_fcn_f] = ...
                eskf_ankle.calculate_imu_equations('foot', mag_north_true, grav, dt);
            [X_s, dX_s, Xt_s, dXk_s, Xk_s, y_s, U_s, Wn_s, reset_fcn_s] = ...
                eskf_ankle.calculate_imu_equations('shank', mag_north_true, grav, dt);
            
            % additional states
            
            % geometrical constraints
            p1 = Xt_f(1:3);
            v1 = Xt_f(4:6);
            Rf = quat2rotm_sym(Xt_f(7:10).');
            p2 = Xt_s(1:3);
            Rs = quat2rotm_sym(Xt_s(7:10).');
            
            pah = p1 + Rf * r1a - (p2 + Rs * r2a);  % translation constraint
            wf = U_f(4:6) - Xt_f(14:16);
            vstance = v1 + Rf * cross(wf, rcop);  % null foot speed at stance
            
            % combine all ==================================
            X = [X_f; X_s; mag_north];  % nominal state
            dX = [dX_f; dX_s; dmag_north];  % error states
            Xt = [Xt_f; Xt_s; mag_north_true];  % true states X (+) dX
            dXk = [dXk_f; dXk_s; dmag_north];
            Xk = [Xk_f; Xk_s; mag_north];
            Y = [y_f; y_s
                 pah
                 vstance
                 wf];  % nominal measurement
            Wn = [Wn_f; Wn_s];   
            U = [U_f; U_s];
            
            nx = size(dX, 1);
            nw = size(Wn, 1);
            
            reset_fcn = @(dx1, dx2) [reset_fcn_f(dx1(1:18), dx2(1:18))
                                     reset_fcn_s(dx1(19:36), dx2(19:36))
                                     dx1(37:39) - dx2(37:39)];
            dX_correct = sym('dX_correct', [nx, 1]);
            reset_dx = subs(jacobian(reset_fcn(dX, dX_correct), dX), dX_correct, dX);

            yout = subs(Y, [dX; Wn], zeros(nx+nw, 1));
            
            % full vector valid only during stance (remove vcop during swing)
            y_dx = subs(jacobian(Y, dX), [dX; Wn], zeros(nx+nw, 1));
            xk_dx = subs(jacobian(dXk, dX), [dX; Wn], zeros(nx+nw, 1));
            xk_w = subs(jacobian(dXk, Wn), [dX; Wn], zeros(nx+nw, 1));
            
            matlabFunction(yout, y_dx, 'Vars', {X, U, grav, dt, r1a, r2a, rcop}, 'file', 'eskf_ankle_measurement')
            matlabFunction(Xk, xk_dx, xk_w, 'Vars', {X, U, grav, dt}, 'file', 'eskf_ankle_prediction')
            matlabFunction(Xt, reset_dx, 'Vars', {X, dX}, 'file', 'eskf_ankle_reset')
        end
      
        function [X, dX, Xt, dXk, Xk, Y, U, Wn, reset_fcn] = ...
                calculate_imu_equations(prefix, mag_north, grav, dt)

            % nominal
            quat = sym([prefix, 'quat'], [1, 4]);
            p = sym([prefix, 'p'], [3, 1]);
            v = sym([prefix, 'v'], [3, 1]);
            wbias = sym([prefix, 'wbias'], [3, 1]);
            abias = sym([prefix, 'abias'], [3, 1]);
            mdist = sym([prefix, 'mdist'], [3, 1]);
            X = [p; v; quat.'; abias; wbias; mdist];

            % error
            dp = sym([prefix, 'dp'], [3, 1]);
            dv = sym([prefix, 'dv'], [3, 1]);
            dq = sym([prefix, 'dq'], [3, 1]);
            dwbias = sym([prefix, 'dwbias'], [3, 1]);
            dabias = sym([prefix, 'dabias'], [3, 1]);
            dmdist = sym([prefix, 'dmdist'], [3, 1]);
            dX = [dp; dv; dq; dabias; dwbias; dmdist];

            % IMU measurements
            wm = sym([prefix, 'wm'], [3, 1]);
            am = sym([prefix, 'am'], [3, 1]);

            % noise
            wn = sym([prefix, 'wn'], [3, 1]);
            an = sym([prefix, 'an'], [3, 1]);
            ww = sym([prefix, 'ww'], [3, 1]);
            aw = sym([prefix, 'aw'], [3, 1]);
            mw = sym([prefix, 'mw'], [3, 1]);

            % update nominal equations, Euler integration
            R = quat2rotm_sym(quat);
            w = wm - wbias;
            a = R * (am - abias) + grav;
            quatk = quatmultiply_sym(quat, [1, 0.5 * w.' * dt]);
            quatk = quatk / sqrt(quatk * quatk.');
            pk = p + v * dt + 0.5 * a * dt^2;
            vk = v + a * dt;
            wbiask = wbias;
            abiask = abias;
            mdistk = mdist;
            
            % true, global
            quatt = quatmultiply_sym([1, 0.5*dq.'], quat);
            quatt = quatt / sqrt(quatt * quatt.');
            pt = p + dp;
            vt = v + dv;
            wbiast = wbias + dwbias;
            abiast = abias + dabias;
            mdistt = mdist + dmdist;
            Xt = [pt; vt; quatt.'; abiast; wbiast; mdistt];
            
            % update true equations, Euler integration
            Rt = quat2rotm_sym(quatt);
            wt = wm - wbiast - wn;
            at = Rt * (am - abiast - an) + grav;
            quattk = quatmultiply_sym(quatt, [1, 0.5 * wt.' * dt]);
            quattk = quattk / sqrt(quattk * quattk.');
            ptk = pt + vt * dt + 0.5 * at * dt^2;
            vtk = vt + at * dt;
            wbiastk = wbiast + ww * dt;
            abiastk = abiast + aw * dt;
            mdisttk = mdistt + mw * dt;
            
            % update nominal equations for prediction step, RK4
            fcn = @(x) eskf_ankle.imu_continuous_sym(x, am, wm, grav);
            Xk = eskf_ankle.RK4(fcn, X, dt);
            
            % update equations, error
            dpk = ptk - pk;
            dvk = vtk - vk;
            dabiask = abiastk - abiask;
            dwbiask = wbiastk - wbiask;
            dmdistk = mdisttk - mdistk;
            dqk_aug = quatmultiply_sym(quattk, quatconj_sym(quatk));
            dqk = 2 * dqk_aug(2:4).';
            
            % measurements
            Rt = quat2rotm_sym(quatt);
            mh = Rt.' * (mag_north) + mdistt;

            % define system vectors
            Wn = [an; wn; aw; ww; mw];  % process noise
            
            dXk = [dpk; dvk; dqk; dabiask; dwbiask; dmdistk];
            Y = mh;
            U = [am; wm];
            
            reset_fcn = @(dx1, dx2)  [dx1(1:6) - dx2(1:6)
                     2 * [zeros(3,1), eye(3)] * quatnormalize_sym(quatmultiply_sym(quatconj_sym([1, 0.5*dx2(7:9).']), [1, 0.5*dx1(7:9).'])).'
                     dx1(10:18) - dx2(10:18)];
        end

        function xk = RK4(fcn, x, dt)
            k1 = fcn(x);
            k2 = fcn(x + k1*dt/2);
            k3 = fcn(x + k2*dt/2);
            k4 = fcn(x + k3*dt);
            xk = x + dt*(k1 + 2*k2 + 2*k3 + k4)/6;
        end
        
        function xd = imu_continuous(x, am, wm)
           % x: [p, v, quat, abias, wbias, mdist]
           pd = x(4:6);
           quat = x(7:10)';
           pdd = quatrotate(quatinv(quat), am' - x(11:13)')' + eskf_ankle.GRAV;
           quatd = 0.5 * quatmultiply(quat, [0; wm - x(14:16)]');
           xd = [pd; pdd; quatd'; zeros(9,1)];
        end
        
        function xd = imu_continuous_sym(x, am, wm, grav, decay)
           % x: [p, v, quat, abias, wbias, mdist]
           pd = x(4:6);
           quat = x(7:10).';
           R = quat2rotm_sym(quat);
           pdd = R * (am - x(11:13)) + grav;
           quatd = 0.5 * quatmultiply_sym(quat, [0; wm - x(14:16)].');
           xd = [pd; pdd; quatd.'; zeros(9,1)];
        end
    end
end

