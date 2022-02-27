classdef eskf_cdprosthesis < handle
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
        MAG_NORTH = Ry(0*pi/180) * [1, 0, 0; 0, 0, 1; 0, -1, 0] * [-1.5, 19.9, -48.6]' % QTM, West Lafayette
        GRAV = [0; -9.82; 0]
        
        % mocap extrinsic calib
        BLOCK_YAW = 0.0055 % heading of the CNC block
        rfa = [0.000481025186908457;0.00104955038446195;-0.00186800748296271]
        rsa = [-0.00213087026849489;-0.101001412610616;0.000779715520280068]
        rs2 = [-0.042883813790953;-0.023992214323566;0.002762408594112]
        
        Ta1 = [0.011323950340988,-0.016939494746075,1.012707236066504;0.935161296729301,0.005267291078370,-0.007070506492664;-0.008920418375473,1.026338177998765,0.018630511643839]
        Tw1 = [-0.001275146115968,-0.022862075633627,0.990674651695081;1.010846175757672,0.017889462773747,-0.006408442508446;-0.018031823450745,0.994265373776011,0.022056667801162]
        wbias1 = [0.009467439417936,-0.015054962876172,0.004634200094612]
        abias1 = [0.384721052261843,-0.457519405371102,0.415731828505475]
        
        % mag1 intrinsic calib
        mbias1 = [23.7386562435748,-17.1345013587654,-30.8978243021282]
        Am1 = eye(3)
        expmfs1 = 45.7393193183346
        
        % IMU2 extrinsic calib //shank data
        rf1 = [0.004322051878450;0.003384083223609;0.024829444850935]
        Ta2 = [0.982044909548764,-0.050156440745373,-1.898640891108553e-04;0.041821128521680,0.952917087317640,0.031887607789532;0.012467547622602,-0.032481642029070,1.023180592777490]
        Tw2 =[1.009146514797720,-0.052065512761995,-0.017693980892473;0.046916021407917,0.998718463136889,0.033366571338913;0.013622785132938,-0.035834488336799,0.998664250224121]
        wbias2 = [0.014944959276279,-0.001512418745301,0.016279034895197]
        abias2 = [0.329802179379453,-0.536283620711767,0.251053428471307]
        
        % mag1 intrinsic calib
        mbias2 = [-95.9651549017379,-103.511791651512,-40.9586003285336]
        Am2 = [1.04992208372345,0.0278139786763245,0.0356850265397141;0.0278139786763245,0.942095177623534,0.0854073703239322;0.0356850265397141,0.0854073703239322,1.02057544389801]
        expmfs2 = 43.6228328944328 % mm_norm = (mm - b) * A        
        
        % IMU intrinsic calib (5h calib), 1/2 slope
        Sa = ones(1,3) * 3.7926e-02
        Sw = ones(1,3) * 1.1221e-02
        Sm = ones(1,3) * 3.8385e+00 * 5  % multiply to compensate mag disturbance
        Sab = ones(1,3) * 2.3639e-03
        Swb = ones(1,3) * 4.0155e-04
        Smb = ones(1,3) * 5.6858e-01
        
        DT = 1 / 400
        R1A = -eskf_cdprosthesis.rf1 + eskf_cdprosthesis.rfa
        R2A = -eskf_cdprosthesis.rs2 + eskf_cdprosthesis.rsa
        
        nx = 18 * 2 + 1 + 3
        nw = 15 * 2
        
        ny = 3 * 2 + 3 + 1 + 3
        nv = 3 * 2  % additive noise
    end
    
    methods
        function obj = eskf_cdprosthesis()
            
            obj.processNoise_ = diag([obj.Sa, obj.Sw, obj.Sab, obj.Swb, obj.Smb / obj.expmfs1, ...
                                      obj.Sa, obj.Sw, obj.Sab, obj.Swb, obj.Smb / obj.expmfs2].^2);

            obj.measurementCov_ = diag([obj.Sm / obj.expmfs1, obj.Sm / obj.expmfs2, ...
                ones(1,3)*1e-3, 1e-3, ones(1,3)*1e-2].^2);

            StdQuat0 = [1, 1, 1];  % 57 deg std error
            StdTrans0 = [0.5, 0.5, 0.5];  % shank trans error irt foot
            StdV0 = [2,2,2];  % use peak vel during walk
            StdWb = obj.Swb * 3;  % should be proportional to bias walk noise
            StdAb = obj.Sab * 3;
            StdMb = 3 * obj.Smb / obj.expmfs1;
            StdYaw = 0.1 * pi/180 * 0;  % make null if known
            StdNorth = [1,1,1]*0.2 * 0;  % make null if known
            obj.errorCov_ = diag([StdTrans0*0, StdV0, StdQuat0, StdAb, StdWb, StdMb, ...
                                  StdTrans0, StdV0, StdQuat0, StdAb, StdWb, StdMb, ...
                                  StdYaw, StdNorth].^2);  % foot trans is unobservable, so make it null
        end
        
        function xproj = project_state_to_constrain(obj, s1, am1, mm1, am2, mm2)

            yaxis = am1 / sqrt(am1' * am1);
            mag_north = -mm1 / sqrt(mm1' * mm1);  % towards (-Y-Z), X~=0
            xaxis = cross(yaxis, mag_north);
            xaxis = xaxis / sqrt(xaxis' * xaxis);
            zaxis = cross(xaxis, yaxis);
            R = [xaxis, yaxis, zaxis]';
            
            yaxis = am2 / sqrt(am2' * am2);
            mag_north = -mm2 / sqrt(mm2' * mm2);  % towards (-Y-Z), X~=0
            xaxis = cross(yaxis, mag_north);
            xaxis = xaxis / sqrt(xaxis' * xaxis);
            zaxis = cross(xaxis, yaxis);
            R2 = [xaxis, yaxis, zaxis]';
            
%             fquat = rotm2quat(R + R2);
%             squat = fquat;
            
            fquat = rotm2quat(R);
            squat = rotm2quat(R2);
            [ry, rz, rx] = quat2angle(quatmultiply(quatinv(squat), fquat), 'YZX');
            aquat = angle2quat(obj.BLOCK_YAW, rz, rx, 'YZX');
            squat = quatmultiply(quatinv(aquat), fquat);
            
            s2 = s1 + quatrotate(quatinv(fquat), -obj.rf1' + obj.rfa')' + ...
                       quatrotate(quatinv(squat), -obj.rsa' + obj.rs2')';
                   
            xproj = [s1; zeros(3,1); fquat'; zeros(9,1)
                     s2; zeros(3,1); squat'; zeros(9,1)
                     obj.BLOCK_YAW; obj.MAG_NORTH / sqrt(sum(obj.MAG_NORTH.^2))];
        end
        
        function set_nominal_state(obj, varargin)
            if nargin == 2
                x = varargin{1};
            elseif nargin == 13
                x = [varargin{1}; varargin{2}; varargin{3}'; varargin{4}; varargin{5}; varargin{6}
                     varargin{7}; varargin{8}; varargin{9}'; varargin{10}; varargin{11}; varargin{12}];
            else
                error("Set nominal state as array or [p, v, quat, abias, wbias, mdist]")
            end
            obj.nominalState = x;
        end
        
        function x = get_nominal_state(obj)
            x = obj.nominalState;
        end
        
        function [dt, grav, processNoise, mag_north, measurementCov, r1a, r2a] = get_params(obj)
            dt = obj.DT;
            grav = obj.GRAV;
            processNoise = obj.processNoise_;
            mag_north = obj.MAG_NORTH / sqrt(obj.MAG_NORTH' * obj.MAG_NORTH);
            measurementCov = obj.measurementCov_;
            r1a = obj.R1A;
            r2a = obj.R2A;
        end
        
        
        function [yh, H] = estimate_measurement(obj, U)
            [dt, grav, ~, mag_north, ~, r1a, r2a] = get_params(obj);
            x_nominal = get_nominal_state(obj);
            [yh, H] = eskf_cdprosthesis_measurement(x_nominal, U, grav, dt, r1a, r2a);
        end
        
        function obj = predict(obj, U)
            
            x = get_nominal_state(obj);
            [dt, grav, processNoise, mag_north] = get_params(obj);
            
            % error state ================================================
            [xk, A, G] = eskf_cdprosthesis_prediction(x, U, grav, dt);  % before nominal?
            set_nominal_state(obj, xk);

            errorCov = get_error_state(obj);
            errorCov = A * errorCov * A' + G * processNoise * G';
            set_error_state(obj, errorCov);
        end


        function obj = correct_measurement(obj, mm1, mm2, is_stance, U)

            [yh, H] = estimate_measurement(obj, U);
            ym = [mm1; mm2; zeros(3,1); 0.0; zeros(3,1)];
            
            rows = [1:6, 7:9, 10, 11:13];  % atrans and vcop
            
            
            mag_tol = 0.2;  
            mm1_abs = sqrt(mm1' * mm1);
            
%             mm1_abs = sqrt(complex(mm1' * mm1));
            if abs(mm1_abs - 1) > mag_tol
                rows = setdiff(rows, 1:3);
            end
            
            mm2_abs = sqrt(mm2' * mm2);
            if abs(mm2_abs - 1) > mag_tol
                rows = setdiff(rows, 4:6);
            end
            
            if ~is_stance
                rows = setdiff(rows, 11:13);
            end
            
            measurementCov = obj.measurementCov_(rows,rows);
            H = H(rows,:);
            res = ym(rows) - yh(rows);

            correct(obj, res, measurementCov, H);
        end
        
        function obj = correct(obj, res, measurementCov, H)
            
            errorCov = get_error_state(obj);
            
            I = eye(obj.nx);
            K = errorCov * H' / (H * errorCov * H' + measurementCov);
            dx = K * res;
            errorCov = (I - K*H) * errorCov * (I - K*H)' + K * measurementCov * K';
            
            % inject error into nominal states ===========================
            x_nominal = get_nominal_state(obj);
            [x_post, G] = eskf_cdprosthesis_reset(x_nominal, dx); % X, dX 
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
%             qEI = sym('qEI', 1);
            
            [X_f, dX_f, Xt_f, dXk_f, Xk_f, y_f, U_f, Wn_f, reset_fcn_f] = ...
                eskf_cdprosthesis.calculate_imu_equations('foot', mag_north_true, grav, dt);
            [X_s, dX_s, Xt_s, dXk_s, Xk_s, y_s, U_s, Wn_s, reset_fcn_s] = ...
                eskf_cdprosthesis.calculate_imu_equations('shank', mag_north_true, grav, dt);
            
            % additional states
            block_yaw = sym('block_yaw', 1);  % block yaw error on pylon
            dblock_yaw = sym('dblock_yaw', 1);
            block_yaw_true = block_yaw + dblock_yaw;
            
            % geometrical constraints
            p1 = Xt_f(1:3);
            v1 = Xt_f(4:6);
            Rf = quat2rotm_sym(Xt_f(7:10).');
            p2 = Xt_s(1:3);
            Rs = quat2rotm_sym(Xt_s(7:10).');
            Ra = Rs.' * Rf;
            
            pah = p1 + Rf * r1a - (p2 + Rs * r2a);  % translation constraint
            qEI = atan(-Ra(3,1) / Ra(1,1)) - block_yaw_true;  % EI rotation constraint
            vstance = v1;  % null foot speed at stance
            
            % combine all ==================================
            X = [X_f; X_s; block_yaw; mag_north];  % nominal state
            dX = [dX_f; dX_s; dblock_yaw; dmag_north];  % error states
            Xt = [Xt_f; Xt_s; block_yaw_true; mag_north_true];  % true states X (+) dX
            dXk = [dXk_f; dXk_s; dblock_yaw; dmag_north];
            Xk = [Xk_f; Xk_s; block_yaw; mag_north];
            Y = [y_f; y_s
                 pah
                 qEI
                 vstance];  % nominal measurement
            Wn = [Wn_f; Wn_s];   
            U = [U_f; U_s];
            
            nx = size(dX, 1);
            nw = size(Wn, 1);
            
            reset_fcn = @(dx1, dx2) [reset_fcn_f(dx1(1:18), dx2(1:18))
                                     reset_fcn_s(dx1(19:36), dx2(19:36))
                                     dx1(37:40) - dx2(37:40)];
            dX_correct = sym('dX_correct', [nx, 1]);
%             reset_dx = subs(jacobian(reset_fcn(dX, dX_correct), dX_correct), dX_correct, dX);
            reset_dx = subs(jacobian(reset_fcn(dX, dX_correct), dX), dX_correct, dX);

            yout = subs(Y, [dX; Wn], zeros(nx+nw, 1));
            
            % full vector valid only during stance (remove vcop during swing)
            y_dx = subs(jacobian(Y, dX), [dX; Wn], zeros(nx+nw, 1));
            xk_dx = subs(jacobian(dXk, dX), [dX; Wn], zeros(nx+nw, 1));
            xk_w = subs(jacobian(dXk, Wn), [dX; Wn], zeros(nx+nw, 1));
            
            matlabFunction(yout, y_dx, 'Vars', {X, U, grav, dt, r1a, r2a}, 'file', 'eskf_cdprosthesis_measurement')
            matlabFunction(Xk, xk_dx, xk_w, 'Vars', {X, U, grav, dt}, 'file', 'eskf_cdprosthesis_prediction')
            matlabFunction(Xt, reset_dx, 'Vars', {X, dX}, 'file', 'eskf_cdprosthesis_reset')
        end
        
        function generate_tuning_functions()
            
            % unknown states
            %    cqf, csf, cqs, css
     
            % measurements
            %    wm1, am1, mm1, wm2, am2, mm2, mftrans, mstrans
            
            % constraints
            %    atrans_constr, qEI_constr, vstance_constr
            
            Ta1 = sym('Ta1', 3);
            Tw1 = sym('Tw1', 3);
            Tm1 = sym('Tm1', 3);
            Ta2 = sym('Ta2', 3);
            Tw2 = sym('Tw2', 3);
            Tm2 = sym('Tm2', 3);
            rf1 = sym('rf1', [3, 1]);
            rs2 = sym('rs2', [3, 1]);
            abias1 = sym('abias1', [3, 1]);
            wbias1 = sym('wbias1', [3, 1]);
            mbias1 = sym('mbias1', [3, 1]);
            abias2 = sym('abias2', [3, 1]);
            wbias2 = sym('wbias2', [3, 1]);
            mbias2 = sym('mbias2', [3, 1]);
            grav = sym('grav', [3, 1]);
            mag_north = sym('mag_north', [3, 1]);
            r1a = sym('r1a', [3, 1]);
            r2a = sym('r2a', [3, 1]);
            block_yaw = sym('block_yaw', 1);
            tshift = sym('tshift', 1);
            
            % fixed params
            rfi = sym('rfi', [3, 1]);
            rsi = sym('rsi', [3, 1]);
            
            % compute derivated states
            %    Rf, wf, p1, v1, a1, Rs, ws, p2, v2, a2
            Rf = sym('Rf', 3);
            wf = sym('wf', [3, 1]);
            p1 = sym('p1', [3, 1]);
            v1 = sym('v1', [3, 1]);
            a1 = sym('a1', [3, 1]);
            Rs = sym('Rs', 3);
            ws = sym('ws', [3, 1]);
            p2 = sym('p2', [3, 1]);
            v2 = sym('v2', [3, 1]);
            a2 = sym('a2', [3, 1]);
            
            Ra = Rs.' * Rf;
            
            wm1 = Tw1 * wf + wbias1;
            wm2 = Tw2 * ws + wbias2;
            mm1 = Tm1 * Rs.' * mag_north + mbias1;
            mm2 = Tm2 * Rf.' * mag_north + mbias2;
            am1 = Ta1 * Rf.' * (a1 - grav) + abias1;
            am2 = Ta2 * Rs.' * (a2 - grav) + abias2;
            
            mftrans = p1 + Rf * (-rf1 + rfi);
            mstrans = p2 + Rs * (-rs2 + rsi);
            
            atrans_constr = p1 + Rf * (r1a) - (p2 + Rs * (r2a));
            qEI_constr = atan(-Ra(3,1) / Ra(1,1)) - block_yaw;
            % vstance_constr = v1;  % redundant given ftrans as measurement
            
            % unknown parameters
            params = [Ta1(:); Tw1(:); Tm1(:); Ta2(:); Tw2(:); Tm2(:); rf1; rs2
                      abias1; wbias1; mbias1; abias2; wbias2; mbias2
                      grav; mag_north; r1a; r2a; block_yaw];
                  
            Yimu = [am1; wm1; mm1; am2; wm2; mm2];
            Yomc = [mftrans; mstrans];
            Yconstr = [atrans_constr; qEI_constr];  % TODO add biasd
            
%             Jimu = jacobian(Yimu, params);
            matlabFunction(Yimu, 'File', 'cdprosthesis_tune_imu', 'Vars', {params, Rf, wf, p1, v1, a1, Rs, ws, p2, v2, a2, rfi, rsi})
            matlabFunction(Yomc, 'File', 'cdprosthesis_tune_omc', 'Vars', {params, Rf, wf, p1, v1, a1, Rs, ws, p2, v2, a2, rfi, rsi})
            matlabFunction(Yconstr, 'File', 'cdprosthesis_tune_constr', 'Vars', {params, Rf, wf, p1, v1, a1, Rs, ws, p2, v2, a2, rfi, rsi})
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
            fcn = @(x) eskf_cdprosthesis.imu_continuous_sym(x, am, wm, grav);
            Xk = eskf_cdprosthesis.RK4(fcn, X, dt);
            
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
        
        function xd = cdprosthesis_continuous(x, U)
            xd = [eskf_cdprosthesis.imu_continuous(x(1:19), U(1:3), U(4:6))
                  eskf_cdprosthesis.imu_continuous(x(20:38), U(7:9), U(10:12))];
        end
        
        function xd = imu_continuous(x, am, wm)
           % x: [p, v, quat, abias, wbias, mdist]
           pd = x(4:6);
           quat = x(7:10)';
           pdd = quatrotate(quatinv(quat), am' - x(11:13)')' + eskf_cdprosthesis.GRAV;
           quatd = 0.5 * quatmultiply(quat, [0; wm - x(14:16)]');
           xd = [pd; pdd; quatd'; zeros(9,1)];
        end
        
        function xd = imu_continuous_sym(x, am, wm, grav, decay)
           pd = x(4:6);
           quat = x(7:10).';
           R = quat2rotm_sym(quat);
           pdd = R * (am - x(11:13)) + grav;
           quatd = 0.5 * quatmultiply_sym(quat, [0; wm - x(14:16)].');
           xd = [pd; pdd; quatd.'; zeros(9,1)];
        end
    end
end

