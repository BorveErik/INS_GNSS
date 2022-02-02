%% INS/GNNS implementation for a circular motion in 2D
clear; clc; format compact; close all;

% - Seed -
seed = 10;
rng(seed);

% - Physical constants -
r = 500;
w = pi/100;
g = 9.81;

% - System specifications -
% Sensor sampling
dt_INS = 1/100;
dt_GNSS = 2;

% dt_INS = 2;
% dt_GNSS = 2;
N_between = dt_GNSS/dt_INS;

% - Noise specifications -
% Gyroscope
bias_W = -400 * (pi/180) * (1/3600);
bias_W_GM = 0;

std_W_GM = 0.01 * (pi/180);                             % std of underlying white noise
std_W_GM = 0.01 * (pi/180) * (1/sqrt(dt_INS));
std_W_RW = 0.1 * (pi/180) * (1/60) * (1/sqrt(dt_INS));  % std of underlying white noise

beta_W = 1/30;
% Accelerometer
bias_A_GM1 = -100 * (1e-6*g);
bias_A_GM2 = 200 * (1e-6*g);

std_A_WN = 50 * (1e-6*g) * (1/sqrt(dt_INS));
std_A_GM = 200 * (1e-6*g);
std_A_GM = 200 * (1e-6*g)* (1/sqrt(dt_INS));

beta_A = 1/60;

% GPS
std_GPS = 1;

%% Disturbance model

fm = @(az,fb) [cos(az) -sin(az);
              sin(az) cos(az)]*fb;

F = @(az,f) [[0 0 0 0 0 1 1 0 0;
               -f(2) 0 0 0 0 0 0 cos(az) -sin(az);
               f(1) 0 0 0 0 0 0 sin(az) cos(az);
               0 1 0 0 0 0 0 0 0;
               0 0 1 0 0 0 0 0 0]; ...
               zeros(4,5) diag([0 -beta_W -beta_A -beta_A])];
G = @(az) [[1 0 0;
           0 cos(az) -sin(az);
           0 sin(az) cos(az);
           0 0 0;
           0 0 0] zeros(5,3);
           zeros(4,3) ...
           [0 0 0;
            1 0 0;
            0 1 0;
            0 0 1]];

%% Kalman filter specification

std_W_GM_WN = sqrt( (1 - exp(-2*beta_W*dt_INS)) .* std_W_GM.^2 );
std_A_GM_WN = sqrt( (1 - exp(-2*beta_A*dt_INS)) .* std_A_GM.^2 );

W = diag([std_W_RW.^2 std_A_WN.^2 std_A_WN.^2 ...
        2*std_W_GM_WN.^2*beta_W 2*std_A_GM_WN.^2*beta_A 2*std_A_GM_WN.^2*beta_A]);
H = [0 0 0 1 0 0 0 0 0; 0 0 0 0 1 0 0 0 0];

R = diag([std_GPS.^2 std_GPS.^2]);

%% Simulation
% - Reference trajectory -              % Should depend on first gyro measurement!
P0_ref = [r,0];
V0_ref = [0,w*r];
R0_ref = [0 -1;1 0];
O0_ref = [0 -w;w 0];

% INS
phi_ref_INS = 0:w*dt_INS:2*pi;
N_INS = length(phi_ref_INS);
P_ref_INS = pos_sim(phi_ref_INS,r);
V_ref_INS = vel_sim(phi_ref_INS,r,w);

% GNSS
phi_ref_GNSS = 0:w*dt_GNSS:2*pi;
N_GNSS = length(phi_ref_GNSS);
P_ref_GNSS = pos_sim(phi_ref_GNSS,r);
V_ref_GNSS = vel_sim(phi_ref_GNSS,r,w);

% - Gyroscope measurements -
W_std = [0 std_W_RW std_W_GM];
W_bias = [bias_W 0 0 bias_W_GM];
W_noise_toggle = [1 0 1 1];
W_sim = gyro_sim([1, N_INS],w,W_noise_toggle,W_std,W_bias,beta_W,dt_INS);

% - Accelerometer measurements -
A_std = [std_A_WN std_A_WN; 0 0;std_A_GM std_A_GM];
A_bias = [0 0; 0 0; bias_A_GM1 bias_A_GM2];
A_noise_toggle = [1 0 1];
A_sim = acc_sim([2, N_INS],r,w,A_noise_toggle,A_std,A_bias,beta_A,dt_INS);

% - GPS measurements - 
GPS_sim = P_ref_GNSS + std_GPS.* randn(2,N_GNSS);


%% INS for comparison
% - Initialize vectors -
R_INS = zeros(2,2,N_INS);
V_INS = zeros(2,N_INS);
P_INS = zeros(2,N_INS);

R_INS(:,:,1) = R0_ref;
V_INS(:,1) = V0_ref;
P_INS(:,1) = P0_ref;

% - Iterate INS -
for i = 2:N_INS
    R_INS(:,:,i) = R_INS(:,:,i-1) * expm(omg_sim(W_sim(:,i)).*dt_INS);
    V_INS(:,i) = V_INS(:,i-1) + 1/2 * (R_INS(:,:,i)*A_sim(:,i)+R_INS(:,:,i-1)*A_sim(:,i-1))*dt_INS;
    P_INS(:,i) = P_INS(:,i-1) + 1/2 * (V_INS(:,i)+V_INS(:,i-1))*dt_INS;
end
%% INS/GNSS

% - Initialize vectors -
% Error model
E = zeros(4,N_GNSS);
dx_hat = zeros(9,N_GNSS);
W_adj = zeros(1,N_INS);
A_adj = zeros(2,N_INS);
W_adj(1) = W_sim(1);
A_adj(:,1) = A_sim(:,1);

% Internal
X_IN_GN = zeros(5,N_INS);
dX_0 = [3*pi/180,-2,-1,GPS_sim(1,1)-r,GPS_sim(2,1)]';
X_IN_GN(:,1) = [pi/2, 0,w*r,r,0]' + dX_0;
R_IN_GN = zeros(2,2,N_INS);
R_IN_GN(:,:,1) = rotMat(X_IN_GN(1,1));
X_IN_GN_iter = zeros(5,N_GNSS);

% Kalman
P_hat = zeros(9,9,N_GNSS);
P_hat(:,:,1) = diag([2*pi/180 5 5 ...
                    10 10 0.05*pi/180 ...
                    0.01*pi/180 3e-4*g 3e-4*g].^2);
P_iter = P_hat(:,:,1);              % Most recent update of P_tile
dx_iter = zeros(9,1);               % Most recent update of dx_tilde


% - Navigation Solution - 
iter = 1;
% Adjust = 1;
for i = 2:N_GNSS
    for j = 1:N_between
        % Increment iterable
        iter = iter +1;

        % Adjust Sensor readings
        W_adj(:,iter) = W_sim(:,iter) + ( E(1,i-1) + E(2,i-1) );
        A_adj(:,iter) = A_sim(:,iter) +  E(3:4,i-1);

        % Calculate INS
        R_IN_GN(:,:,iter) = R_IN_GN(:,:,iter-1)*expm(omg_sim(W_adj(:,iter).*dt_INS));
        X_IN_GN(1,iter) = X_IN_GN(1,iter-1) + 1/2 *(W_adj(:,iter)+W_adj(:,iter-1))*dt_INS;
        X_IN_GN(2:3,iter) = X_IN_GN(2:3,iter-1) + ...
                        + 1/2 * (R_IN_GN(:,:,iter)*A_adj(:,iter)+R_IN_GN(:,:,iter-1)*A_adj(:,iter-1))*dt_INS;
        X_IN_GN(4:5,iter) = X_IN_GN(4:5,iter-1) + 1/2 * (X_IN_GN(2:3,iter)+X_IN_GN(2:3,iter-1))*dt_INS;

        % Kalman prediction
        fm_iter = fm(X_IN_GN(1,iter),A_adj(:,iter));
        F_iter = F(X_IN_GN(1,iter),fm_iter);
        G_iter = G(X_IN_GN(1,iter));
        [Phi_iter,Q_iter] = Cont2Disc(F_iter,G_iter,W,dt_INS);
        dx_iter = Phi_iter * dx_iter;
        
        P_iter = Phi_iter * P_iter * Phi_iter' + Q_iter;
    end
    
    % Kalman gain
    K = P_iter * H' * (H*P_iter*H' + R)^(-1);
    
    % Kalman update
    X_IN_GN_iter(:,i) = X_IN_GN(:,iter);

    dz = GPS_sim(:,i)-X_IN_GN(4:5,iter);
    dx_hat(:,i) = dx_iter + K*(dz - H*dx_iter);
    P_hat(:,:,i) = (eye(9) - K*H)*P_iter;

    % Update Error and INS
    E(:,i) = E(:,i-1) + dx_hat(6:9,i);
    X_IN_GN(:,iter) = X_IN_GN(:,iter) + dx_hat(1:5,i);

    % Update iterables
    P_iter = P_hat(:,:,i);
    dx_iter= zeros(9,1);                   % Since in the EKF we are linearizing the current state
    
end

%% Plots 
%% Azimuth
figure
plot(phi_ref_INS+pi/2);hold on;
plot(X_IN_GN(1,:))
legend('Reference','Simulation')
title('Azimuth simulation')

%% Veloity
figure
plot(V_ref_INS(2,:),V_ref_INS(1,:));hold on;
plot(X_IN_GN(3,:),X_IN_GN(2,:))

%% Postition
figure
plot(P_ref_INS(2,:),P_ref_INS(1,:),'r');hold on
plot(GPS_sim(2,:),GPS_sim(1,:),'.k')
plot(X_IN_GN(5,:),X_IN_GN(4,:),'g','LineWidth',1)
plot(P_INS(2,:),P_INS(1,:))
plot(X_IN_GN_iter(5,:),X_IN_GN_iter(4,:),'.b','LineWidth',1)
axis([-700 700 -700 700])
legend('Reference','GPS','INS/GNSS','INS')
title('Trajectory simulation')
%% Gyro
figure
plot(ones(1,N_INS)*w);hold on
plot(W_sim);
plot(W_adj);
legend('Truth,','Simulation','Adjusted')

%% Acc
figure
subplot(2,1,1)
plot(zeros(1,N_INS));hold on;
plot(A_sim(1,:));
plot(A_adj(1,:));
legend('Truth,','Simulation','Adjusted')

subplot(2,1,2)
plot(r*w.^2*ones(1,N_INS));hold on;
plot(A_sim(2,:));
plot(A_adj(2,:));
legend('Truth,','Simulation','Adjusted')

%% Functions 
function P = pos_sim(Az,r)
    P = r*[cos(Az);sin(Az)];
end

function V  = vel_sim(Az,r,w)
    V = w*r*[-sin(Az);cos(Az)];
end

function Rmb = rotMat(az);
    Rmb = [cos(az) -sin(az);sin(az) cos(az)];
end

function O = omg_sim(gyro)
    O = [0 -gyro;gyro 0];
end

function W = gyro_sim(DIM,w,noise,std,bias,beta,dt)
    % noise - binary vector determining added noise
    % (1) - white noise
    % (2) - random walk 
    % (3) - gauss markov

    W = w.*ones(DIM);
    if nargin > 2
        B  = bias(1).*ones(DIM);
        WN = std(1) * randn(DIM) + bias(2);
        RW = randomWalk(DIM,std(2),bias(3));
        GM = gaussMarkov(DIM,std(3),bias(4),beta,dt);
        W = W + sum(noise'.*[B;WN; RW; GM]);
    end
end

function A = acc_sim(DIM,r,w,noise,std,bias,beta,dt)
    % noise - binary vector determining added noise
    % (1) - white noise
    % (2) - random walk 
    % (3) - gauss markov

    A = r*w.^2.*[zeros(1,DIM(2));ones(1,DIM(2))];
    if nargin > 3 
        WN = diag(std(1,:))*randn(DIM);
        RW = randomWalk(DIM,std(2,:),bias(2,:));
        GM = gaussMarkov(DIM,std(3,:),bias(3,:),beta,dt);

        A = A + noise(1).*WN + noise(2).*RW + noise(3).*GM;
    end

end

function RW = randomWalk(DIM,std,bias)
    % DIM - dimension of simulated signal
    % std - standard deviation of underlying white noise
    % bias - bias term of RW

    WN = diag(std)*randn(DIM);
    RW = zeros(DIM);
    if bias ~= 0
        RW(:,1) = bias';
    else
        RW(:,1) = WN(:,1);
    end
    
    for i = 2:DIM(2)
        RW(:,i) = RW(:,i-1) + WN(:,i-1);
    end
end

function GM = gaussMarkov(DIM,std,bias,beta,dt)
    % DIM - dimension of signal
    % beta - correlation coefficient
    % std - standard deviation of underlying white noise
    % bias - bias term of RW

    std_GM = sqrt( (1 - exp(-2*beta*dt)) .* std.^2 );

    WN = diag(std_GM) * randn(DIM);
    GM = zeros(DIM);

    if bias ~= 0
        GM(:,1) = bias';
    else
        GM(:,1) = WN(:,1);
    end
    
    for i = 2:DIM(2)
        GM(:,i) = exp(-beta.*dt).*GM(:,i-1)+WN(:,i-1);
    end
end 

function [Phi,Q] = Cont2Disc(F,G,W,dt)
    A = [-F G*W*G'; zeros(size(F)) F'].*dt;
    B = expm(A);
    B22 = B(10:18,10:18);
    B12 = B(1:9,10:18);
    Phi = B22';
    Q = B22'*B12;
    
end