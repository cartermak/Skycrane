clearvars; close all; clc;

% Definitions
dt = 0.1;
N = 500;
dx0 = [0; 0.1; 5; 0; 0; 0.001];
Q = [10,5,0;    % xi_ddot
    5,100,0;    % z_ddot
     0,0,1]; % theta_ddot
R = [25,0,0,0;       % xi
    0,5,0,0;         % z
    0,0,2.5e-05,0;   % theta_dot
    0,0,0,0.000225]; % xi_ddot

% Instantiate system and simulate truth data
sys = SkycraneSystem(dt,N,dx0);

%% Run Linearized Kalman filter
[x,P] = run_lkf(sys,Q,R);

% Plot results
x_labels = ["$\xi$"; "$\dot{\xi}$"; "$z$";
            "$\dot{z}$"; "$\theta$"; "$\dot{\theta}$"];
x_units = ["m"; "m/s"; "m";
            "m/s"; "rad"; "rad/s"];
name = "Linearized Kalman Filter";
plot_filter_performance(x,P,sys.xs,sys.ts,x_labels,x_units,name);

t = sys.ts;
x_filter = x;
P_filter = P;
x_truth = sys.xs;
y_truth = sys.ys;
[~,~,~,H] = sys.get_lin_matrices(0);

save('lkf_run.mat','x_filter','P_filter','x_truth','y_truth','H','t')

%% Run Extended Kalman filter
[x,P] = run_ekf(sys,Q,R);

% Plot results
x_labels = ["$\xi$"; "$\dot{\xi}$"; "$z$";
            "$\dot{z}$"; "$\theta$"; "$\dot{\theta}$"];
x_units = ["m"; "m/s"; "m";
            "m/s"; "rad"; "rad/s"];
name = "Extended Kalman Filter";
plot_filter_performance(x,P,sys.xs,sys.ts,x_labels,x_units,name);

t = sys.ts;
x_filter = x;
P_filter = P;
x_truth = sys.xs;
y_truth = sys.ys;
[~,~,~,H] = sys.get_lin_matrices(0);

save('lkf_run.mat','x_filter','P_filter','x_truth','y_truth','H','t')