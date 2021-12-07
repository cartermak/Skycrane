%% Data needed 
% x_err = xk - x_kp
% y_err = yk - H*x_km
% Sk = H*P_km*H' + R
clear;close all;clc;

% Definitions
dt = 0.1;
N = 500;
dx0 = [0; 0.1; 1; 0; 0; 0.001];
Q = [10,0,0;
     0,100,0;
     0,0,1000];
R = [250,0,0,0;
    0,2.5,0,0;
    0,0,2.5e-04,0;
    0,0,0,0.00225];

t = N+1;
NTMT = 50;
NEES_data = zeros(NTMT,t);
NIS_data  = zeros(NTMT,t);


for i = 1:NTMT % for each round of testing
    
    

% Instantiate system and simulate truth data
sys = SkycraneSystem(dt,N,dx0);

% Run Kalman filter
[x,P,P_pri,dx_Pri] = run_lkf(sys,Q,R);
dx_Pri = dx_Pri + sys.x_noms;
t = sys.ts;
x_filter = x;
P_filter = P;
x_truth = sys.xs;
y_truth = sys.ys;
[~,~,~,H] = sys.get_lin_matrices(0);
    
    
    
    for j = 1:N % for each time step
    x_err(:,j) = x_truth(:,j+1) - x(:,j+1);
    P_kp = P(:,:,j+1);
    NEES_data(i,j)  = NEES(x_err(:,j), P_kp);
    
    P_km = P_pri(:,:,j+1);
    
    y_err(:,j) = y_truth(:,j+1) - H*dx_Pri(:,j+1);
    Sk = H*P_km*H' + R;
    NIS_data(i,j) = NIS(y_err(:,j), Sk);
    
    
    end
end

alpha = 0.05;

% Mean across each test
exbar = mean(NEES_data,1);
eybar = mean(NIS_data,1);

% N = number of rounds of testing
% n,p = degrees of freedom
n = 4;
p = 2;
r1NEES = chi2inv(alpha/2, NTMT*n)/NTMT;
r2NEES = chi2inv(1-alpha/2, NTMT*n)/NTMT;

r1NIS = chi2inv(alpha/2,NTMT*p)/NTMT;
r2NIS = chi2inv(1-alpha/2,NTMT*p)/NTMT;



%% Plots

% NEES Plot


figure()
hold on
plot(r1NEES*ones(1,N),'-r')
plot(r2NEES*ones(1,N),'-r')
plot(exbar,'.')
xlabel('Time')
ylabel('Average NEES $\bar{\epsilon}_x$','Interpreter','latex')
title('Average NEES Statistis for TMT')
hold off


% NIS Plot
figure()
hold on
plot(eybar,'.')
plot(r1NIS*ones(1,N),'-r')
plot(r2NIS*ones(1,N),'-r')
xlabel('Time')
ylabel('Average NIS $\bar{\epsilon}_y$','Interpreter','latex')
title('Average NIS Statistis for TMT')

hold off