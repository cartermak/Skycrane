function [x,P,P_pri,dx_pri] = run_lkf(sys,Q,R)
%run_lkf Run Linearized Kalman Filter given input parameters
% sys [BaseSystem]: Dynamical system to filter
% Q [n-by-n matrix]: Estimate of process noise covariance
% R [p-by-p matrix]: Estimate of measurement noise covariance

% Get system dimensions
n = sys.n;
N = sys.N;

% Common values
I = eye(n);

% Preallocate and initialize output matrices
dx = zeros(n,N+1);
dx_pri = zeros(n,N+1);
P = zeros(n,n,N+1);
P_pri = zeros(n,n,N+1);
P(:,:,1) = eye(n);

% Initialize terms before first iteration
[F,G,Omega,~,~] = sys.get_lin_matrices(0);
du = sys.get_ctrl_perturbation(0);

for k = 1:N
    
    % Propagate previous state estimate with previous control inputs
    dx_pri = F*dx(:,k) + G*du;
    
    % Update to get current control perturbation
    du = sys.get_ctrl_perturbation(k);
    
    % Propagate previous state cov. through dynamics with process noise
    P_pri(:,:,k+1) = F*P(:,:,k)*F' + Omega*Q*Omega';
    
    % Get current measurement perturbation
    dy = sys.get_meas_perturbation(k);
    
    % Update DT SS matrices to current timestep
    [F,G,Omega,H,M] = sys.get_lin_matrices(k);
    
    % Calculate Kalman gain
    K = P_pri(:,:,k+1)*H'/(H*P_pri(:,:,k+1)*H' + R);
    
    % Innovation vector
    dy_hat = dy-H*dx_pri-M*du;
    
    % Correct state estimate with measurement
    dx(:,k+1) = dx_pri + K*(dy_hat);
    
    % Update covariance to consider measurement
    P(:,:,k+1) = (I - K*H)*P_pri(:,:,k+1);
end

x = dx + sys.x_noms;

end