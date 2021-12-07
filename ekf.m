function [x,P,P_pri,x_pri] = ekf(sys,Q,R)
%run_ekf Run Extended Kalman Filter given input parameters
% sys [BaseSystem]: Dynamical system to filter
% Q [n-by-n matrix]: Estimate of process noise covariance
% R [p-by-p matrix]: Estimate of measurement noise covariance

% Get system dimensions
n = sys.n;
N = sys.N;

% Common values
I = eye(n);

% Preallocate and initialize output matrices
x = sys.x_noms;
x_pri = zeros(n,N+1);
P = zeros(n,n,N+1);
P_pri = zeros(n,n,N+1);
P(:,:,1) = eye(n);

% Initialize terms before first iteration
[F,~,Omega,~,~] = sys.get_lin_matrices(0);
u = sys.get_ctrl(0);

for k = 1:N
    
    % Propagate previous state estimate with previous control inputs
    x_pri(:,k+1) = sys.integrate_nl_dynamics(x(:,k),u);
    
    % Update to get current control perturbation
    u = sys.get_ctrl(k);
    
    % Propagate previous state cov. through dynamics with process noise
    P_pri(:,:,k+1) = F*P(:,:,k)*F' + Omega*Q*Omega';
    
    % Get current measurement perturbation
    y = sys.get_meas(k);
    
    % Update DT SS matrices to current timestep
    [F,~,Omega,H,~] = sys.get_lin_matrices(k);
    
    % Calculate Kalman gain
    K = P_pri(:,:,k+1)*H'/(H*P_pri(:,:,k+1)*H' + R);
    
    % Innovation vector
    e_y = y - sys.h(x_pri(:,k+1),u);
    
    % Correct state estimate with measurement
    x(:,k+1) = x_pri(:,k+1) + K*(e_y);
    
    % Update covariance to consider measurement
    P(:,:,k+1) = (I - K*H)*P_pri(:,:,k+1);
end

end