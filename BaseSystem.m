classdef (Abstract) BaseSystem
    %BaseSystem Base class for dynamical systems 
    % Manage access to data and EOM for a dynamical system as it pertains
    % to filter design and tuning.

    properties(Access=public)
        f % Nonlinear dynamics equations. f(x,u)
        h % Nonlinear measurement equations. h(x,u)
%         K_ctrl % Feedback control gain matrix
        A % A(x,u)
        B % B(x,u)
        Gamma % static matrix
        C % C(x,u)
        D % D(x,u)
        Q % Process noise covariance matrix
        R % Measurement noise covariance matrix
        
        dt
        N
        ts
        x_noms
        u_noms
        y_noms
        xs
        us
        ys
        dxs
        dus
        dys
        n
        m
        p
    end

    methods(Access=public,Abstract)
        generate_data(obj,dx0)
    end
    
    methods(Access=public)
        
        function obj = BaseSystem(dt,N,dx0)
            obj.dt = dt;
            obj.N = N;
            obj.ts = linspace(0,N*dt,N+1);
            obj = obj.generate_data(dx0);
        end

        function [F,G,Omega,H,M] = get_lin_matrices(obj,k)
            %get_lin_matrices Get SS matrices for Linearized KF
            % Get the linearized DT SS matrices for use with the linearized
            % Kalman filter.
            
            x_nom = obj.x_noms(:,k+1);
            u_nom = obj.u_noms(:,k+1);
            
            F = eye(obj.n) + obj.dt*obj.A(x_nom,u_nom);
            G = obj.dt*obj.B(x_nom,u_nom);
            Omega = obj.dt*obj.Gamma;
            H = obj.C(x_nom,u_nom);
            M = obj.D(x_nom,u_nom);
            
            F(isnan(F)) = 0;
            G(isnan(G)) = 0;
            H(isnan(H)) = 0;
            M(isnan(M)) = 0;
            
        end
        
        function [F,G,Omega,H,M] = get_nl_matrices(obj,x_k,u_k)
            %get_nl_matrices Get SS matrices linearized about given state
            % This method evaluates the CT Jacobians about the given state
            % x_k at timestep k and converts to DT using the first-order
            % matrix Euler approximation over the configured timestep.
            
            F = eye(obj.n) + obj.dt*obj.A(x_k,u_k);
            G = obj.dt*obj.B(x_k,u_k);
            Omega = obj.dt*obj.Gamma;
            H = obj.C(x_k,u_k);
            M = obj.D(x_k,u_k);
            
        end
        
        function x_kp1 = integrate_nl_dynamics(obj,x_k,u_k)
            %integrate_nl_dynamics Propagate full nonlinear dynamics
            % Numerically integrates from the given state x_k to the next
            % timestep x_{k+1} using the full nonlinear dynamics and a
            % single RK4 iteration using the configured timestep.
            
%             k1 = obj.f(x_k,u_k);
%             k2 = obj.f(x_k+obj.dt*k1/2,u_k);
%             k3 = obj.f(x_k+obj.dt*k2/2,u_k);
%             k4 = obj.f(x_k+obj.dt*k3,u_k);
%             
%             x_kp1 = x_k + (obj.dt/6)*(k1+2*k2+2*k3+k4);

            [~,x] = ode45(@(~,x)obj.f(x,u_k),[0,obj.dt],x_k);
            
            x_kp1 = x(end,:)';
            
        end
        
        function yhat_k = predict_nl_measurement(obj,x_k,u_k)
            %predict_nl_measurement Nonlinear expected meas. based on state
            % Evaluate the full nonlinear measurement equations on the
            % given state estimate in order to predict the expected
            % measurement vector.
            
            yhat_k = obj.h(x_k,u_k);
            
        end
        
        function du_k = get_ctrl_perturbation(obj,k)
            %get_ctrl_perturbation Get control perturb. from nominal at k
            % Return the control perturbation from the nominal solution at
            % time k, du_k. 
            
            du_k = obj.dus(:,k+1);
            
        end
        function u_k = get_ctrl(obj,k)
            %get_ctrl Get control input at time k
            % Return the full control input at timestep k, u_k.
            
            u_k = obj.us(:,k+1);
            
        end
        function dy_k = get_meas_perturbation(obj,k)
            %get_meas_perturbation Get measurement perturb. from nominal
            % Return the difference between the measurement and the nominal
            % measurement at time k.
            
            dy_k = obj.dys(:,k+1);
            
        end
        function y_k = get_meas(obj,k)
            %get_meas Get absolute measurement at time k
            % Return the full measurement vector at time k.
            
            y_k = obj.ys(:,k+1);
            
        end
    end
end