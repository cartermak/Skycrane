classdef SkycraneSystem < BaseSystem
    methods
        
        function obj = SkycraneSystem(dt,N,dx0)
            obj = obj@BaseSystem(dt,N,dx0);
        end

        function obj = generate_data(obj,dx0)

            syms xi xi_dot
            syms z z_dot
            syms theta theta_dot
            syms T_1 T_2
            syms t
            assume(xi_dot,'real')
            assume(z_dot,'real')
            
            rho = 0.02;
            g = 3.711;
            beta = pi/4;
            C_D = 0.2;
            m_f  = 390;
            m_b = 1510;
            w_b  = 3.2;
            h_cm = 0.9421;
            A_side = 7.75;
            A_bot = 10.28;
            I_eta = (1/12)*(1510*(3.2^2+2.5^2) + 390*(1+0.25));
            w_tilde_1 = 0;
            w_tilde_2 = 0;
            w_tilde_3 = 0;

            x = [xi;xi_dot;z;z_dot;theta;theta_dot];
            u = [T_1;T_2];

            alpha = atan2(z_dot,xi_dot);

            F_D_xi = (1/2)*C_D*rho*(A_side*cos(theta-alpha) + ...
                A_bot*sin(theta-alpha))*xi_dot*sqrt(xi_dot^2 + z_dot^2);

            F_D_z = (1/2)*C_D*rho*(A_side*cos(theta-alpha) + ...
                A_bot*sin(theta-alpha))*z_dot*sqrt(xi_dot^2 + z_dot^2);

            xi_ddot = (T_1*(cos(beta)*sin(theta)+sin(beta)*cos(theta)) +...
                T_2*(cos(beta)*sin(theta)-sin(beta)*cos(theta)) - ...
                F_D_xi)/(m_b+m_f) + w_tilde_1;

            z_ddot = (T_1*(cos(beta)*cos(theta)-sin(beta)*sin(theta)) + ...
                T_2*(cos(beta)*cos(theta)+sin(beta)*sin(theta)) - ...
                F_D_z)/(m_b+m_f) - g + w_tilde_2;

            theta_ddot = (1/I_eta)*((T_1 - T_2)*cos(beta)*w_b/2 + ...
                (T_2 - T_1)*sin(beta)*h_cm) + w_tilde_3;
            
            % Nonlinear dynamics function
            f = [xi_dot;xi_ddot;
                z_dot;z_ddot;
                theta_dot;theta_ddot];
            df_dx = jacobian(f,x);
            df_du = jacobian(f,u);
            
            % Measurement function
            h = [xi;z;theta_dot;xi_ddot];
            dh_dx = jacobian(h,x);
%             dh_du = jacobian(h,u);

            % Function handles for later evaluations
            obj.f = matlabFunction(f,'Vars',{x,u});
            obj.h = matlabFunction(h,'Vars',{x,u});
            obj.A = matlabFunction(df_dx,'Vars',{x,u});
            obj.B = matlabFunction(df_du,'Vars',{x,u});
            obj.C = matlabFunction(dh_dx,'Vars',{x,u});
            obj.Gamma = [0,0,0;1,0,0;0,0,0;0,1,0;0,0,0;0,0,1];
            
            % Matrix dimensions
            obj.n = 6;
            obj.m = 3;
            obj.p = 4;
            
            % Static nominal solution
            x_nom = [0; 0; 20; 0; 0; 0];
            u_nom = [0.5*3.711*(1510+390)/cos(pi/4);...
                     0.5*3.711*(1510+390)/cos(pi/4)];
            y_nom = obj.h(x_nom,u_nom);
            obj.x_noms = repmat(x_nom,1,obj.N+1);
            obj.u_noms = repmat(u_nom,1,obj.N+1);
            obj.y_noms = repmat(y_nom,1,obj.N+1);
            
            % Evaluate mapping from process noise to DT state space
            Omega = obj.dt*obj.Gamma;
            
            % Load given matrices for control and noise
            K_ctrl = 1e-7*load('skycrane_finalproj_KFdata.mat','Klin').Klin;
            Q = 1e-7*load('skycrane_finalproj_KFdata.mat','Qtrue').Qtrue;
            R = load('skycrane_finalproj_KFdata.mat','Rtrue').Rtrue;
            
            % Matrix square roots via Cholesky decomp. for noise
            S_w = chol(Q,'lower');
            S_v = chol(R,'lower');
            
            % Preallocate output matrices
            obj.xs = zeros(6,obj.N+1);
            obj.us = zeros(2,obj.N+1);
            obj.ys = zeros(4,obj.N+1);
            
            % Initialize state and control
            obj.xs(:,1) = dx0 + x_nom;
            obj.us(:,1) = u_nom - K_ctrl*(obj.xs(:,1)-x_nom);
            
            for k = 1:obj.N
                % Propagate state from k-1 to k with simulated noise
                obj.xs(:,k+1) = obj.integrate_nl_dynamics(...
                    obj.xs(:,k),obj.us(:,k)) + Omega*S_w*randn(3,1);
                
                % Evalulate new control input at time k
                obj.us(:,k+1) = u_nom-K_ctrl*(obj.xs(:,k+1)-x_nom);
                
                % Evaluate measurement with simulated noise
                obj.ys(:,k+1) = obj.h(obj.xs(:,k+1),obj.us(:,k+1)) + ...
                    S_v*randn(4,1);
            end
            
            % Calculate delta- terms for x, u, and y
            obj.dxs = obj.xs - obj.x_noms;
            obj.dus = obj.us - obj.u_noms;
            obj.dys = obj.ys - obj.y_noms;
            
        end
        
    end
end