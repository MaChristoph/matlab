%
%     This file is part of CasADi.
%
%     CasADi -- A symbolic framework for dynamic optimization.
%     Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
%                             K.U. Leuven. All rights reserved.
%     Copyright (C) 2011-2014 Greg Horn
%
%     CasADi is free software; you can redistribute it and/or
%     modify it under the terms of the GNU Lesser General Public
%     License as published by the Free Software Foundation; either
%     version 3 of the License, or (at your option) any later version.
%
%     CasADi is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%     Lesser General Public License for more details.
%
%     You should have received a copy of the GNU Lesser General Public
%     License along with CasADi; if not, write to the Free Software
%     Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
%

% An implementation of direct multiple shooting
% Joel Andersson, 2016

clc;
close all;


import casadi.*

T = 1.0; % Time horizon
N = 10; % number of control intervals

% Declare model variables
x1 = SX.sym('x1');
x2 = SX.sym('x2');
x3 = SX.sym('x3');
x4 = SX.sym('x4');
x5 = SX.sym('x5');
x6 = SX.sym('x6');
x7 = SX.sym('x7');
x8 = SX.sym('x8');
x9 = SX.sym('x9');
x10 = SX.sym('x10');
x = [x1; x2; x3; x4; x5; x6; x7; x8; x9; x10];

u1 = SX.sym('u1');
u2 = SX.sym('u2');
u3 = SX.sym('u3');
u4 = SX.sym('u4');
u5 = SX.sym('u5');
u6 = SX.sym('u6');
u7 = SX.sym('u7');
u8 = SX.sym('u8');
u9 = SX.sym('u9');
u10 = SX.sym('u10');
u = [u1; u2; u3; u4; u5; u6; u7; u8; u9; u10];

nu = 10;
nx = 10;

q0 = [0;0;0;0.01;0.01;0.01;0.01;0.01;0.01;0.01];  %% Due to singularity of the manipulator

x_des = [0.5;0.5;0.2];
err = 10000;

min_dist_link0_ee = 0.6;
min_dist_link0_5 = 0.05;
min_dist_link3_ee = 0.05;

error_old = 1;
% for kk = 1:22
ii = 1;

set = 1;
for set = 1:3
    if set == 1
        x_des = [0.5;0.5;0.2];
        x_rot_des = [0; 0; 0];
    elseif set == 2
        x_des = [-0.5;0.3;0.2];
        x_rot_des = [0; 0; 0];
    elseif set == 3
        x_des = [5;5;0.5];
        x_rot_des = [0; 0; 0];
    end

%     while err > 1e-4
    for ii = 1:50
        error_old = err;
        % Model equations
        xdot = [u1; u2; u3; u4; u5; u6; u7; u8; u9; u10];

        Q = 1*eye(3,3);
        Qr = 1*eye(3,3);
        R = 0.05*eye(nu,nu);
        %% Objective term

        % Add Barrier for distance
        d0_ee = (forward_kinematics(x))' * (forward_kinematics(x));
%         d0_5 = (forward_kinematics_link5(x))' * (forward_kinematics_link5(x));
%         d3_ee = (forward_kinematics_link3(x) -forward_kinematics(x))' * (forward_kinematics_link3(x) -forward_kinematics(x));
%         
% 
        dist_barrier = exp((min_dist_link0_ee - sqrt(d0_ee))/0.001);
%         dist_barrier2 = exp((min_dist_link0_5 - sqrt(d0_5))/0.001);
%         dist_barrier3 = exp((min_dist_link3_ee - sqrt(d3_ee))/0.001);

    %         dist_barrier = -log(sqrt(d0_ee) - min_dist_link0_ee);
        L = u' * R * u + ...
            (forward_kinematics(x)-x_des)' * Q * (forward_kinematics(x)-x_des) + ...
             dist_barrier;
% + dist_barrier2 + dist_barrier3; % + ...
%             (getRot(x) - x_rot_des)' * Qr * (getRot(x) - x_rot_des);


        % Continuous time dynamics
%         f = Function('f', {x, u}, {xdot, L});


        %% Euler discretize
        DT = T/N;
        f = Function('f', {x, u}, {xdot, L})
        X0 = MX.sym('X0', nx);
        U = MX.sym('U',nu);
        X = X0;
        Q = 0;

        [k1, k1_q] = f(X, U);      %[x, cost] = f(X,U)

        X=X+DT * k1;        % X_new
        Q = Q + DT * k1_q;  % Cost_new

        F = Function('F', {X0, U}, {X, Q}, {'x0','p'}, {'xf', 'qf'});
                          %Inout,  output,  Input(syms), output(syms)

%         dae = struct('x',x,'p',u,'ode',xdot,'quad',L);
%         opts = struct('tf',T/N);
%         F = integrator('F', 'cvodes', dae, opts);
        % Start with an empty NLP
        w={};
        w0 = [];
        lbw = [];
        ubw = [];
        J = 0;
        g={};
        lbg = [];
        ubg = [];

        % "Lift" initial conditions
        X0 = MX.sym('X0', nx);
        w = {w{:}, X0};
        lbw = [lbw; q0];
        ubw = [ubw; q0];
        w0 = [w0; q0];


        % Formulate the NLP
        Xk = X0;
        for k=0:N-1
            % New NLP variable for the control
            Uk = MX.sym(['U_' num2str(k)], nu);
            w = {w{:}, Uk};
            lbw = [lbw; -2; -2; -2; -2; -2; -2; -2;-2; -2; -2]; % 7 inputs, the first element is to continue the matrix
            ubw = [ubw; 2; 2; 2; 2; 2; 2; 2; 2; 2; 2];
            w0 = [w0;  0; 0; 0; 0; 0; 0; 0; 0; 0; 0];

            % Integrate till the end of the interval
            Fk = F('x0', Xk, 'p', Uk);
            Xk_end = Fk.xf;
            J=J+Fk.qf;

            % New NLP variable for state at end of interval
            Xk = MX.sym(['X_' num2str(k+1)], nx);
            w = {w{:}, Xk};
            lbw = [lbw; -inf;-inf; -pi; -4.712; -1.92; -2.967; -2.618; -2.967; -2.443; -2.967]; % 7 states, the first element is to continue the matrix
            ubw = [ubw;  inf; inf; pi; 4.712; 1.92;  2.967; 2.618; 2.967;  2.443; 2.967];
            w0 = [w0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];

            % Add equality constraint
            g = {g{:}, Xk_end-Xk};
            lbg = [lbg; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
            ubg = [ubg; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
        end

        % Create an NLP solver
        prob = struct('f', J, 'x', vertcat(w{:}), 'g', vertcat(g{:}));
        
        options = struct;
%         options.ipopt.hessian_approximation = 'limited-memory';
        solver = nlpsol('solver', 'ipopt', prob, options);


        % Solve the NLP
        sol = solver('x0', w0, 'lbx', lbw, 'ubx', ubw,...
                    'lbg', lbg, 'ubg', ubg);
        w_opt = full(sol.x);

        % Plot the solution
        n_var = nu + nx;
        x1_opt = w_opt(1:n_var:end);
        x2_opt = w_opt(2:n_var:end);
        x3_opt = w_opt(3:n_var:end);
        x4_opt = w_opt(4:n_var:end);
        x5_opt = w_opt(5:n_var:end);
        x6_opt = w_opt(6:n_var:end);
        x7_opt = w_opt(7:n_var:end);
        x8_opt = w_opt(8:n_var:end);
        x9_opt = w_opt(9:n_var:end);
        x10_opt = w_opt(10:n_var:end);


        u1_opt = w_opt(11:n_var:end);
        u2_opt = w_opt(12:n_var:end);
        u3_opt = w_opt(13:n_var:end);
        u4_opt = w_opt(14:n_var:end);
        u5_opt = w_opt(15:n_var:end);
        u6_opt = w_opt(16:n_var:end);
        u7_opt = w_opt(17:n_var:end);
        u8_opt = w_opt(18:n_var:end);
        u9_opt = w_opt(19:n_var:end);
        u10_opt = w_opt(20:n_var:end);

        q_opt = [x1_opt x2_opt x3_opt x4_opt x5_opt x6_opt x7_opt x8_opt x9_opt x10_opt];

        x_cartesian=[];

        for i=1:size(x1_opt)
        %     q_opt(1,:)
            x_cartesian = [x_cartesian forward_kinematics(q_opt(i,:))];
        %     q_opt(i,:)
        end

        hold all;
        plot3(x_cartesian(1,1), x_cartesian(2,1),  x_cartesian(3,1),'.')
        view(3)
        plot3(x_cartesian(1,:), x_cartesian(2,:),  x_cartesian(3,:),'g')

        axis([-5 5 -5 5 0 1])
        grid
        plot3(x_des(1),x_des(2),x_des(3) ,'x')
        q0 = [x1_opt(2); x2_opt(2); x3_opt(2); x4_opt(2); x5_opt(2); x6_opt(2); x7_opt(2); x8_opt(2); x9_opt(2); x10_opt(2)];  %% Optimal State after aplling optimal input

        plot3(x1_opt, x2_opt, zeros(size(x1_opt)), 'r')
        
                [xs,ys,zs] = sphere;
        a=0;
        b=0;
        c=0;
        surf(xs*min_dist_link0_ee+a, ys*min_dist_link0_ee+b, zs*min_dist_link0_ee+c,'FaceColor', [1 0 0],'FaceAlpha',1.0,'EdgeColor','none','LineStyle','none');
    
%         x_3 = forward_kinematics_link3(q0);
%         x_5 = forward_kinematics_link5(q0);
%         
%         a = x_3(1); b=x_3(2); c=x_3(3);
%         surf(xs*min_dist_link3_ee+a, ys*min_dist_link3_ee+b, zs*min_dist_link3_ee+c,'FaceColor', [1 0 0],'FaceAlpha',1.0,'EdgeColor','none','LineStyle','none');
%         
%         a = x_5(1); b=x_5(2); c=x_5(3);
%         surf(xs*min_dist_link0_5+a, ys*min_dist_link0_5+b, zs*min_dist_link0_5+c,'FaceColor', [1 0 0],'FaceAlpha',1.0,'EdgeColor','none','LineStyle','none');
% 
%          
%         err = (forward_kinematics(q0) - x_des)' * (forward_kinematics(q0) - x_des);
        forward_kinematics(q0)

        e_delta = err - error_old;
    %     if e_delta == 0
    %         disp('Break due error condition')
    %         break;
    %     end

%         if ii > 50
%             break;
%         end
        ii = ii + 1;
    end
end


function ee_euler = getRot(q)
    q1 = q(1);
    q2 = q(2);
    q3 = q(3);
    q4 = q(4);
    q5 = q(5);
    q6 = q(6);
    q7 = q(7);
rotm= [cos(q7)*(sin(q6)*(sin(q4)*(sin(q1)*sin(q3) - cos(q1)*cos(q2)*cos(q3)) - cos(q1)*cos(q4)*sin(q2)) - cos(q6)*(cos(q5)*(cos(q4)*(sin(q1)*sin(q3) - cos(q1)*cos(q2)*cos(q3)) + cos(q1)*sin(q2)*sin(q4)) + sin(q5)*(cos(q3)*sin(q1) + cos(q1)*cos(q2)*sin(q3)))) + sin(q7)*(sin(q5)*(cos(q4)*(sin(q1)*sin(q3) - cos(q1)*cos(q2)*cos(q3)) + cos(q1)*sin(q2)*sin(q4)) - cos(q5)*(cos(q3)*sin(q1) + cos(q1)*cos(q2)*sin(q3))), cos(q7)*(sin(q5)*(cos(q4)*(sin(q1)*sin(q3) - cos(q1)*cos(q2)*cos(q3)) + cos(q1)*sin(q2)*sin(q4)) - cos(q5)*(cos(q3)*sin(q1) + cos(q1)*cos(q2)*sin(q3))) - sin(q7)*(sin(q6)*(sin(q4)*(sin(q1)*sin(q3) - cos(q1)*cos(q2)*cos(q3)) - cos(q1)*cos(q4)*sin(q2)) - cos(q6)*(cos(q5)*(cos(q4)*(sin(q1)*sin(q3) - cos(q1)*cos(q2)*cos(q3)) + cos(q1)*sin(q2)*sin(q4)) + sin(q5)*(cos(q3)*sin(q1) + cos(q1)*cos(q2)*sin(q3)))),   cos(q6)*(sin(q4)*(sin(q1)*sin(q3) - cos(q1)*cos(q2)*cos(q3)) - cos(q1)*cos(q4)*sin(q2)) + sin(q6)*(cos(q5)*(cos(q4)*(sin(q1)*sin(q3) - cos(q1)*cos(q2)*cos(q3)) + cos(q1)*sin(q2)*sin(q4)) + sin(q5)*(cos(q3)*sin(q1) + cos(q1)*cos(q2)*sin(q3)));
 -cos(q7)*(sin(q6)*(sin(q4)*(cos(q1)*sin(q3) + cos(q2)*cos(q3)*sin(q1)) + cos(q4)*sin(q1)*sin(q2)) - cos(q6)*(cos(q5)*(cos(q4)*(cos(q1)*sin(q3) + cos(q2)*cos(q3)*sin(q1)) - sin(q1)*sin(q2)*sin(q4)) + sin(q5)*(cos(q1)*cos(q3) - cos(q2)*sin(q1)*sin(q3)))) - sin(q7)*(sin(q5)*(cos(q4)*(cos(q1)*sin(q3) + cos(q2)*cos(q3)*sin(q1)) - sin(q1)*sin(q2)*sin(q4)) - cos(q5)*(cos(q1)*cos(q3) - cos(q2)*sin(q1)*sin(q3))), sin(q7)*(sin(q6)*(sin(q4)*(cos(q1)*sin(q3) + cos(q2)*cos(q3)*sin(q1)) + cos(q4)*sin(q1)*sin(q2)) - cos(q6)*(cos(q5)*(cos(q4)*(cos(q1)*sin(q3) + cos(q2)*cos(q3)*sin(q1)) - sin(q1)*sin(q2)*sin(q4)) + sin(q5)*(cos(q1)*cos(q3) - cos(q2)*sin(q1)*sin(q3)))) - cos(q7)*(sin(q5)*(cos(q4)*(cos(q1)*sin(q3) + cos(q2)*cos(q3)*sin(q1)) - sin(q1)*sin(q2)*sin(q4)) - cos(q5)*(cos(q1)*cos(q3) - cos(q2)*sin(q1)*sin(q3))), - cos(q6)*(sin(q4)*(cos(q1)*sin(q3) + cos(q2)*cos(q3)*sin(q1)) + cos(q4)*sin(q1)*sin(q2)) - sin(q6)*(cos(q5)*(cos(q4)*(cos(q1)*sin(q3) + cos(q2)*cos(q3)*sin(q1)) - sin(q1)*sin(q2)*sin(q4)) + sin(q5)*(cos(q1)*cos(q3) - cos(q2)*sin(q1)*sin(q3)));
 cos(q7)*(cos(q6)*(cos(q5)*(cos(q2)*sin(q4) + cos(q3)*cos(q4)*sin(q2)) - sin(q2)*sin(q3)*sin(q5)) + sin(q6)*(cos(q2)*cos(q4) - cos(q3)*sin(q2)*sin(q4))) - sin(q7)*(sin(q5)*(cos(q2)*sin(q4) + cos(q3)*cos(q4)*sin(q2)) + cos(q5)*sin(q2)*sin(q3)),                                                                                                                                                                   - cos(q7)*(sin(q5)*(cos(q2)*sin(q4) + cos(q3)*cos(q4)*sin(q2)) + cos(q5)*sin(q2)*sin(q3)) - sin(q7)*(cos(q6)*(cos(q5)*(cos(q2)*sin(q4) + cos(q3)*cos(q4)*sin(q2)) - sin(q2)*sin(q3)*sin(q5)) + sin(q6)*(cos(q2)*cos(q4) - cos(q3)*sin(q2)*sin(q4))),                                                                                                       cos(q6)*(cos(q2)*cos(q4) - cos(q3)*sin(q2)*sin(q4)) - sin(q6)*(cos(q5)*(cos(q2)*sin(q4) + cos(q3)*cos(q4)*sin(q2)) - sin(q2)*sin(q3)*sin(q5))];

x_e = atan2(-rotm(3,1), rotm(1,1));
y_e = asin(rotm(2,1));
z_e = atan2(-rotm(2,3),rotm(2,2));

ee_euler = [x_e; y_e; z_e];

end


function x_ee = forward_kinematics(q)
    q1 = q(1);
    q2 = q(2);
    q3 = q(3);
    q4 = q(4);
    q5 = q(5);
    q6 = q(6);
    q7 = q(7);
    q8 = q(8);
    q9 = q(9);
    q10 = q(10);
    
%     x_ee =  [(1493*sin(q4)*(sin(q1)*sin(q3) - cos(q1)*cos(q2)*cos(q3)))/5000 - (567*cos(q1)*sin(q2))/1250 - (1493*cos(q1)*cos(q4)*sin(q2))/5000;
%              -(567*sin(q1)*sin(q2))/1250 - (1493*sin(q4)*(cos(q1)*sin(q3) + cos(q2)*cos(q3)*sin(q1)))/5000 - (1493*cos(q4)*sin(q1)*sin(q2))/5000;
%                (567*cos(q2))/1250 + (1493*cos(q2)*cos(q4))/5000 - (1493*cos(q3)*sin(q2)*sin(q4))/5000 + 563/2000];

    x_ee =  [q1 + (1493*sin(q7)*(sin(q3 + q4)*sin(q6) - cos(q3 + q4)*cos(q5)*cos(q6)))/5000 - (567*cos(q3 + q4)*sin(q5))/1250 - (1493*cos(q3 + q4)*cos(q7)*sin(q5))/5000;
             q2 - (1493*sin(q7)*(cos(q3 + q4)*sin(q6) + sin(q3 + q4)*cos(q5)*cos(q6)))/5000 - (567*sin(q3 + q4)*sin(q5))/1250 - (1493*sin(q3 + q4)*cos(q7)*sin(q5))/5000;
              (567*cos(q5))/1250 + (1493*cos(q5)*cos(q7))/5000 - (1493*cos(q6)*sin(q5)*sin(q7))/5000 + 563/2000];
 

end


function x = forward_kinematics_link3(q)
    q1 = q(1);
    q2 = q(2);
    q3 = q(3);
    q4 = q(4);
    q5 = q(5);
    q6 = q(6);
    q7 = q(7);
  x =  [-(567*cos(q1)*sin(q2))/1250;
           -(567*sin(q1)*sin(q2))/1250;
         (567*cos(q2))/1250 + 563/2000];
end


function x = forward_kinematics_link5(q)
    q1 = q(1);
    q2 = q(2);
    q3 = q(3);
    q4 = q(4);
    q5 = q(5);
    q6 = q(6);
    q7 = q(7);
  x =  [ (1493*sin(q4)*(sin(q1)*sin(q3) - cos(q1)*cos(q2)*cos(q3)))/5000 - (567*cos(q1)*sin(q2))/1250 - (1493*cos(q1)*cos(q4)*sin(q2))/5000;
         -(567*sin(q1)*sin(q2))/1250 - (1493*sin(q4)*(cos(q1)*sin(q3) + cos(q2)*cos(q3)*sin(q1)))/5000 - (1493*cos(q4)*sin(q1)*sin(q2))/5000;
         (567*cos(q2))/1250 + (1493*cos(q2)*cos(q4))/5000 - (1493*cos(q3)*sin(q2)*sin(q4))/5000 + 563/2000];
end

  