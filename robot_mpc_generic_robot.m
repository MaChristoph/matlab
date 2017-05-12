

function MPC_Exercise2_Solution

% clear workspace, close open figures
clear all
close all
clc

% Due to some versions, there might be some (neglectable) warnings
warning off

% Set options 
tol_opt       = 1e-8;
options = optimset('Display','off',...
                'TolFun', tol_opt,...
                'MaxIter', 10000,...
                'Algorithm', 'sqp',...
                'FinDiffType', 'forward',...
                'RelLineSrchBnd', [],...
                'RelLineSrchBndDuration', 1,...
                'TolConSQP', 1e-6);

% system dimensions
n = 2; % state dimension
m = 3; % input dimension 


% Number of MPC iterations
mpciterations = 100;

% Horizon
N             = 5;

% Discretization steps
delta             = 0.1;

% initial states
tmeasure      = 0.0;
Theta1=pi/6;
Theta2=pi/6;
Theta3 = pi/6;
q = [Theta1; Theta2; Theta3];


% Initial guess for input
u0 = 0*ones(m,N);

% Weighting Matrices
Q = 1*eye(2,2);
R = 1*eye(3,3);
P = 11*eye(2,2);

% Set variables for output
t = [];
x = [];
u = [];

% Distance Constraints
min_dist_link1_ee = 0.2;
min_dist_link0_2 = 0.2;

min_dist_link0_ee = 1;

dist1 =[];
dist2 =[];
dist3 = [];
cost = [];

L1=1;
L2=1;
L3=1;
counter = 0;


% Print Header
fprintf('   k  |      u(k)        x(1)        x(2)     Time \n');
fprintf('---------------------------------------------------\n');

q = [Theta1; Theta2; Theta3];

% Desired Position
x_d = [-1;-1];

for ii = 1:mpciterations % maximal number of iterations
    Theta1 = q(1);
    Theta2 = q(2);
    Theta3 = q(3);

    %Forward kinematics
    xo = f(q);
    xmeasure = xo;
    qmeasure = q;

%     error = x_d - xo
    
%     if(sqrt(error' * error) < 0.05)
%         error
%         P
%        alpha = alpha*1.1; 
%        P = alpha*eye(2,2);
%     end
    %%Plot Robot
    %Calculate the location of the middle two joints
    pointl1 = [L1*cos(Theta1) ; L1*sin(Theta1)];
    pointl2 = pointl1 + [L2*cos(Theta1+Theta2);L2*sin(Theta1+Theta2)];


    
    if (mod(counter,1)==0) %plots every 10 iterations
        f1 = figure(1);
        axis([-3 3 -3 3])
        axis square
        line([0,pointl1(1)],[0,pointl1(2,1)])
        hold on
        line([pointl1(1),pointl2(1)],[pointl1(2,1),pointl2(2,1)])
        line([pointl2(1),xo(1)],[pointl2(2,1),xo(2,1)])
        plot(xo(1),xo(2),'o')
%         ang=0:0.01:2*pi; 
%         xp=min_dist_link0_2*cos(ang);
%         yp=min_dist_link0_2*sin(ang);
%         plot(pointl1(1)+xp,pointl1(2)+yp);
        
         ang=0:0.01:2*pi; 
        xp=min_dist_link0_ee*cos(ang);
        yp=min_dist_link0_ee*sin(ang);
        plot(xp,yp);
 
        pause(0.01)
    end

    %% MPC Stuff
    % Derive the open loop solution using the initial guess for the input/
    % old iteration for the input
    xOL = computeOpenloopSolution(@system, N, delta, tmeasure, xmeasure, qmeasure, u0);
    % not needed here, since we have no linear constraints on the state.

    % Set control and linear bounds for time instance
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    lb = [];
    ub = [];
    for k=1:N
        [Anew, bnew, Aeqnew, beqnew, lbnew, ubnew] = ...
               linearconstraints(tmeasure+k*delta,xOL(:,k),u0(:,k));
           
        A = blkdiag(A,Anew);
        b = [b, bnew];
        Aeq = blkdiag(Aeq,Aeqnew);
        beq = [beq, beqnew];
        lb = [lb, lbnew];
        ub = [ub, ubnew];
    end
 
    % Solve open loop optimization problem
    [u_OL, V, exitflag, output] = fmincon(@(u) costfunction(@runningcosts, ...
        @terminalcosts, @system, N, delta, tmeasure, xmeasure, qmeasure, u,...
        xOL, Q, R, P, x_d, min_dist_link1_ee, min_dist_link0_2,min_dist_link0_ee), ...
        u0, A, b, Aeq, beq, lb, ub, ...
        @(u) nonlinearconstraints(@constraints, @terminalconstraints, ...
        @system, N, delta, tmeasure, xmeasure, qmeasure, u, P, x_d), options);
 
    % Prepare restart (take the optimal open loop solution as initial
    % guess and add a last piece (LQR))
%     u0 = [u_OL(:,2:size(u_OL,2)) -[2.118 2.118]*x_OLapllied(end,:)'];
    u0 = u_OL;

    %% Apply first optimal input
    q = q + u_OL(:,1) * delta;


    %% Stuff for performance analysis
        % Store closed loop data
    t = [ t; tmeasure ];
    x = [ x; xmeasure ];
    u = [ u; u_OL(:,1) ];
    
    tmeasure = tmeasure + delta;
    
    % Plot closed-loop trajectories
    % Apply control to process for the first iteration
    q_OLapplied = computeOpenloopSolution(@system, N, delta, tmeasure, xmeasure, qmeasure, u_OL);
 
    x_OLapllied = [];
    for i=1:size(q_OLapplied)
        x_OLapllied = [x_OLapllied f(q_OLapplied(:,i))];
    end
%     x_OLapllied = f(q_OLapplied)
    f2 = figure(2);
    title('Open-loop trajectories');
    plot(x_OLapllied(1,:),x_OLapllied(2,:),'g')
    hold all;
    plot(x_OLapllied(1,1),x_OLapllied(2,1),'o')
    xlabel('x(1)')
    ylabel('x(2)')  
    hold all;
    ang=0:0.01:2*pi; 
    xp=min_dist_link0_ee*cos(ang);
    yp=min_dist_link0_ee*sin(ang);
    plot(xp,yp);
    axis([-3 3 -3 3])
    axis square
    
    dd1 = (f(q) - fk_link1(q))' * (f(q) - fk_link1(q));
    dd2 = fk_link1(q)' * fk_link1(q);
    dd3 = f(q)' *  f(q);
    dist1 = [dist1; sqrt(dd1)];
    dist2 = [dist2; sqrt(dd2)]; 
    dist3 = [dist3; sqrt(dd3)]; 

     % Define Barrier Functions
    barrier_link0_2 =exp((min_dist_link0_2 - sqrt(dd2))/0.05);
    barrier_link1_ee = exp((min_dist_link1_ee - sqrt(dd1))/0.05);
    barrier_link0_ee = exp((min_dist_link0_ee - sqrt(dd3))/0.05);
    
    % Provide the running cost
    cl_cost = 0.5*(f(q)-x_d)'*Q*(f(q)-x_d) + u_OL(:,1)'*R*u_OL(:,1)...
        +barrier_link0_2 + barrier_link1_ee + barrier_link0_ee;
    
    cost = [cost; cl_cost];

  
  counter = counter +1;
end
figure

subplot(3,1,1)
plot(t,dist1, t,min_dist_link1_ee * ones(1,mpciterations),'r--');
legend('Distance','Threshold');
title('Distance from Link 1 to End-effector')
axis([-0.01 mpciterations*delta -0.01 5])
xlabel('t [s]')
ylabel('distance [m]')

subplot(3,1,2)
plot(t,dist2, t,min_dist_link0_2 * ones(1,mpciterations),'r--');
legend('Distance','Threshold');
title('Distance from Link 0 to 2')
axis([-0.01 mpciterations*delta -0.01 5])
xlabel('t [s]')
ylabel('distance [m]')

subplot(3,1,3)
plot(t,dist3, t,min_dist_link0_ee * ones(1,mpciterations),'r--');
legend('Distance','Threshold');
title('Distance from Link 0 to End-effector')
axis([-0.01 mpciterations*delta -0.01 5])
xlabel('t [s]')
ylabel('distance [m]')

figure
plot(t,cost);
title('Cost')


end

function qdot = system(t, u, delta)
    qdot = u;
end

function cost = costfunction(runningcosts, terminalcosts, system, ...
                    N, delta, t0, x0, q0, u, xOL, Q, R, P, x_d, min_dist_link1_ee, min_dist_link0_2, min_dist_link0_ee)
    % Formulate the cost function due to which the open loop optimization
    % will be minimized
    cost = 0;
  
    % First, derive the open loop state values for the open loop input of
    % the former iteration
%     x = zeros(N+1, length(x0));
    q = computeOpenloopSolution(system, N, delta, t0, x0, q0, u);
                            
    %Second, build the overall cost by means of the running and the
    %terminal cost
%     for k=1:N
%         cost = cost+runningcosts(t0+k*delta, x(:,k), u(:,k), Q, R, x_d);
%     end

    %% With penalty functions
    for k=1:N

        cost = cost+runningcosts(t0+k*delta, q(:,k), u(:,k), Q, R, x_d, min_dist_link1_ee, min_dist_link0_2, min_dist_link0_ee);
    end

    cost = cost+terminalcosts(t0+(N+1)*delta, q(:,N+1), P, x_d, min_dist_link1_ee, min_dist_link0_2, min_dist_link0_ee);
end

function cost = runningcosts(t, q, u, Q, R, x_d, min_dist_link1_ee, min_dist_link0_2, min_dist_link0_ee)
   
    
    % Get Distances
    d1_ee = (f(q) - fk_link1(q))' * (f(q) - fk_link1(q) );
    d0_2 = (fk_link2(q))' * (fk_link2(q));
    d0_ee = (f(q))' * (f(q));
    
    % Define Barrier Functions
    barrier_link0_2 =exp((min_dist_link0_2 - sqrt(d0_2))/0.05);
    barrier_link1_ee = exp((min_dist_link1_ee - sqrt(d1_ee))/0.05);
    barrier_link0_ee = exp((min_dist_link0_ee - sqrt(d0_ee))/0.05);
    
    % Get EE Position
    x = f(q);
    
    % Provide the running cost
    cost = 0.5*(x-x_d)'*Q*(x-x_d) + u'*R*u + barrier_link0_2 + barrier_link1_ee + barrier_link0_ee;
%     cost = 0.5*(x-x_d)'*Q*(x-x_d) + u'*R*u;
end

function [A, b, Aeq, beq, lb, ub] = linearconstraints(t, x, u)
    % Introduce the linear constraints
    A   = [];
    b   = [];
    Aeq = [];
    beq = [];
    lb  =  [];
    ub  =  [];
    
end

function [c,ceq] = nonlinearconstraints(constraints, ...
    terminalconstraints, system, N, delta, t0, x0, q0, u, P, x_d)
    
    % Introduce the nonlinear constraints also for the terminal state
%     x = zeros(N+1, length(x0));
    x = computeOpenloopSolution(system, N, delta, t0, x0, q0, u);
    c = [];
    ceq = [];
    for k=1:N
        [cnew, ceqnew] = constraints(t0+k*delta,x(:,k),u(:,k));
        c = [c cnew];
        ceq = [ceq ceqnew];
    end
    [cnew, ceqnew] = terminalconstraints(t0+(N+1)*delta,x(:,N+1), P, x_d);
    c = [c cnew];
    ceq = [ceq ceqnew];
    
end

function cost = terminalcosts(t, q, P, x_d, min_dist_link1_ee, min_dist_link0_2, min_dist_link0_ee)
 
    % Get Distances
    d1_ee = (f(q) - fk_link1(q))' * (f(q) - fk_link1(q) );
    d0_2 = (fk_link2(q))' * (fk_link2(q));
    d0_ee = (f(q))' * (f(q));
    
    % Define Barrier Functions
    barrier_link0_2 = exp((min_dist_link0_2 - sqrt(d0_2))/0.05);
    barrier_link1_ee = exp((min_dist_link1_ee - sqrt(d1_ee))/0.05);
    barrier_link0_ee = exp((min_dist_link0_ee - sqrt(d0_ee))/0.05);

    % Get EE Position
    x = f(q);
    
    % Intro the terminal cost
    cost = (x-x_d)'* P *(x-x_d) + barrier_link0_2 + barrier_link1_ee + barrier_link0_ee;
% cost = (x-x_d)'* P *(x-x_d);
end

function [c,ceq] = constraints(t, x, u)


%     d1 = (f(q) - fk_link1(q))' * (f(q) - fk_link1(q) );
%     c = -sqrt(d1) + 0.6;
    c = [];
    ceq = [];

end

function [c,ceq] = terminalconstraints(t, x, P, x_d)

%     d1 = (f(q) - fk_link1(q))' * (f(q) - fk_link1(q) );
%     x = x(1:2);
    % Intro the terminal constraint
%     c   = (x-x_d)'* P *(x-x_d);
%     c = -sqrt(d1) + 1.6;
    c = [];
    ceq = [];

end

function q = computeOpenloopSolution(system, N, delta, t0, x0, q0, u)
                                
    % Apply open loop input u to the system
    q(:,1) = q0;
  
    
    for k=1:N
        q(:,k+1) = q(:,k) + system(0, u(:,k), 0) * delta;
    end 
    
end

function x_ee = f(q)
    L1=1;
    L2=1;
    L3=1;
    x_ee = [L1*cos(q(1))+L2*cos(q(1)+q(2))+L3*cos(q(1)+q(2)+q(3));...
            L2*sin(q(1))+L2*sin(q(1)+q(2))+L3*sin(q(1)+q(2)+q(3))];
end

function x = fk_link2(q)
    L1=1;
    L2=1;
    dummy = q(3);
    x = [L1*cos(q(1))+L2*cos(q(1)+q(2));...
            L2*sin(q(1))+L2*sin(q(1)+q(2))];

end


function x = fk_link1(q)
    L1=1;
    L2=1;
    dummy = q(2);
    dummy = q(3);
    x = [L1*cos(q(1));...
         L2*sin(q(1))];
end