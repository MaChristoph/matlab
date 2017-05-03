%% Define Symbolic Variables
syms q1 q2 q3 q4 q5 q6 q7 q8 q9 q10
syms T1(q1)

syms x y

syms T

syms theta d a alpha

%% Symbolic function which represents a homogeneous transform for Denavit-Hartenberg
T1(theta, d, a, alpha) = [  cos(theta)   , -sin(theta) * cos(alpha)   , sin(theta) * sin(alpha) , a * cos(theta);
                            sin(theta)   , cos(theta) * cos(alpha)    , -cos(theta) * sin(alpha), a * sin(theta);
                            0            ,   sin(alpha)               , cos(alpha)              , d;
                            0            , 0                          , 0                       , 1];

T_base(x, y, theta) = [ cos(theta)   , -sin(theta) , 0 , x;
                        sin(theta)   ,  cos(theta) , 0 , y;
                        0            ,  0          , 1 , 0;
                        0            ,  0          , 0 , 1];
%% Define Forward-Kinematics with DH-parameters
% T(q1,q2,q3,q4,q5,q6,q7) = T1(q1, 0.2815, 0, pi/2) * ...
%                           T1(q2, 0, 0, -pi/2) * ...
%                           T1(q3, 0.4536, 0, pi/2) * ...
%                           T1(q4, 0, 0, -pi/2) * ...
%                           T1(q5, 0.2986, 0, pi/2) * ...
%                           T1(q6, 0, 0, -pi/2) * ...
%                           T1(q7, 0, 0, 0);

% T(q1,q2,q3,q4,q5,q6,q7) = T1(q1, 0.2815, 0, pi/2) * ...
%                           T1(q2, 0, 0, -pi/2);


% T(q1,q2,q3) = T1(q1,0,1,0) * ...
%               T1(q2,0,1,0) * ...
%               T1(q3,0,1,0);

%               (theta, d, a, alpha)

% COB                                    T1(q4, 0.15314, 0, 0) * ...          % arm_left_base_link -> arm_left_1_link
% T(q1,q2,q3,q4,q5,q6,q7,q8,q9,q10) = T_base(q1, q2, q3) * ...            %base_link
%                                     T1(pi +q4 , 1.0065, 0, pi/2) * ...   % base_link          -> arm_left_base_link
%                                     T1(q5, 0.0915+0.15314 , 0, pi/2) * ...      % arm_left_1_link    -> arm_left_2_link
%                                     T1(q6, 0, 0, -pi/2) * ...           % arm_left_2_link    -> arm_left_3_link
%                                     T1(q7, 0.4536, 0, pi/2) * ...       % arm_left_3_link    -> arm_left_4_link
%                                     T1(q8, 0, 0, -pi/2) * ...           % arm_left_4_link    -> arm_left_5_link
%                                     T1(q9, 0.2986, 0, pi/2) * ...       % arm_left_5_link     -> arm_left_6_link
%                                     T1(q10, 0, 0, -pi/2)                                


% Arm
% T(q1,q2,q3,q4,q5,q6,q7,q8,q9,q10) =   T_base(q1, q2, q3) * ...
%                                       T1(q4, 0.2815, 0, pi/2) * ...
%                                       T1(q5, 0, 0, -pi/2) * ...
%                                       T1(q6, 0.4536, 0, pi/2) * ...
%                                       T1(q7, 0, 0, -pi/2) * ...
%                                       T1(q8, 0.2986, 0, pi/2) * ...
%                                       T1(q9, 0, 0, -pi/2) * ...
%                                       T1(q10, 0, 0, 0);



T(q1,q2,q3,q4,q5,q6,q7) = T1(q1, 0.2815, 0, pi/2)


%% Get Endeffector Position
pos = T * [0, 0, 0, 1]';

%% Simplify analytic expression
T_simplified = simplify(T,'Steps',10);

T_simplified = T_simplified(q1,q2,q3,q4,q5,q6,q7);
% T_simplified = T_simplified(q1,q2,q3);

% 0.191, 0.760, 0.835
% hehe = T_simplified(0,0,0,0,0,0,0,0,0,0);

% hehe = T_simplified(0,0,0, -0.9655394142680995, -1.2848527807718577, -0.4245127097461321, 1.7204204604569675, -0.7581014523483809, -0.6142954362219815, 1.7479758113542214);
% hehe = T_simplified(0,0,0, 0.0295, -0.0788, 0.591, 1.426, 0.044, -1.359, -0.621)
%  hehe = T_simplified(0,0,0 ,0.4197, 0.755, -0.696, -1.945, -0.462, 1.349, 0.237)

%     pos = double(hehe) * [0,0,0,1]';
%     pos = pos(1:3)
%  rotm = hehe(1:3, 1:3);
%  rotm = double(rotm);
 
 % quat = [w, x,y,z]
%  quat = rotm2quat(rotm)
%  T_simplified = T_simplified(q1,q2,q3,q4,q5,q6,q7,q8,q9,q10);

%% Transform to Position Vector and Rotation Matrix
tmp = T_simplified * [0, 0, 0, 1]';
tmp = simplify(tmp,'Steps',10);
p = tmp(1:3)

R = T_simplified(1:3,1:3);
R = simplify(R,'Steps',10)

     
