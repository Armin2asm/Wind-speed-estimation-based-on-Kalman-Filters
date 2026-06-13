%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Effective Wind Speed Estimation Using Kalman-Based Filters
%
% MATLAB implementation for:
%   1- Linear Kalman Filter (KF)
%   2- Extended Kalman Filter (EKF)
%   3- H-infinity Kalman Filter (HKF)
%
% Wind turbine model : NREL 5 MW
% Wind profiles      : 10.5 m/s
%
% Author : Armin Sarkoobi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
close all
clear all

beta=0;%=============================================================
N=3;            %number of blades
mT=436865;         % mass of the tower and nacelle
mB=4435;    %mass of each blade
rB=21.975;

Jr=38759227;             %Inertia of the rotor
Jg=5025347;              %inertia of generator
Tg=0.02;               %Applied generator torque

M=[mT+(N*mB)   N*mB   0   0
    N*mB       N*mB   0   0
    0           0     Jr  0
    0           0     0   Jg];
%===============================================================

dT=70000  ;       %Damping of tower
dS=6215000;        %Damping of transmission
dB=20000;         %Damping of blade
ng=97;              %gearbox ratio

D=[ dT    0      0      0 
    0    N*dB    0      0 
    0     0      dS    -dS/ng 
    0      0     -dS/ng    dS/(ng)^2 ];

%================================================================
KT=1962000;      %stiffness of each tower
KB=40000;      %stiffness of each blade
KS=867637000;      %stiffness of each transmission
% K is Stiffness matrix for model with torsion angle
K= [KT   0     0       0 
    0   N*KB   0       0
    0    0     KS    -KS/ng
    0    0    -KS/ng     KS/(ng)^2];
%Stiffness matrix

Kmad=[KT      0      0 
      0      N*KB    0
      0       0      KS
      0       0     -KS];
%Stiffness matrix for model with torsion angle

%=================================================================

Z33=[0 0 0 ;0 0 0 ;0 0 0 ];
Z31=[0; 0; 0];
Z13=[0 0 0];
Z14=[0 0 0 0];
Z41=[0;0;0;0];
Z71=[0;0;0;0;0;0;0];
Z2=[0 0 0;0 0 0 ;0 0 0 ;0 0 0 ];
L34=[1 0 0 0;0 1 0 0 ;0 0 1 -1];
O1=[1 0 0 0;0 1 0 0 ;0 0 1 0;0 0 0 1];
W=[Z31;-1/Jg];

T=0.1;           %Delay time constant for pitch dynamics
Tg=0.02;        %Delay time constant for torque dynamics

A=[Z33                  L34                Z31        Z31
    -(M^(-1))*Kmad   -(M^(-1))*D           Z41         W 
    Z13          Z14              -1/T        0
    Z13          Z14               0       -1/Tg];

B=[Z71     Z71
   1/T     0
    0    1/Tg];

Q=[1   0   0 
    1  0   0 
    0  1   0 
    0  0  -1];
C=eye(9);



ro=1.225;         %Air density
R=63;          %Rotor radius
v=8;      %Wind speed
landa=2;

c1=0.005;
c2=1.53;
c3=0.5;
c4=0.18;
c5=121;
c6=27.9;
c7=198;
c8=2.36;
c9=5.74;
c10=11.35;
c11=16.1;
c12=201;

landai=1/(landa+0.08*beta)-(0.035/c11+c12*beta^3);

a1=0.006;
a2=0.095;
a3=-4.15;
a4=2.75;
a5=0.001;
a6=7.8;
a7=-0.00016;
a8=-8.88;

CTmad=a1+[a2*(landa-a3*beta)*exp(-a4*beta)]+[a5*(landa^2)*exp(-a6*beta)]+[a7*(landa^3)*exp(-a8*beta)];
CT=CTmad*((1+sign(CTmad))/2);



CQmad=c1*(1+c2*((beta+c3)^1/2))+(c4/landa)*(c5*landai-c6*beta-c7*beta)*exp(-c10*landai);
CQ=CQmad*((1+sign(CQmad))/2);

FT=((ro *pi*R^2)/2)*CT*v^2;   
Ta=((ro*pi*R^3)/2)*CQ*v^2;  
Bd=0.1745;                %Demanded pitch angle


u=[Bd Tg ]';
q=eye(9);

r=0.1;

wk1=[Z41 ; (1/(N*mB))*FT; (1/(Jr))*Ta; Z31 ];
    vk=[randn];
wk2=sqrt(q)*wk1*[randn]';





%%  input data
wdata=xlsread('winddata10.5.xlsx')
index=1;
%%
    X=wdata(:,2); % time wind 

% kalman Filter

u=[Bd Tg ]';
q=eye(9);

r=0.1;

wk1=[Z41 ; (1/(N*mB))*FT; (1/(Jr))*Ta; Z31 ];
    vk=[randn];
wk2=sqrt(q)*wk1*[randn]';
%=============================================================
%=============================================================

x = linspace(0, 5, 666);  
y =X; % Y-axis

% Number of data points
num_points = length(x);  

% State matrix
dt = 1/20;   
A = [1 dt; 0 1]; % State transition model
  
H = [1 0]; % Observation matrix
  

% Initial state
x_est = [0; 0]; % position and velocity  
P = eye(2);  

% Process and measurement noise covariance matrices
Q = [0.01 0; 0 0.01];   
R = 5;   

% Arrays for storing the results
filtered_estimates = zeros(2, num_points);  

%% ================= Linear Kalman Filter =====================
for i = 1:num_points  
    % prediction
    x_pred = A * x_est;   
    P_pred = A * P * A' + Q;   

    % update 
    z = y(i); %   
    y_tilde = z - H * x_pred; 
    S = H * P_pred * H' + R;   
    K = P_pred * H' / S; % Kalman gain  
    x_est = x_pred + K * y_tilde;  
    P = (eye(2) - K * H) * P_pred; 

    % Save result  
    filtered_estimates(:, i) = x_est;  
end 

%output(index:index+8,1)=y(:,1)


% Estimation of position extraction 
true_state=X;



% plot  Kalman filter
%=============================================

time=linspace(0,5,666);

figure;  
plot(time,  filtered_estimates(1, :), 'DisplayName', 'Kalman filter','color',[0.8 0.7 0],'lineWidth',1.2);  

hold on;  

estimated_position1=filtered_estimates(1, :);

%
%
%==============================================================
%==============================================================
%% ================= Extended Kalman Filter ==================

wdata=xlsread('winddata10.5.xlsx')
    X=wdata(:,2);
% Number of points  
N = 666;   

time = linspace(0, 5, N);   

% Define true values (randomly generated for demonstration)  
true_position = X;   
measurements = X;   

% Define measurements (with added noise)  
noise = 0.5 * randn(1, N);   
measurements = true_position + noise;   

% Initialize Extended Kalman Filter  
x_hat = [true_position(1); 0];   
P = eye(2); 

% Model parameters  
Q = [0.01, 0; 0, 0.01];   
R = 0.9;  

% Estimation history  
x_hat_history = zeros(2, N);   

for k = 1:N  
% For prediction step  
    if k > 1  
        dt = time(k) - time(k - 1);   
    else  
        dt = 0;  
    end  
    
% State transition matrix  
    F = [1, dt; 0, 1];  
    x_hat_pred = F * x_hat;   
    P_pred = F * P * F' + Q;   

% Measurement function  
    H = [1, 0];  
    z_pred = H * x_hat_pred;   

% Compute nonlinear symmetry  
    K = P_pred * H' / (H * P_pred * H' + R);   

% Update step  
    x_hat = x_hat_pred + K * (measurements(k) - z_pred);   
    P = (eye(size(K, 1)) - K * H) * P_pred;   

% Save estimate  
    X_hat_history(:, k) = x_hat;   
end  

% plot Extended Kalman filter
%=============================================
plot(time, X_hat_history(1, :), 'g', 'DisplayName', 'Extended Kalman filter','LineWidth',1.2);  
hold on
legend;  
xlim([0 5]);   
ylim([8 14]);  
hold on
% Compute difference between true state and estimate  true_state=X;
difference2 =  X_hat_history(1,:)'- true_state ;  

% Position estimate extraction  


estimated_position2 = X_hat_history(1, :);  



%==================================================
%===================================================
% H-infinity kalman filter
wdata=xlsread('winddata10.5.xlsx')
    X=wdata(:,2);
A6=X;
N = 666;   

% Create time array  time = linspace(0, 5, N);   

% Create arrays for true states and measurements  
true_state = X 
measurements = X 

% Initialize H-infinity Kalman Filter  
x_hat = true_state(1);   
P = 1;   

% Model matrices  
A = 1;  
H = 1;
Q = 0.01; 
R = 0.51;   

% Estimation history  
x_hat_history = zeros(1, N);   

for k = 1:N  
% prediction step 
    x_hat_pred = A * x_hat;   
    P_pred = A * P * A' + Q;   

% Update step  
    K = P_pred * H' / (H * P_pred * H' + R);   
    x_hat = x_hat_pred + K * (measurements(k) - H * x_hat_pred);   
    P = (1 - K * H) * P_pred;   

% Save estimate  
    x_hat_history(k) = x_hat;   
    estimated_position3= x_hat_history;
end  

% plot H- infinity Kalman filter
%=============================================
plot(time, x_hat_history, 'b', 'DisplayName', 'H-infinity Kalman filter','LineWidth',1.2)

% plot x-true
%=============================================
plot(time, true_state , 'r', 'DisplayName', 'Actual speed','LineWidth',1)

ylim([8 14])
   xlim([0 5])
 %% ================= Performance Evaluation ===================

 % Compute difference between true state and estimate    
difference3 =  x_hat_history'- true_state ; 

%  MSE  
mse1 = mean((true_position -  estimated_position1').^2);  
mse2 = mean((true_position - estimated_position2').^2);  
mse3 = mean((true_position - estimated_position3').^2);  


%  MSE  
disp(['Mean Squared Error (MSE kf): ', num2str(mse1)]) 
disp(['Mean Squared Error (MSE EKF): ', num2str(mse2)]) 
disp(['Mean Squared Error (MSE hKF): ', num2str(mse3)]) 
%% ================= Plot Results =============================

% Plot difference  
difference1 =  filtered_estimates(1, :)'- true_state ;
 ylabel('wind speed (m/s)')
 xlabel( 'Time(sec)')
hold off

figure;   
plot(time, difference1,'color',[0.8 0.7 0],'DisplayName', 'Kalman filter','lineWidth',1.2);  
hold on
plot(time, difference2,'g', 'DisplayName', 'Extended Kalman filter','LineWidth',1.2);  
hold on
plot(time, difference3,'b', 'DisplayName', 'H-infinity Kalman filter','LineWidth',1.2)
hold on
 
ylim([-2 2])
   xlim([0 5])



% fprintf('Mean Squared Error (MSE): %.4f\n', mse);  
 ylabel(' Estimation error (m/s)')
 xlabel( 'Time(sec)')
 grid on
hold off
   

