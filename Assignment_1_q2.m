% Name: Evan Chiu
% Date: 20240319-20240331
% Description: Assignment 1 - Question 2
% Input(s): angle, frequency of reps, mass of weight
% Output(s):
%
clc
clear
close all
%% Question 2 [28 Marks]
% Part C - Coding
%
% Constants
g = 9.81; % gravitational acceleration [m/s^2]
m_arm= 3.5; % [kg]
a = 15e-2; % [m]
b = 45e-2; % [m]
mu_u = 20e-2; % upper arm muscle attachment distance [m]
mu_d = 3e-2; % forearm muscle attachment distance [m]
t = linspace(0,30,2000); % [s]
%
% Input code for angle in degrees
angle = input('Give an angle in degrees between 0 and 70 [degrees theta].\n'); % Requesting an angle input in degrees
    % Error Code
    while isempty(angle) || angle < 0 || angle > 70
        disp('Invalid angle, try again')
        angle = input('Give an angle in degrees.\n');
    end
    %
% Input code for reps per minute
f = input('Input the number of reps per minute.\n'); % Requesting frequency of reps in a minute
    % Error Code
    while isempty(f) || f < 0
        disp('Invalid number of reps, try again')
        f = input('Give an number of reps.\n');
    end
    %
% Input code for mass of weight
m_weight = input('Give the mass of the carried weight [kg].\n'); % Requesting frequency of reps in a minute
    % Error Code
    while isempty(m_weight) || m_weight < 0
        disp('Invalid mass, try again')
        m_weight = input('Give a mass.\n');
    end
    %
% Expressions
theta_0 = deg2rad(angle); % inputted angle in [rad]
f_seconds = f/60; % frequency per second rather than per minute
phi=2*pi*f_seconds; % angular frequency of the motion
I_o = 0.08875 + 0.2065*m_weight; % mass moment of inertia
r_com = (0.525+0.45*m_weight)/(3.5+m_weight); % location of the centre of mass of the system
theta = theta_0*cos(phi*t); % angular displacement in [rad]
theta_deg = rad2deg(theta); % angular displacement in [deg]
omega = phi*theta_0.*-sin(phi*t); % angular velocity [rad/s]
alpha = phi^2*theta_0.*-cos(phi*t); % angular acceleration [rad/s^2]
% Muscle force
F_muscle1 = sqrt((mu_u^2)+(mu_d^2)-2*mu_u*mu_d*sin(theta));
F_muscle2 = (mu_u*mu_d*cos(theta));
F_muscle3 = I_o*alpha + (a*m_arm + b*m_weight)*g*cos(theta);
F_muscle = (F_muscle1./F_muscle2).*F_muscle3;
% Parallel reaction force
R_par1 = F_muscle .* ((mu_d-mu_u*sin(theta))./F_muscle1);
R_par2 = (m_arm+m_weight)*(g*sin(theta)-(omega.^2)*r_com);
R_par = R_par1 + R_par2;
% Perpendicular reaction force
R_perp1 = (m_arm+m_weight)*(g*cos(theta)+alpha*r_com);
R_perp2 = F_muscle.*((mu_u*cos(theta))./F_muscle1);
R_perp = R_perp1 - R_perp2;
%
% Plotting
% Coordinates for arm, weight and muscle
forearm_x=b*cos(theta); % defines x-position of forearm
forearm_y=b*sin(theta); % defines y-position of forearm
weight_x=b*cos(theta); % defines x-position of weight
weight_y=b*sin(theta); % defines y-position of weight
muscle_x=mu_d*cos(theta); % defines x-position of muscle
muscle_y=mu_d*sin(theta); % defines y-position of muscle
% Subplotting scheme
subplot(4,2,[1,3,5,7]) % Large subplot on the left for bicep curl animation
% Setting up 2 animatedline objects, corresponding to the forearm and the bicep muscle.
upperarm=animatedline('Color','k','LineWidth',1.5);
forearm=animatedline('Color','K','LineWidth',1.5);
muscle=animatedline('Color','r','LineWidth',2.0);
weight=animatedline('Marker','o','MarkerSize',20,'MarkerFaceColor','w','MarkerEdgeColor','k','LineWidth',1.5);
title('Animation of bicep curl')
xlabel('')
ylabel('')
axis equal % sets axis to equal unit sizes so animation does not look weird
axis([-0.2 b+0.2 -b-0.2 b+0.2])
%
subplot(4,2,2) % first righthand subplot describing angular displacement
thetaline=animatedline('Color','b','LineWidth',1.5);
title('Graph of \theta [deg] vs t [s]')
xlabel('t [s]')
ylabel('\theta [deg]')
axis([0 30 -angle angle])
%
subplot(4,2,4) % second righthand subplot describing parallel reaction force vs time
parline=animatedline('Color','r','LineWidth',1.5);
title('Graph of R_{par} [N] vs t [s]')
xlabel('t [s]')
ylabel('R_{par} [N]')
axis([0 30 min(R_par) max(R_par)])
%
subplot(4,2,6) % third righthand subplot describing perpendicular reaction force vs time
perpline=animatedline('Color','g','LineWidth',1.5);
title('Graph of R_{perp} [N] vs t [s]')
xlabel('t [s]')
ylabel('R_{perp} [N]')
axis([0 30 min(R_perp) max(R_perp)])
%
subplot(4,2,8) % fourthe righthand subplot describing perpendicular reaction force vs time
muscleline=animatedline('Color','m','LineWidth',1.5);
title('Graph of F_{muscle} [N] vs t [s]')
xlabel('t [s]')
ylabel('F_{muscle} [N]')
axis([0 30 min(F_muscle) max(F_muscle)])
%
% Animations
n=length(t); % determines number of points in t
for i=1:n
    subplot(4,2,[1,3,5,7])
    clearpoints(forearm) % clears previous points in forearm
    clearpoints(weight) % clears previous points in weight
    clearpoints(muscle) % clears previous points in muscle
    addpoints(upperarm,[0,0],[0,mu_u+0.2]);
    addpoints(forearm,[0,forearm_x(i)],[0,forearm_y(i)]);
    addpoints(muscle,[0,muscle_x(i)],[mu_u,muscle_y(i)]);
    addpoints(weight,weight_x(i),weight_y(i));
    drawnow
    % animates plot of angular displacement
    subplot(4,2,2)
    addpoints(thetaline,t(i),theta_deg(i))
    %
    % animates plot of angular velocity
    subplot(4,2,4)
    addpoints(parline,t(i),R_par(i))
    %
    % animates plot of angular acceleration
    subplot(4,2,6)
    addpoints(perpline,t(i),R_perp(i))
    %
    % animates plot of tension
    subplot(4,2,8)
    addpoints(muscleline,t(i),F_muscle(i))
    drawnow
end