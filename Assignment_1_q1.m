% Name: Evan Chiu
% Date: 20240319-20240331
% Description: Assignment 1 - Question 1
% Input(s): theta
% Output(s):R_parallel, R_perp, m0
%
clc
clear
close all
%% Question 1 [22 Marks]
% Part C - Coding
%
% Constants
g = 9.81;
m_arm= 3.5; % [kg]
m_weight=20; % [kg]
a = 15e-2; % [m]
b = 45e-2; % [m]
%
% Input Code
theta = input('Give an non-zero angle in degrees.\n'); % Requesting an angle input in degrees
% Error Code
while isempty(theta) | theta == 0; 
    disp('Invalid angle, try again')
    theta = input('Give an angle in degrees.\n');
end
%
fprintf('Your angle is %.2f degrees.\n',theta)
%
% Reporting the three quantities
R_parallel=abs(-(g*m_arm)*sind(theta) - (g*m_weight)*sind(theta)); % Expression for parallel force
fprintf('The reaction force at Point O parallel to the arm segment is %.2f N.\n',R_parallel)
R_perp=abs((g*m_arm)*cosd(theta) + (g*m_weight)*cosd(theta)); % Expression for perpendicular force
fprintf('The reaction force at Point O perpendicular to the arm segment is %.2f N.\n',R_perp)
m0=abs(-a*m_arm*g*cosd(theta)-b*m_weight*g*cosd(theta)); % Expression for moment
fprintf('The moment M0 at Point O is %.2f N.\n',m0)
%