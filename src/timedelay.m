clear;clc;close all

throttle_previous = sym('tp');
throttle_optimized = sym('to');
steering_previous = sym('sp');
steering_optimized = sym('so');
v_previous = sym('vp');
v_current = sym('vc');
psi_previous = sym('pp');
psi_current = sym('pc');
% syms throttle_optimized throttle_optimized
% syms steering_previous steering_optimized
% syms v_previous v_current
% syms psi_previous psi_current
syms time_delay time_until_next_control
syms Lf



v_intermediate = v_previous + throttle_previous * time_delay;
v_current_ = v_intermediate + throttle_optimized * time_until_next_control;

psi_intermediate = psi_previous + v_previous * steering_previous / Lf * time_delay;
psi_current_ = psi_intermediate + v_intermediate * steering_optimized / Lf * time_until_next_control;

% Now solve this system of equations for the times 
solution = solve( v_current - v_current_, psi_current - psi_current_, time_delay, time_until_next_control );
td = simplify( solution.time_delay )