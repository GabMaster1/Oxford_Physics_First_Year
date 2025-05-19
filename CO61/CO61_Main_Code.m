%% CO61 Rocket Science ====================================================

% Cleaning up terminal
clear;
clc;
close all;

%% Functions ==============================================================

function [a] = acc_calc(p_r,p_m)
% Author: Gabriel Bristot, Date: 05/02/2025
% Defines the acceleration (2D vector) given a certain position of the moon
% and rocket the earth is considered stationary at the origin and scaled in Moon−radius.
% Input:
% * p_r: Position of the rocket (2D vector) with respect to the earth.
% * p_m: Position of the moon (2D vector) with respect to the earth.
% Output:
% * a: Acceleration vector (2D vector)

G = 9.63*10^(-7);  % Gravitational Constant
M_e = 83.3;        % Mass of the Earth
M_m = 1;           % Mass of the Moon

d_e = sqrt(sum(p_r.*p_r));                  % Distance from the rocket to the earth
d_m = sqrt(sum((p_r-p_m).*(p_r-p_m)));      % Distance from the rocket to the moon

a = - G*((M_e*p_r)/d_e^3) - G*((M_m*(p_r-p_m))/d_m^3); % Newton's gravitation formula

end

function [tout, pos] = simulate_rocket(init_pos , init_vel , moon_pos , t)
% Author: Gabriel Bristot, Date: 05/02/2025
% Simulate the rocket trajectory with the Earth and Moon influence. The coordinate
% used in this function is centred at Earth’s centre (i.e. Earth centre at (0,0) )
% and scaled in Moon−radius.
% The simulation finishes when it simulates for the whole t, or the rocket landed
% on the Moon.
% Input:
% ∗ init_pos: 2−elements vector (x, y) indicating the initial position of the rocket.
% ∗ init_vel: 2−elements vector (vx, vy) of the initial velocity of the rocket.
% ∗ moon_pos: a function that receives time, t, and return a 2−elements vector (x, y)
% indicating the Moon position relative to Earth.
% ∗ t: an N−elements vector of the time step where the position of the rocket will be
% returned.
%
% Output:
% ∗ tout: an M−elements vector of the time step where the position is described,
% if the rocket does not land on the Moon, M = N.
% ∗ pos: (M x 2) matrix indicating the positions of the rocket as function of time,
% with the first column is x and the second column is y.
      
    rocket_pos = init_pos; % Initial conditions for
    vel = init_vel;        % velocity and position

    collision = false;     % Collision set to false

    step = 1;              % Initiating counter to 1
    N = length(t);         % Maximum number of iterations

    tout = [];             % Defining vector containing the times of all positions calculated
    pos = [];              % Defining position matrix
    
    d_t = t(2)-t(1);       % Defining the delta t used

    while collision == false & step <= N % Looping over all N or until collision with the moon

        d_m = sqrt(sum((rocket_pos-moon_pos(t(step))).*(rocket_pos-moon_pos(t(step))))); % Distance from the rocket to the moon

        if d_m <= 1 % Collision check

            collision = true; 

        end

        % Start of Advanced Euler Method
        
        rocket_pos_intermidiate = rocket_pos + vel*d_t;

        vel_f = vel + (acc_calc(rocket_pos,moon_pos(t(step)))+acc_calc(rocket_pos_intermidiate,moon_pos(t(step)+t(1))))*(d_t/2);

        rocket_pos_f = rocket_pos + (vel_f+vel)*(d_t/2);

        rocket_pos = rocket_pos_f;
        vel = vel_f;
        
        % End of Andvanced Euler Method

        tout = [tout;t(step)];  % Updating the list of times
        pos = [pos;rocket_pos]; % Updating the matrix of positions
        
        step = step + 1; % Increasing the counter by 1

    end

end

%% Main code ==============================================================

% Creating the first figure of the stationary Moon case ===================
figure(1) 

% Setting the initial conditions of the first simulation
init_pos = [0, 3.7];
init_vel = 0.0066*[cosd(89.9), sind(89.9)];
t = linspace(0, 350000, 35000);

% Setting the position of the Moon over time
moon_pos = @(t) [0, 222];

% Running the simulation
[tout1, pos] = simulate_rocket(init_pos, init_vel, moon_pos, t);

% Plotting data
plot(pos(:,2), pos(:,1));

hold on

% Setting the initial conditions of the second simulation
init_pos = [3.7, 0];
init_vel = 0.0066*[cosd(51.5), sind(51.5)];
t = linspace(0, 350000, 35000);

% Setting the position of the Moon over time
moon_pos = @(t) [0, 222];

% Running the simulation
[tout2, pos] = simulate_rocket(init_pos, init_vel, moon_pos, t);

% Plotting data
plot(pos(:,2), pos(:,1));

% Making the graph
xlabel('Y position (Moon-radii)')
ylabel('X position (Moon-radii)')

title('The two different paths taken by the rocket, with the sationary Moon (\Deltat = 10s)')
legend('a) starting at [0, 3.7] Moon-radii, v_0 = 0.0066 Moon-radii/s, and \theta = 89.9^{\circ}', ...
       'b) starting at [3.7, 0] Moon-radii, v_0 = 0.0066 Moon-radii/s, and \theta = 51.5^{\circ}', ...
       'Location', 'southwest')

txt1 = ['\uparrow Earth'];
text(-2,-7,txt1)
txt2 = ['Moon \downarrow'];
text(199,10,txt2)

grid on
ylim([-125,125])
xlim([-10,240])


% Creating the second figure of the orbiting Moon case ====================
figure(2)

% Setting the initial conditions of the simulation
init_pos = [0, 3.7];
init_vel = 0.0066*[cosd(52.2), sind(52.2)];
w = 2.6615*10^(-6);
t = linspace(0, 350000, 35000);

% Setting the position of the Moon over time
moon_pos = @(t) 222*[cos(w*t),sin(w*t)];

% Running the simulation
[tout3, pos] = simulate_rocket(init_pos, init_vel, moon_pos, t);

% Plotting data
plot(pos(:,2), pos(:,1));

% Making the graph
xlabel('Y position (Moon-radii)')
ylabel('X position (Moon-radii)')

title('The path taken by the rocket, with the orbiting Moon (\Deltat = 10s)')
legend('a) starting at [0, 3.7] Moon-radii, v_0 = 0.0066 Moon-radii/s, and \theta = 52.2^{\circ}' ...
        ,'Location', 'southwest')

txt1 = '\uparrow Earth';
text(-2,-7,txt1)
txt2 = 'Moon \downarrow';
text(65,210,txt2)

grid on
ylim([-50,220])
xlim([-10,260])