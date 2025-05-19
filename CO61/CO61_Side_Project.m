%% CO61 Side project =================================================

% Cleaning up terminal
clear;
clc;
close all;

%% Functions =========================================================

function[acc] = acc_calc(p_r, p_m)
% Author: Gabriel Bristot , Date: 05/02/2025
% Defines the acceleration (2D vector) given a certain position of the moon
% and rocket the earth is considered stationary at the origin and scaled in Moon−radius.
% p_r: Position of the rocket (2D vector) with respect to the earth.
% p_m: Position of the moon (2D vector) with respect to the earth.
% a: Acceleration vector (2D vector)

G = 9.63*10^(-7);  % Gravitational Constant
M_e = 83.3;        % Mass of the Earth
M_m = 1;           % Mass of the Moon

d_e = sqrt(sum(p_r.*p_r));                  % Distance from the rocket to the earth
d_m = sqrt(sum((p_r-p_m).*(p_r-p_m)));      % Distance from the rocket to the moon

acc = - G*((M_e*p_r)/d_e^3) - G*((M_m*(p_r-p_m))/d_m^3); % Newton's gravitation formula

end

function[d_min] = simulate_rocket(init_pos ,init_vel ,moon_pos ,t)
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
% ∗ d_min: a scalar representing the minimum distance between the centre of
% the Moon and the rocket
      
    rocket_pos = init_pos; % Initial conditions for
    vel = init_vel;        % velocity and position

    collision = false;     % Collision set to false

    step = 1;              % Initiating counter to 1
    N = length(t);         % Maximum number of iterations
    
    d_t = t(2)-t(1);       % Defining the delta t used

    d_min = 1000;          % Defining the minimu idstance to the moon

    while collision == false & step <= N % Looping over all N or until collision with the moon
        
        d_e = sqrt(sum(rocket_pos.*rocket_pos));
        d_m = sqrt(sum((rocket_pos-moon_pos(t(step))).*(rocket_pos-moon_pos(t(step))))); % Distance from the rocket to the moon

        if d_m <= 1 || d_e <= 3.65  % Collision check

            collision = true; % Uptade collision

        end

        if d_m < d_min  % Check if new distance to the Moon is  smaller then minimum distance
                        
            d_min = d_m; % Update d_min

        end

        % Start of Advanced Euler Method
        
        rocket_pos_intermidiate = rocket_pos + vel*d_t;

        vel_f = vel + (acc_calc(rocket_pos,moon_pos(t(step)))+acc_calc(rocket_pos_intermidiate,moon_pos(t(step)+t(1))))*(d_t/2);

        rocket_pos_f = rocket_pos + (vel_f+vel)*(d_t/2);

        rocket_pos = rocket_pos_f;
        vel = vel_f;
        
        % End of Andvanced Euler Method
        
        % Increasing the counter by 1
        step = step + 1; 

    end

end

%% Main program ======================================================

% Number of points sampled in each axis
N = 100;

% Initial position
init_pos = [0, 3.7];

% Variables to be used in program
w = 2.6615*10^(-6);
t = linspace(0, 350000, 35000);

% Minimum and maximum values to be simulated
v_max = 0.005;
v_min = 0.010;

theta_max = 80;
theta_min = 90;

% Small change in variables every simulation
d_v = (v_max-v_min)/(N-1);
d_theta = (theta_max-theta_min)/(N-1);

% Position of the moon over time
moon_pos = @(t) 222*[0,1];

% defining the initial magnitude of the velocity
vel_mag = v_min;

% Pre defining a matrix M thta contains all the simulation data
Data = zeros(N*N,3);

% Keeping track of the number of rows
row = 1;

% Double for loop

% Cycles through v
for counter_1 = 1:N 
    
    % Starting theta at 0 in every inner loop
    theta = theta_min;

    % Cycles through theta
    for counter_2 = 1:N
        
        % Defining the velocity vector
        initial_vel = vel_mag*[cosd(theta), sind(theta)];
        
        % Simulating to find d_min
        d_min = simulate_rocket(init_pos, initial_vel, moon_pos, t);
        
        % Entering the data into data matrix 
        Data(row,1) = vel_mag;
        Data(row,2) = theta;
        Data(row,3) = d_min;
        
        % Updating counters
        theta = theta + d_theta;

        row = row + 1;

    end

    % Updating counter
    vel_mag = vel_mag + d_v;
    
    % Cosmetic detail of how far into the process the computation is
    fprintf('\r%3.0f%% completed', (counter_1/N) * 100);

end

%% Extract columns ===================================================

mag_v = 1000*Data(:,1);
theta = Data(:,2);
min_distance = Data(:,3);

% Get unique values for the grid
mag_v_vals = unique(mag_v);
theta_vals = unique(theta);

% Reshape the data into a grid
[X, Y] = meshgrid(mag_v_vals, theta_vals);
Z = griddata(mag_v, theta, min_distance, X, Y, 'cubic');

% Plot the surface
figure;
surf(X, Y, Z);
set(gca, 'ZScale', 'log')
xlabel('Speed (Moon-radii/s x 10^{-3})');
ylabel('Theta (degrees)');
zlabel('D_{min} (Moon-radii)');
title('Surface Plot of Minimum Distance to the Moon (D_{min})');
colorbar;
set(gca, 'ColorScale', 'log')
c = colorbar;
c.Label.String = 'D_{min} (Moon-radii)';
shading interp;

% Add contour lines on the surface
hold on;
contour3(X, Y, Z, 20, 'k');