clc; clear; close all;
% Simulates the four-link snake robot with one of the eight methods and
% plays a 3D animation of its motion. Set saveVideo = 1 to also save it as
% an MP4 file.

global m1 l1 m2 l2 m3 l3 m4 l4 k1 k2 k3
m1 = 3.5; m2 = 3; m3 = 2.5; m4 = 2;
l1 = 2; l2 = 1.5; l3 = 1; l4 = 0.5;
k1 = 3; k2 = 4; k3 = 2;

%_________________________Settings_________________________________________
% Method: 'IntMul', 'Augmented', 'Elimination', 'Embedding',
% 'ModifiedLagrange', 'Greenwood', 'GibbsAppell' or 'Kane'
method = 'ModifiedLagrange';
T = 30;                  % simulation time (s)
fps = 30;                % frames per second of the animation
speed = 2;               % playback speed (2 = twice as fast as real time)
saveVideo = 0;           % 1 to save the animation as an MP4 file

% Initial conditions (same as RunMe.m): the robot starts at rest with
% twisted springs; [x; y] is the center of mass of the first link.
q0 = [1; 1.5; pi/12; pi/4; pi/10; pi/6];

%_________________________Simulation_______________________________________
t = (0:speed/fps:T).';
options = odeset('RelTol', 1e-9, 'AbsTol', 1e-10);
switch method
    case 'IntMul',           [~, z] = ode45(@IntMulDyn, t, [q0; zeros(10, 1)], options);
    case 'Augmented',        [~, z] = ode45(@AugmntDyn, t, [q0; zeros(10, 1)], options);
    case 'Elimination',      [~, z] = ode45(@ElimDyn, t, [q0; zeros(6, 1)], options);
    case 'Embedding',        [~, z] = ode45(@EmbdDyn, t, [q0; 0; 0], options);
    case 'ModifiedLagrange', [~, z] = ode45(@ModlagDyn, t, [q0; 0; 0], options);
    case 'Greenwood',        [~, z] = ode45(@GreenDyn, t, [q0; zeros(6, 1)], options);
    case 'GibbsAppell',      [~, z] = ode45(@GibssDyn, t, [q0; 0; 0], options);
    case 'Kane',             [~, z] = ode45(@KaneDyn, t, [q0; 0; 0], options);
end
q = z(:, 1:6);

%_________________________Animation________________________________________
% Speed of the head link along its own axis (the robot moves head first,
% toward the free end of the first link).
v = zeros(numel(t), 1);
for k = 1:numel(t) - 1
    d = (q(k+1, 1:2) - q(k, 1:2))/(t(k+1) - t(k));
    v(k) = -d*[cos(q(k, 3)); sin(q(k, 3))];
end
v(end) = v(end-1);
names = struct('IntMul', 'integrated multipliers', 'Augmented', 'augmented', 'Elimination', 'elimination', ...
    'Embedding', 'embedding', 'ModifiedLagrange', 'modified Lagrange', 'Greenwood', 'Greenwood', ...
    'GibbsAppell', 'Gibbs-Appell', 'Kane', 'Kane');
readout = cell(numel(t), 1);
for k = 1:numel(t)
    readout{k} = sprintf('Snake robot (%s method, %g\\times speed)\nt = %4.1f s,  speed = %4.2f m/s', names.(method), speed, t(k), v(k));
end
opts.speed = speed; opts.readout = readout;
if saveVideo, opts.video = sprintf('snake_robot_%s.mp4', method); end
animateSnake(t, q, [l1 l2 l3 l4], opts);
