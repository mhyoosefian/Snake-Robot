function animateSnake(t, q, L, opts)
% ANIMATESNAKE  3D animation of the four-link wheeled snake robot.
%
%   animateSnake(t, q, L)        plays the motion in a figure window.
%   animateSnake(t, q, L, opts)  with options (all optional):
%       opts.video   = 'snake.mp4'  also saves the animation as an MP4 file
%       opts.speed   = 1            playback speed (2 = twice as fast)
%       opts.readout = {}           N x 1 cell array of text shown at each frame
%
%   t : N x 1 time vector (uniformly sampled, s)
%   q : N x 6 generalized coordinates [x y theta1 theta2 theta3 theta4], where
%       (x, y) is the center of mass of the first link (the head link)
%   L : link lengths [l1 l2 l3 l4] (m)

if nargin < 4, opts = struct; end
if ~isfield(opts, 'video'), opts.video = ''; end
if ~isfield(opts, 'speed'), opts.speed = 1; end
if ~isfield(opts, 'readout'), opts.readout = {}; end
N = numel(t); dt = t(2) - t(1);

rad = [.17 .15 .12 .08];          % body radius of each link (m)
wheelR = .09;                     % wheel radius (m)
zc = wheelR + .06;                % height of the body axis above the floor
green = [.33 .62 .30; .26 .50 .24; .33 .62 .30; .26 .50 .24];

% ---------------- Figure ----------------
fig = figure('Color', 'w', 'Position', [100 100 1280 720], 'Name', 'Snake robot');
ax = axes(fig, 'Position', [0 0 1 1]); hold(ax, 'on'); axis(ax, 'equal', 'off');
P = joints(q(1, :), L);
xs = zeros(N, 1); ys = xs;
for k = 1:N, Pk = joints(q(k, :), L); xs(k) = mean(Pk(1, :)); ys(k) = mean(Pk(2, :)); end
pad = 6;
xl = [min(xs) - pad, max(xs) + pad]; yl = [min(ys) - pad, max(ys) + pad];
axis(ax, [xl yl 0 1.5]);
view(ax, -20, 38); camva(ax, 3.1);
% The camera follows the snake at a fixed offset.
target = @(k) [xs(k) ys(k) .2];
camtarget(ax, target(1)); offset = campos(ax) - camtarget(ax);

% Floor with a 1 m grid
surf(ax, xl([1 2; 1 2]), yl([1 1; 2 2]), zeros(2), 'FaceColor', [.95 .95 .93], 'EdgeColor', 'none', 'FaceLighting', 'none');
for g = ceil(xl(1)):floor(xl(2)), plot3(ax, [g g], yl, [1 1]*1e-3, 'Color', [.86 .86 .84]); end
for g = ceil(yl(1)):floor(yl(2)), plot3(ax, xl, [g g], [1 1]*1e-3, 'Color', [.86 .86 .84]); end

% Wheel tracks: the path of each link's center of mass (where its wheels are)
trackColor = [.55 .45 .35];
for i = 1:4
    track(i) = plot3(ax, nan, nan, nan, 'Color', trackColor, 'LineWidth', 1); %#ok<AGROW>
end

% Links: tapered capsules with wheels at the center of mass, drawn in each
% link's own frame (x along the link) and moved with hgtransform.
[cx, cy, cz] = cylinder(1, 32);
[sx, sy, sz] = sphere(24);
for i = 1:4
    tf(i) = hgtransform(ax); %#ok<AGROW>
    r = rad(i);
    surf(ax, cz*L(i), cx*r, cy*r + zc, 'Parent', tf(i), 'FaceColor', green(i, :), 'EdgeColor', 'none');
    surf(ax, sx*r, sy*r, sz*r + zc, 'Parent', tf(i), 'FaceColor', green(i, :), 'EdgeColor', 'none');
    % Two wheels on an axle normal to the link, at its center of mass
    [wx, wy, wz] = cylinder(wheelR, 24);
    for s = [-1 1]
        surf(ax, wx + L(i)/2, s*(r + .02) + s*.06*wz, wy + wheelR, 'Parent', tf(i), 'FaceColor', [.15 .15 .17], 'EdgeColor', 'none');
        fill3(ax, wheelR*cos(linspace(0, 2*pi, 25)) + L(i)/2, s*(r + .08)*ones(1, 25), wheelR*sin(linspace(0, 2*pi, 25)) + wheelR, ...
            [.45 .45 .5], 'Parent', tf(i), 'EdgeColor', 'none');
    end
    if i < 4
        % Joint with a torsional spring: a gold coil on top of the joint
        a = linspace(0, 8*pi, 200); rc = .07;
        plot3(ax, L(i) + rc*cos(a), rc*sin(a), zc + rad(i+1) + .01 + .1*a/(8*pi), 'Parent', tf(i), ...
            'Color', [.85 .62 .12], 'LineWidth', 2.5);
    else
        surf(ax, sx*r + L(i), sy*r, sz*r + zc, 'Parent', tf(i), 'FaceColor', green(i, :), 'EdgeColor', 'none');
    end
end
% Head: eyes at the free end of the first link
for s = [-1 1]
    surf(ax, sx*.045 - .07, sy*.045 + s*.09, sz*.045 + zc + .12, 'Parent', tf(1), 'FaceColor', 'w', 'EdgeColor', 'none');
    surf(ax, sx*.022 - .1, sy*.022 + s*.095, sz*.022 + zc + .135, 'Parent', tf(1), 'FaceColor', 'k', 'EdgeColor', 'none');
end
light(ax, 'Position', [1 -1 2], 'Style', 'infinite'); light(ax, 'Position', [-1 1 .6], 'Style', 'infinite');
lighting(ax, 'gouraud'); material(ax, 'dull');
label = annotation(fig, 'textbox', [.03 .82 .6 .15], 'String', '', 'EdgeColor', 'none', ...
    'FontSize', 15, 'FontName', 'Helvetica', 'VerticalAlignment', 'top', 'Color', [.15 .15 .2]);

% ---------------- Playback / recording ----------------
if ~isempty(opts.video)
    vw = VideoWriter(opts.video, 'MPEG-4'); vw.FrameRate = round(opts.speed/dt); vw.Quality = 95; open(vw);
end
cm = zeros(2, N, 4);
for k = 1:N
    P = joints(q(k, :), L);
    for i = 1:4
        th = q(k, 2 + i);
        M = makehgtform('translate', [P(:, i); 0], 'zrotate', th);
        set(tf(i), 'Matrix', M);
        cm(:, k, i) = P(:, i) + L(i)/2*[cos(th); sin(th)];
        set(track(i), 'XData', cm(1, 1:k, i), 'YData', cm(2, 1:k, i), 'ZData', 2e-3*ones(1, k));
    end
    camtarget(ax, target(k)); campos(ax, target(k) + offset);
    if ~isempty(opts.readout), set(label, 'String', opts.readout{k}); end
    drawnow;
    if isempty(opts.video)
        pause(max(0, dt/opts.speed - .002));
        if ~isvalid(fig), return; end
    else
        writeVideo(vw, getframe(fig));
    end
end
if ~isempty(opts.video), close(vw); fprintf('Saved %s\n', opts.video); end
end

function P = joints(q, L)
% Start point of each link: the head S, then the joints A, B and C.
x = q(1); y = q(2); th = q(3:6);
P = zeros(2, 4);
P(:, 1) = [x; y] - L(1)/2*[cos(th(1)); sin(th(1))];
for i = 2:4
    P(:, i) = P(:, i-1) + L(i-1)*[cos(th(i-1)); sin(th(i-1))];
end
end
