clc; clf; clear all; close all;

% =======================
% SIM CONFIG
% =======================

% snake params
global Nl l m J;
Nl = 10;
l = 0.09;
m = 1.56;
J = (1/3)*m*l^2;

% sim params
global T dt N_steps;
T = 60;
dt = 0.01;
N_steps = round(T/dt);
t = linspace(0,T,N_steps);

global int_method ref;
online_vis = true;
dynamics = @dyn635;                 % dyn635 / dyn815
controller = @co_fblin;              % co_fblin / co_adap
int_method = @RK4;                  % RK4 / euler_step
ref = @lat_ond;                     % lat_ond / eel_like
heading_law = @no;          % no / pd_heading / casc_heading / adap_heading
rand_init_pose = false;
rand_c_hat0 = false;
rand_l_hat0 = false;

% gait pattern
global alpha omega delta;
alpha = 0.052;                      % 5.2cm approx 16° in the complex model
omega = 1*120*pi/180;
delta = 50*pi/180;

% los settings
global lh_dist
lh_dist = 1*Nl*l;

% terrain
global ct cn cp c_hat maxc minc
minc = 0.01;
maxc = 100;
ct = 0.5;
cn = 3;
c_hat = [minc;minc]; %cn;cp
if rand_c_hat0
    c_hat(1) = rand()*(maxc - minc) + minc;     %cn
    c_hat(2) = rand()*maxc/(2*l);               %cp
end

% lyap analysis
global v vdot;
[v, vdot] = deal([0,0]);

% path
yref_fun = @(x) -15+15*sin(0.05*x+pi/2) + 1*sin(0.5*x);
yref_fun = @(x) heaviside(x-4);

% heading controller params
global vt_min iss_min;
vt_min = 0.02;
iss_min = 0.1;
cond = false;

% CO model params
global l1 l2 l1_hat l2_hat a a_hat;
l1 = 0.5; l2 = 20;
a = [1/l2; l1/l2];
a_hat = [0;0];
[l1_hat, l2_hat] = deal(NaN);

%% initial conditions ------------------------------------

init.angs = ones(Nl-1,1)*0;
init.theta = 0;
init.pos = [0;0];
if rand_init_pose
    init.angs = -l/4 + rand(Nl-1,1)*2*l/4;
    init.theta = -pi/4 + rand()*2*pi/4;
    init.pos = [0; rand()*2];
end
init.ang_vel = zeros(Nl-1,1);
init.theta_vel = 0;
init.vel = [0;0]; 

init.x0 = [init.angs; init.theta; init.pos; init.ang_vel; init.theta_vel; init.vel];
x = init.x0;


%% derived params

global nx nu e e_bar A D D_bar K V H;
nx = 2*Nl+4;
nu = Nl-1;
e = ones(Nl,1);
e_bar = ones(Nl-1,1);

%(2.34)
A = diag(ones(Nl, 1)) + diag(ones(Nl-1, 1), 1);
A = A(1:Nl-1,:);
D = diag(ones(Nl, 1)) + diag(-ones(Nl-1, 1), 1);
D = D(1:Nl-1,:);
D_bar = D'*inv(D*D');
K = A'*inv(D*D')*D;
V = A'*inv(D*D')*A;
H = triu(ones(Nl,Nl));


%% savearr sv

sv.x = zeros(N_steps,nx);
[sv.phi_ref, sv.phi_dot_ref, sv.phi, sv.phi_dot] = deal(zeros(N_steps,nu));
[sv.v, sv.vdot] = deal(zeros(N_steps,2));
sv.target = zeros(N_steps,2);
sv.heading = zeros(N_steps,1);
[sv.threfs, sv.phi0s] = deal(zeros(N_steps,4));
terrain.c_trues = zeros(N_steps,2);     %cn;cp
sv.c_hat = zeros(N_steps,2);            %cn;cp
[sv.l_hat, sv.a_hat] = deal(zeros(N_steps,2));
sv.cond = zeros(N_steps,1);

% =======================
% SIM LOOP
% =======================

last_percent = 0;
for k = 1:N_steps
    percent = floor(k/N_steps * 100);
    if percent > last_percent
        fprintf('Processing...%d%%\n', percent);
        last_percent = percent;
    end

    px = x(Nl+1);
    if px >= 30         %change terrain
        ct = 4;
        cn = 17;
    end
    assert(cn >= ct)    %for propulsion
    cp = (cn-ct)/(2*l);
    
    % control loop
    target = [px + lh_dist, yref_fun(px+lh_dist)];
    [heading, phi0, phi0_bar, threfs] = heading_law(t(k), x, target(2), cond);
    [phi_ref, phi_ref_dot, phi_ref_ddot] = ref(t(k));
    [phi_ref, phi_ref_dot, phi_ref_ddot, phi0s] = ...
        hybrid_add_phi0(phi_ref, phi_ref_dot, phi_ref_ddot, phi0, phi0_bar);
    u = controller(x, phi_ref, phi_ref_dot, phi_ref_ddot);
    x = int_method(x, u, dynamics, dt);

    % save stuff
    sv.x(k,:) = x;
    sv.phi(k,:) = x(1:Nl-1);
    sv.phi_dot(k,:) = x(Nl+3:end-3);
    sv.phi_ref(k,:) = phi_ref;
    sv.phi_dot_ref(k,:) = phi_ref_dot;
    sv.heading(k) = heading;
    sv.threfs(k,:) = threfs;
    terrain.path(k,:) = [px, yref_fun(px)];
    sv.target(k,:) = target;
    sv.phi0s(k,:) = phi0s;
    terrain.c_trues(k,:) = [cn,cp];
    sv.c_hat(k,:) = c_hat;
    sv.a_hat(k,:) = a_hat;
    sv.l_hat(k,:) = [l1_hat; l2_hat];
    sv.v(k,:) = v;
    sv.vdot(k,:) = vdot;
    sv.cond(k) = cond;
    
    % tresholds to run the heading controller
    vt = x(end-1);
    err = vecnorm(x(1:Nl-1) - phi_ref, 2, 1);
    cond = (err <= iss_min && vt > vt_min);
end

% =======================
% VIS LOOP
% =======================

speed_mult = 1;
target_speed = speed_mult;
speed_corrected = false;

if online_vis
    figure(1);
    handles.title = title('', 'FontSize', 8);
    hold on;
    
    % for closing figure
    set(gcf, 'KeyPressFcn', @(src, event) setappdata(src, 'keyPressed', event.Key));
    setappdata(gcf, 'keyPressed', '');

    % terrain
    handles.path = plot(terrain.path(:,1), terrain.path(:,2), 'k--');
    handles.trace = plot(NaN, NaN, '-', 'Color', [0.5 0.5 0.5]);

    % ref snek
    handles.reflinks = plot(NaN, NaN, '-', 'Color', [0.3 0.3 0.3]);
    handles.target_dir = plot(NaN, NaN, 'm-');
    handles.target_spot = scatter(NaN, NaN, 'm>');

    % snek
    handles.pos = scatter(NaN, NaN, 'r+');
    handles.thetaref_filt_dir = plot(NaN, NaN, 'm--');
    handles.heading_dir = plot(NaN, NaN, 'g-');
    handles.phi0_filt_dir = plot(NaN, NaN, 'b-');
    handles.links = plot(NaN, NaN, 'b-');
    handles.dots = scatter(NaN, NaN, 6, 'ro');
    
    grid minor; axis equal; hold off;

    for j=1:N_steps
        
        % handle eventual closing
        key = getappdata(gcf, 'keyPressed');
        if strcmp(key, 'q')
            disp('<q> pressed. Exiting visualization');
            break;
        end
        setappdata(gcf, 'keyPressed', '');
        
        % speed regulation
        tic;
        if speed_mult < target_speed
            speed_mult = speed_mult - dt*(speed_mult - target_speed);
        end
        k=max(1,min(round(j*speed_mult),N_steps));

        % snek
        phi = sv.x(k,1:Nl-1)';
        theta = sv.x(k,Nl);
        px = sv.x(k,Nl+1); 
        py = sv.x(k,Nl+2);

        ti = px*e - l*D_bar*e_bar;
        ni = py*e - D_bar*phi;
        
        % rotation matrix
        R = [cos(theta), -sin(theta);
             sin(theta),  cos(theta)];
        
        % pivot around px and py and re-offset back
        temp = [ti-px, ni-py]*R';
        X = temp(:,1) + px;
        Y = temp(:,2) + py;

        % plot the links
        pis = zeros(2*Nl,2);
        for i = 0:Nl-1
            pis(2*i+1,:) = [X(i+1), Y(i+1)] - l/2*[cos(theta), sin(theta)];
            pis(2*i+2,:) = [X(i+1), Y(i+1)] + l/2*[cos(theta), sin(theta)];
        end

        % update plots
        set(handles.trace, 'XData', sv.x(1:k,Nl+1), 'YData', sv.x(1:k,Nl+2));
        set(handles.links, 'XData', pis(:,1), 'YData', pis(:,2));
        set(handles.dots, 'XData', X, 'YData', Y);
        set(handles.pos, 'XData', px, 'YData', py);
        set(handles.heading_dir, 'XData', [px px + l*Nl*cos(theta)], ...
                                 'YData', [py py + l*Nl*sin(theta)]);
        set(handles.thetaref_filt_dir, 'XData', [px, px + lh_dist*cos(sv.threfs(k,2))], ...
                                       'YData', [py, py + lh_dist*sin(sv.threfs(k,2))]);
        set(handles.phi0_filt_dir, 'XData', px + l*Nl*cos(theta) + [0, +20*sv.phi0s(k,1)*cos(theta+pi/2)], ...
                                   'YData', py + l*Nl*sin(theta) + [0, +20*sv.phi0s(k,1)*sin(theta+pi/2)]);

        % ref snek
        phiref = sv.phi_ref(k,:)';
        ni = py*e - D_bar*phiref; % <-- instead of phi

        % pivot around px and py and re-offset back
        temp = [ti-px, ni-py]*R';
        X = temp(:,1) + px;
        Y = temp(:,2) + py;

        % plot the links
        pis = zeros(3*Nl,2);
        for i = 0:Nl-1
            pis(3*i+1,:) = [X(i+1), Y(i+1)] - l/2*[cos(theta), sin(theta)];
            pis(3*i+2,:) = [X(i+1), Y(i+1)] + l/2*[cos(theta), sin(theta)];
            pis(3*i+3,:) = [NaN NaN];
        end

        % update plot
        set(handles.reflinks, 'XData', pis(:,1), 'YData', pis(:,2));
        set(handles.target_dir, 'XData', [px, sv.target(k,1)], ...
                                'YData', [py, sv.target(k,2)]);
        set(handles.target_spot, 'XData', sv.target(k,1), 'YData', sv.target(k,2));
        set(handles.title, 'String', sprintf('>>%.2fx, t: %.2f s, c: [%.2f %.2f], chat: [%.2f %.2f], model: %s, contr: %s, heading:%s, int: %s', speed_mult, t(k), [terrain.c_trues(k,1); terrain.c_trues(k,2)], sv.c_hat(k,:), func2str(dynamics), func2str(controller), func2str(heading_law), func2str(int_method)));

        pause(dt);
        if j*dt >= 0.1 && ~speed_corrected % a little to stabilize
            target_speed = toc/dt;
            speed_corrected = true;
        end

        if k>=N_steps
            break;
        end        
    end
end

% =======================
% PLOTS
% =======================

figure(10);
subplot(2,1,1);
plot(t, sv.threfs);
legend('thetaref','filt','dot','ddot');
title('thetaref');

subplot(2,1,2);
plot(t, sv.phi0s);
legend('phi0','filt','dot','ddot');
title('phi0');

%phi, phiref
figure(2);
[r,c] = find_best_subplot_ratio(2*Nl-1);
for i=1:Nl-1
    subplot(r,c,i);
    plot(t, sv.phi_ref(:,i)); hold on;
    plot(t, sv.phi(:,i)); hold off;
    title(sprintf('phi%d',i))
    set(gca, 'XTick', []);
end
for i=1:Nl-1
    subplot(r,c,i+Nl);
    plot(t, sv.phi_dot_ref(:,i)); hold on;
    plot(t, sv.phi_dot(:,i)); hold off;
    title(sprintf('phidot%d',i))
    set(gca, 'XTick', []);
end

%phi - phiref
figure(3)
subplot(2,2,3)
tracking_rms = 0;
[r,c] = find_best_subplot_ratio(2*Nl-1);
for i=1:Nl-1
    subplot(r,c,i);
    diff = sv.phi(:,i) - sv.phi_ref(:,i);
    tracking_rms = tracking_rms + norm(diff);
    plot(t, diff);
    title(sprintf('phi%d, rms: %.3f', i, norm(diff)))
    set(gca, 'XTick', []);
end
for i=1:Nl-1
    subplot(r,c,i+Nl);
    diff = sv.phi_dot(:,i) - sv.phi_dot_ref(:,i);
    tracking_rms = tracking_rms + norm(diff);
    plot(t, diff);
    title(sprintf('phidot%d, rms: %.3f', i, norm(diff)))
    set(gca, 'XTick', []);
end

fprintf("tracking_rms: %f\n", tracking_rms);

%friction coeffs
figure(4);
subplot(2,3,1)
plot(t, terrain.c_trues(:,1), 'r--', DisplayName='cn'); hold on;
plot(t, terrain.c_trues(:,2), 'g--', DisplayName='cp');
plot(t, sv.c_hat(:,1), 'r-', DisplayName='cn\^');
plot(t, sv.c_hat(:,2), 'g-', DisplayName='cp\^'); hold off;
xlim([t(1) t(end)]);
xlabel('t');
legend();
title('friction coefficients');

%lyapunov
subplot(2,3,2)
plot(t, sv.v(:,1), DisplayName='v');
hold on;
plot(t, sv.vdot(:,1), DisplayName='vdot');
legend();
yline(0,'k--');
hold off;
ylim([-10,50]); %zoom
xlim([t(1) t(end)]);
xlabel('t');
title('friction lyapunov');

fprintf("[cn cp] tilde: %.2f\n", c_hat-[cn;cp]);
fprintf("v(end): %f\n", v(1));

%norm of phis togheter
subplot(2,3,3)
plot(t, vecnorm(sv.phi - sv.phi_ref, 2, 2), DisplayName='phi');
hold on;
plot(t, vecnorm(sv.phi_dot - sv.phi_dot_ref, 2, 2), DisplayName='phidot');
yline(iss_min, 'k--', DisplayName='iss min');
hold off;
xlim([t(1) t(end)]);
xlabel('t');
legend();
title('||error variables||');

%x-y plot + path
subplot(2,3,4)
plot(sv.x(:,Nl+1), sv.x(:,Nl+2), DisplayName='yx');
hold on;
plot(terrain.path(:,1), terrain.path(:,2), 'k--', DisplayName='ref');
hold off;
xlim([min(sv.x(:,Nl+1)), max(sv.x(:,Nl+1))]);
xlabel('x');
legend();
title('position');

%heading
subplot(2,3,5)
plot(t, sv.heading, DisplayName='theta');
hold on;
plot(t, sv.x(:,end-2), DisplayName='thetadot');
plot(t, sv.threfs(:,3), DisplayName='thetaref dot');
plot(t, sv.threfs(:,2), DisplayName='thetaref filt');
hold off;
xlim([t(1) t(end)]);
xlabel('t');
legend();
title('heading');

%vt vn
subplot(2,3,6)
vt = sv.x(:,end-1);
vn = sv.x(:,end);
plot(t, vt, DisplayName='vt');
hold on;
plot(t, vn, DisplayName='vn');
yline(vt_min, 'k--', DisplayName='vt min');
hold off;
xlim([t(1) t(end)]);
xlabel('t');
legend();
title('vt vn');

figure(5)
subplot(3,1,1)
plot(t, sv.l_hat);
hold on;
yline(l1);
yline(l2);
hold off;
xlim([t(1) t(end)]);
xlabel('t');
legend('l1\^','l2\^', 'l1', 'l2');
title('lambdas');

subplot(3,1,2)
plot(t, sv.a_hat);
hold on;
yline(a(1));
yline(a(2));
hold off;
xlim([t(1) t(end)]);
xlabel('t');
legend('a1\^','a2\^', 'a1=1/l2', 'a2=l1/l2');
title('a');

subplot(3,1,3)
area(t, 0.2*sv.cond, 'FaceColor', [0.9 0.9 1], 'EdgeColor', 'none');
hold on;
plot(t, sv.v(:,2), DisplayName='v');
plot(t, sv.vdot(:,2), DisplayName='vdot');
yline(0,'k--');
hold off;
xlim([t(1) t(end)]);
ylim([-0.2,0.2]);
xlabel('t');
legend();
title('lambdas lyapunov');

fprintf("[l1 l2] tilde: %.2f\n", [l1_hat; l2_hat] -[l1;l2]);
fprintf("[a1 a2] tilde: %.2f\n", a_hat-a);
fprintf("v(end): %f\n", v(2));

% subplot tool -----------------------------------------------

function [r,c] = find_best_subplot_ratio(n)
    done = false;
    [r,c] = deal(2);
    while ~done 
        if r*c >= n
            done = true;
        else
            if c >= 2*r  %keep aspect ratio <=2
                r = r+1;
            else
                c=c+1;
            end
        end
    end
end


%% dynamics using (6.35) ----------------------------

function dx = dyn635(x,u)
    global Nl m ct cn cp l1 l2;
    global A D D_bar e_bar nx;
    
    phi = x(1:Nl-1); theta = x(Nl);
    vphi = x(Nl+3:end-3); vtheta = x(end-2);
    vt = x(end-1); vn = x(end);
    
    dx = zeros(nx,1);
    dx(1:Nl-1) = vphi;
    dx(Nl) = vtheta;
    dx(Nl+1) = vt*cos(theta) - vn*sin(theta);
    dx(Nl+2) = vt*sin(theta) + vn*cos(theta);
    dx(Nl+3:end-3) = -cn/m * vphi + cp/m * vt * A*D' * phi + 1/m * (D*D')*u;
    dx(end-2) = -l1*vtheta + l2/(Nl-1) * vt * e_bar' * phi;
    dx(end-1) = -ct/m * vt + 2*cp/(Nl*m) * vn * e_bar' * phi + ...
                -cp/(Nl*m) * phi' * A*D_bar * vphi;
    dx(end) = -cn/m * vn + 2*cp/(Nl*m) * vt * e_bar' * phi;
end


%% dynamics using (8.15) ----------------------------

function dx = dyn815(x,u)
    global Nl m ct cn cp l1 l2;
    global A D D_bar e_bar nx;

    phi = x(1:Nl-1); theta = x(Nl);
    vphi = x(Nl+3:end-3); vtheta = x(end-2);
    vt = x(end-1); vn = x(end);

    eps = -2*(Nl-1)/(Nl*m) * cp/l2;

    dx = zeros(nx,1);
    dx(1:Nl-1) = vphi;
    dx(Nl) = vtheta;
    dx(Nl+1) = vt*cos(theta) - vn*sin(theta);
    dx(Nl+2) = vt*sin(theta) + vn*cos(theta);
    dx(Nl+3:end-3) = -cn/m * vphi + cp/m * vt * A*D' * phi + 1/m * (D*D')*u;
    dx(end-2) = -l1*vtheta + l2/(Nl-1) * vt * e_bar' * phi;
    dx(end-1) = -ct/m * vt + 2*cp/(Nl*m) * vn * e_bar' * phi - cp/(Nl*m) * phi' * A*D_bar * vphi;
    X = eps*(cn/m-l1);
    Y = -cn/m;
    dx(end) = X*vtheta + Y*vn;
end


%% integration formulas --------------------------------------

function x_next = euler_step(x,u,law,dt)
    x_next = x + dt*law(x,u);
end

function x_next = RK4(x,u,law,dt)
    k1 = law(x, u);
    k2 = law(x+dt/2*k1, u);
    k3 = law(x+dt/2*k2, u);
    k4 = law(x+dt*k3, u);
    
    x_next = x + dt/6*(k1+2*k2+2*k3+k4);
end


%% co-fb linearizing controller

function u = co_fblin(x, phi_ref, phi_ref_dot, phi_ref_ddot)
    global Nl m cn cp A D
    
    %TUNING ---------
    Kp = 3; Kd = 1;
    %----------------

    phi = x(1:Nl-1);
    phi_dot = x(Nl+3:end-3);
    vt = x(end-1);
    
    %fb linearizing control law
    mDDTinv = m*inv(D*D');
    u_bar = phi_ref_ddot - Kd * (phi_dot - phi_ref_dot) - Kp * (phi - phi_ref);
    u = mDDTinv * (u_bar + cn/m * phi_dot - cp/m * vt * A*D' * phi);
end


%% adaptive fb-lin controller ----------------------------------

function u = co_adap(x, phi_ref, phi_ref_dot, phi_ref_ddot)
    global Nl l m A D;
    global dt int_method cn cp c_hat maxc minc;

    %TUNING ---------
    Kp = 3; Kd = 1;
    Gamma = diag([30,10000]);
    %----------------

    %(14b-21-22)
    Az = [-Kd*eye(Nl-1), -Kp*eye(Nl-1);
        eye(Nl-1), zeros(Nl-1, Nl-1)];
    Q = [Kd*eye(Nl-1), zeros(Nl-1,Nl-1);
         zeros(Nl-1,Nl-1), Kp*eye(Nl-1)];
    P = lyap(Az', Q);

    %some checks
    assert(isequal(P,P'), 'P is not symmetric');
    assert(isequal(Q,Q'), 'Q is not symmetric');
    assert(all(eig(P)>0), 'P is not positive def');
    assert(all(eig(Q)>0), 'Q is not positive def');
    assert(all(abs(Az'*P+P*Az - -Q) < 1e-12, 'all'), ...
        'riccati equation is not satisfied');

    phi = x(1:Nl-1);
    phi_dot = x(Nl+3:end-3);
    vt = x(end-1);
    
    % adaptive laws
    z = [phi_dot - phi_ref_dot; phi - phi_ref];
    Y = [-1/m * phi_dot, vt/m * A*D' * phi];
    F = [Y; zeros(Nl-1,2)];
    
    adap_law = (@(c,z) 2*Gamma*F'*P*z);
    c_hat = int_method(c_hat, z, adap_law, dt);

    % parameter projection
    c_hat(1) = max(minc, min(c_hat(1), maxc));              %cn
    c_hat(2) = max(0, min(c_hat(2), (maxc-minc)/(2*l)));    %cp

     % lyapunov analysis
    global v vdot;
    c_tilde = [cn;cp] - c_hat;
    Gammainv = inv(Gamma);
    v(1) = z'*P*z + 0.5*c_tilde'*Gammainv*c_tilde;
    vdot(1) = -z'*Q*z;
   
    %fb linearizing control law
    mDDTinv = m*inv(D*D');
    cn_hat = c_hat(1);
    cp_hat = c_hat(2);
    
    u_bar = phi_ref_ddot - Kd * (phi_dot - phi_ref_dot) - Kp * (phi - phi_ref);
    u = mDDTinv * (u_bar + cn_hat/m * phi_dot - cp_hat/m * vt * A*D' * phi);
end


%% pd-heading controllers ----------------------------------

function [heading, phi0, phi0_bar, threfs] = pd_heading(t,x,yref,~)
    global Nl lh_dist ct cn cp;

    % initialization
    persistent thref_filt thref_dot thref_ddot;
    if isempty(thref_filt)
        [thref_filt, thref_dot, thref_ddot] = deal(0);
    end

    %TUNING ---------
    Kp = 0.4;
    Kd = 0.1;
    %----------------

    % LOS
    py = x(Nl+2);
    thref = -atan2(py - yref, lh_dist);

    %TUNING ------
    xi = 1;
    deltaT = 4;
    %-------------

    % theta_ref derivatives
    [thref_filt, thref_dot, thref_ddot] = ...
        ref_filter([thref_filt, thref_dot, thref_ddot]', thref, xi, deltaT);
    threfs = [thref, thref_filt, thref_dot, thref_ddot];
    
    % control law
    heading = x(Nl);
    heading_dot = x(end-2);
    phi0 = - Kp*(heading - thref_filt) - Kd*(heading_dot - thref_dot);
    
    sat = Nl/(2*(Nl-1))*sqrt(ct*cn)/cp; %7.25
    phi0 = max(min(phi0,sat),-sat);
    phi0_bar = [0;0;0];
end


%% cascaded heading controller --------------------------

function [heading, phi0, phi0_bar, threfs] = casc_heading(t,x,yref,~)
    global Nl lh_dist ref;
    global l1 l2 vt_min e_bar ct cn cp;

    % initialization
    persistent thref_filt thref_dot thref_ddot;
    if isempty(thref_filt)
        [thref_filt, thref_dot, thref_ddot] = deal(0);
    end

    %TUNING ---------
    ktheta = 0.05;
    %----------------

    % LOS
    py = x(Nl+2);
    thref = -atan2(py - yref, lh_dist);

    %TUNING ------
    xi = 1;
    deltaT = 4;
    %-------------

    % theta_ref derivatives
    [thref_filt, thref_dot, thref_ddot] = ...
        ref_filter([thref_filt, thref_dot, thref_ddot]', thref, xi, deltaT);
    threfs = [thref, thref_filt, thref_dot, thref_ddot];

    heading = x(Nl);
    heading_dot = x(end-2);
    vt = x(end-1);

    % control law
    phi0 = 1/(vt*l2) * ...
        (thref_ddot + l1*thref_dot - ktheta*(heading - thref_filt));
    % but
    if vt <= vt_min
        phi0 = 0;
    end

    sat = Nl/(2*(Nl-1))*sqrt(ct*cn)/cp; %7.25
    phi0 = max(-sat, min(phi0, sat));

    [phiref_bar, phiref_bar_dot, phiref_bar_ddot] = ref(t);
    phi0_bar(1) = -e_bar'*phiref_bar/(Nl-1);        %phi0_bar
    phi0_bar(2) = -e_bar'*phiref_bar_dot/(Nl-1);    %phi0_bar_dot
    phi0_bar(3) = -e_bar'*phiref_bar_ddot/(Nl-1);   %phi0_bar_ddot
end

%% adaptive heading controller --------------------------

function [heading, phi0, phi0_bar, threfs] = adap_heading(t,x,yref,cond)
    global Nl lh_dist ref;
    global a e_bar;
    global l1_hat l2_hat a_hat dt;

    % initialization
    persistent thref_filt thref_dot thref_ddot;
    if isempty(thref_filt)
        [thref_filt, thref_dot, thref_ddot] = deal(0);
    end

    %TUNING ---------
    Lambda = 0.1;
    Kd = 0.01;
    Gamma = diag([10,1]);
    %----------------

    % LOS
    py = x(Nl+2);
    thref = -atan2(py - yref, lh_dist);

    %TUNING ------
    xi = 1;
    deltaT = 4;
    %-------------
    
    [thref_filt, thref_dot, thref_ddot] = ...
        ref_filter([thref_filt, thref_dot, thref_ddot]', thref, xi, deltaT);
    threfs = [thref, thref_filt, thref_dot, thref_ddot];

    heading = x(Nl);
    heading_dot = x(end-2);
    vt = x(end-1);
        
    % filtered tracking error
    s = (heading_dot - thref_dot) + Lambda*(heading - thref_filt);
    thetaR_dot = thref_dot - Lambda*(heading - thref_filt);
    thetaR_ddot = thref_ddot - Lambda*(heading_dot - thref_dot);
    
    % adaptive part
    global v vdot;
    a_tilde = a - a_hat;
    Gammainv = inv(Gamma);
    v(2) = 0.5*a(1)*s^2 + 0.5*a_tilde'*Gammainv*a_tilde;
    vdot(2) = -Kd*s^2;
    
    Y = [thetaR_ddot, heading_dot];
    a_hat_dot = @(a,s) -Gamma*Y'*s;

    if cond
        a_hat = euler_step(a_hat, s, a_hat_dot, dt);
    
        try
            l2_hat = 1/a_hat(1);        %try to avoid 1/0
            l1_hat = a_hat(2)*l2_hat;
        end
        
        u_bar = Y*a_hat - Kd*s;
        phi0 = 1/vt * u_bar;
    else
        phi0 = 0;
    end

    [phiref_bar, phiref_bar_dot, phiref_bar_ddot] = ref(t);
    phi0_bar(1) = -e_bar'*phiref_bar/(Nl-1);        %phi0_bar
    phi0_bar(2) = -e_bar'*phiref_bar_dot/(Nl-1);    %phi0_bar_dot
    phi0_bar(3) = -e_bar'*phiref_bar_ddot/(Nl-1);   %phi0_bar_ddot
end

%% no heading control ----------------------------------

function [heading, phi0, phi0_bar, threfs] = no(t,x,yref,~)
    global Nl;

    heading = mean(x(1:Nl));
    threfs = [0, 0, 0, 0];
    phi0 = 0;
    phi0_bar = [0;0;0];
end


%% create gait swimming reference

function [phi_ref, phi_ref_dot, phi_ref_ddot] = lat_ond(t)
    global alpha omega delta Nl;
    
    ii = (1:Nl-1)';
    phi_ref = alpha*sin(omega*t + (ii-1)*delta);
    phi_ref_dot = omega*alpha*cos(omega*t + (ii-1)*delta);
    phi_ref_ddot = -(omega^2)*alpha*sin(omega*t + (ii-1)*delta);
end

function [phi_ref, phi_ref_dot, phi_ref_ddot] = eel_like(t)
    global alpha omega delta Nl;
    
    ii = (1:Nl-1)';
    g = (Nl-ii)/(Nl+1);
    phi_ref = g .* (alpha*sin(omega*t + (ii-1)*delta));
    phi_ref_dot = g .* (omega*alpha*cos(omega*t + (ii-1)*delta));
    phi_ref_ddot = g .* (-(omega^2)*alpha*sin(omega*t + (ii-1)*delta));
end

%% reference filters

% 2rd order reference filter ----------------------------------

function [x, xdot] = ref_filter_2nd(x, r, xi, deltaT)
    global dt;

    wc = 2*pi/deltaT;
    F = [0, 1; -wc^2, -2*xi*wc];
    tf = @(x,r) F*x + [0; wc^2]*r;
    
    z = euler_step(x, r, tf, dt);
    x=z(1); xdot=z(2);
end


% 3rd order reference filter ----------------------------------

function [x, xdot, xddot] = ref_filter(x, r, xi, deltaT)
    global dt;
    
    wc = 2*pi/deltaT;
    F = [0, 1, 0; 0, 0, 1; -wc^3, -(2*xi+1)*(wc^2), -(2*xi+1)*wc];
    tf = @(x,r) F*x + [0; 0; wc^3]*r;
    
    z = euler_step(x, r, tf, dt);
    x=z(1); xdot=z(2); xddot=z(3);
end


% phi0, phi0_dot, phi0_ddot ----------------------------------

function [phi_ref, phi_ref_dot, phi_ref_ddot, phi0s] = hybrid_add_phi0(phi_ref, phi_ref_dot, phi_ref_ddot, phi0, phi0_bar)

    % initialization
    persistent phi0filt phi0filt_dot phi0filt_ddot;
    if isempty(phi0filt)
        [phi0filt, phi0filt_dot, phi0filt_ddot] = deal(0);
    end
    
    %TUNING ------
    xi = 1;
    deltaT = 0.1;
    %-------------

    [phi0filt, phi0filt_dot, phi0filt_ddot] = ...
        ref_filter([phi0filt, phi0filt_dot, phi0filt_ddot]', phi0, xi, deltaT);    
    
    phi_ref = phi_ref + phi0filt + phi0_bar(1);
    phi_ref_dot = phi_ref_dot + phi0filt_dot + phi0_bar(2);
    phi_ref_ddot = phi_ref_ddot + phi0filt_ddot + phi0_bar(3);

    phi0s = [phi0, phi0filt, phi0filt_dot, phi0filt_ddot];
end

function [phi_ref, phi_ref_dot, phi_ref_ddot, phi0s] = dont_filter_phi0(phi_ref, phi_ref_dot, phi_ref_ddot, phi0, phi0_bar)
    phi_ref = phi_ref + phi0;
    phi0s = [phi0, 0, 0, 0];
end


% THE END ----------------------------------------------------