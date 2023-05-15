clear;
close all;

P_B = 2250 * 745.7; % Convert HP to W
V_r = 2300; % Rated voltage, V

V_B = V_r/sqrt(3);
I_B = P_B/3/V_B;

P = 4; % Number of poles
omega_b = 60 * 2 * pi;

T_B = P_B / (2/P * omega_b);

rpm = 1786; % Rated mechanical speed, rpm
J = 63.87; % Inertia constant, J*s^2 = kg*m

Z_B = V_B / I_B;

r_s = 0.029/Z_B; % Stator winding resistance
r_r_p = 0.022/Z_B;
X_ls = 0.226/Z_B; % Stator leakage reactance  at rated freq
X_M = 13.04/Z_B; % Stator q-axis reactance at rated freq
X_lr_p = 0.226/Z_B; % Stator d-axis reactance at rated freq

H = 0.5 * (2/P) * J * omega_b / T_B;

v_in = 1;

% Calculate necessary parameters
X_aq = (1/X_M + 1/X_ls + 1/X_lr_p)^-1;
X_ad = X_aq;

%% Free acceleration
param = [
        omega_b
        r_s
        r_r_p
        X_ls
        X_lr_p
        X_M
        H
        0
        v_in
    ];
    
    [st,sy] = rk4(@(t,y) odefun(t,y,param), [0,5],[0,0,0,0,0,0,0,0],1e-4);
        
    y_ss = sy(:,end);
    y_low = zeros(size(y_ss));
    
    
%% Run for all torques


    torques = [0:.05:8];
    slips = zeros(size(torques));
    slips_l = zeros(size(torques));
    ends = false;
    endl = false;

for T_i = 1:numel(torques)
    disp(T_i/numel(torques) * 100 + "%")
    param(8) = torques(T_i);

    if ~ends
        [st,sy] = rk4(@(t,y) odefun(t,y,param), [0,10],y_ss,1e-4);
    end
    
    if ~endl
        [lt,ly] = rk4(@(t,y) odefun(t,y,param), [0,10],y_low,1e-4);
    end

    
    % Extract solutions
    psi_qs = sy(1,:);
    psi_ds = sy(2,:);
    psi_0s = sy(3,:);
    psi_qr_p = sy(4,:);
    psi_dr_p = sy(5,:);
    psi_0r_p = sy(6,:);
    omega_r = sy(7,:);
    theta_r = sy(8,:);
    
    psi_qs_l = ly(1,:);
    psi_ds_l = ly(2,:);
    psi_qr_p_l = ly(4,:);
    psi_dr_p_l = ly(5,:);
    omega_r_l = ly(7,:);

    % Calculate dependent quantities
    psi_mq = X_aq * (psi_qs / X_ls + psi_qr_p / X_lr_p);
    psi_md = X_ad * (psi_ds / X_ls + psi_dr_p / X_lr_p);
    
    psi_mq_l = X_aq * (psi_qs_l / X_ls + psi_qr_p_l / X_lr_p);
    psi_md_l = X_aq * (psi_ds_l / X_ls + psi_dr_p_l / X_lr_p);
    
    i_qs = 1/X_ls * (psi_qs - psi_mq);
    i_ds = 1/X_ls * (psi_ds - psi_md);
    
    i_qs_l = 1/X_ls * (psi_qs_l - psi_mq_l);
    i_ds_l = 1/X_ls * (psi_ds_l - psi_md_l);

    T_e = psi_ds .* i_qs - psi_qs .* i_ds; %Torque eqn in rotor ref. frame
    
    T_e_l = psi_ds_l .* i_qs_l - psi_qs_l .* i_ds_l;
        
    if ~ends
        slip = (omega_b - omega_r)/omega_b;
        if slip(end) < 1
            slips(T_i) = slip(end);
        elseif ~ends
            slips(T_i:end) = NaN;
            ends = true;
        end
    end
    
    if ~endl
        slip_l = (omega_b - omega_r_l)/omega_b;
        if slip_l(end)<1
            slips_l(T_i) = slip_l(end);
        elseif ~endl
            slips_l(T_i:end) = NaN;
            endl = true;
        end
    end
end

figure;
hold on;
plot(slips,torques,'-','linewidth',2)
plot(slips_l,torques,'-','linewidth',2)

format('Slip','Torque')

% Plot
figure;
hold on;
yline(torques(T_i)*T_B,'k--','linewidth',2);
plot(st,T_e*T_B,'linewidth',2);
format("t","Torque")
title("Torque")

%figure;
%plot(omega_r*60/(2*pi)*(2/P), T_e*T_B)

figure;
hold on;
plot(st, slip,'linewidth',2);
plot(lt,(omega_b-ly(7,:))/omega_b,'linewidth',2)

format("t","Slip")
title("Slip")

function dydt = odefun(t,y,param)

    omega_b = param(1);
    r_s = param(2); % Stator winding resistance
    r_r_p = param(3);
    X_ls = param(4); % Stator leakage reactance at rated freq
    X_lr_p = param(5); % Field winding leakage reactance at rated freq (ref to stator side)
    X_M = param(6);
    H = param(7);
    
    T_I = param(8); % mechanical torque    
    
    v_in = param(9);
    
    X_aq = (1/X_M + 1/X_ls + 1/X_lr_p)^-1;
    X_ad = X_aq;

    % def y
    psi_qs = y(1);
    psi_ds = y(2);
    psi_0s = y(3);
    psi_qr_p = y(4);
    psi_dr_p = y(5);
    psi_0r_p = y(6);
    omega_r = y(7);
    theta_r = y(8);
    
    K_s_r = (2/3) * [
        cos(theta_r) cos(theta_r - (2*pi/3)) cos(theta_r + (2*pi/3))
        sin(theta_r) sin(theta_r - (2*pi/3)) sin(theta_r + (2*pi/3))
        0.5 0.5 0.5
    ];
    
    v_abcs = v_in * [cos(t * omega_b); cos(t * omega_b - (2*pi/3)); cos(t*omega_b+(2*pi/3))];
    v_qd0s_r = (K_s_r * v_abcs);
    
    v_qs = v_qd0s_r(1);
    v_ds = v_qd0s_r(2);
    v_0s = v_qd0s_r(3);
    
    %Shorted rotor windings
    v_qr_p = 0;
    v_dr_p = 0;
    v_0r_p = 0;
    
    psi_mq = X_aq * (psi_qs / X_ls + psi_qr_p / X_lr_p);
    psi_md = X_ad * (psi_ds / X_ls + psi_dr_p / X_lr_p);
    
    i_qs = 1/X_ls * (psi_qs - psi_mq);
    i_ds = 1/X_ls * (psi_ds - psi_md);
    i_0s = 1/X_ls * (psi_0s);
    i_qr_p = 1/X_lr_p * (psi_qr_p - psi_mq);
    i_dr_p = 1/X_lr_p * (psi_dr_p - psi_md);
    i_0r_p = 1/X_lr_p * (psi_0r_p);

    
    T_e = psi_ds * i_qs - psi_qs * i_ds; %rotor ref. frame
    %T_e = psi_qr_p * i_dr_p - psi_dr_p * i_qr_p; %other ref frame
    
    % Rotor reference frame
    omega = omega_r;
    %Synchronous ref frame
    %omega = omega_b;
    
    dydt = [
    omega_b * (v_qs - omega / omega_b * psi_ds + r_s / X_ls * (psi_mq - psi_qs))
    omega_b * (v_ds + omega / omega_b * psi_qs + r_s / X_ls * (psi_md - psi_ds))
    omega_b * (v_0s - r_s / X_ls * psi_0s)
    omega_b * (v_qr_p - (omega - omega_r)/omega_b * psi_dr_p + r_r_p / X_lr_p * (psi_mq-psi_qr_p))
    omega_b * (v_dr_p + (omega - omega_r)/omega_b * psi_qr_p + r_r_p / X_lr_p * (psi_md-psi_dr_p))
    omega_b * (v_0r_p - r_r_p / X_lr_p * psi_0r_p)
    omega_b / (2*H) * (T_e - T_I)
    omega_r
    ];
    
end

function [t,y] = rk2(odefun, tspan, y0, h)
% Custom Implementation of rk2

    t = tspan(1):h:tspan(end);
    y = zeros(numel(y0),numel(t));
    y(:,1) = y0;
    
    alpha = 0.5;
    beta = alpha;
    gamma1 = 1-1/(2*alpha);
    gamma2 = 1/(2*alpha);
    for i = 1:numel(t)-1
        k1 = h * odefun(t(i),y(:,i));
        k2 = h * odefun(t(i) + alpha * h,y(:,i) + beta * k1);
        y(:,i+1) = y(:,i) + gamma1 * k1 + gamma2 * k2;
    end
end

function [t,y] = rk4(odefun, tspan, y0, h)
% Custom Implementation of rk4 from PS3
    t = tspan(1):h:tspan(end);
    y = zeros(numel(y0),numel(t));
    y(:,1) = y0;
    
    for i = 1:numel(t)-1
        k1 = h*odefun(t(i),y(:,i));
        k2 = h*odefun(t(i) + h/2,y(:,i)+k1/2);
        k3 = h*odefun(t(i) + h/2,y(:,i)+k2/2);
        k4 = h*odefun(t(i) + h,y(:,i)+k3);
        y(:,i+1) = y(:,i) + k1/6 + 1/3 * (k2+k3) + k4/6;
    end
end

function format(xl,yl)
    xlabel(xl,'Interpreter','Latex','FontSize',20)
    ylabel(yl,'Interpreter','Latex','FontSize',20)
    
    box on;
    set(gcf,'color','w')
    set(gca,'FontName','Times','FontSize',20)
    set(gca,'TickLabelInterpreter','latex')
    xa = gca;
    xa.TickLength = [0.025,0.025];
    xa.LineWidth = 1.5;
end