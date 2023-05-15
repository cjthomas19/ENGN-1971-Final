clear;
close all;

P_B = 0.25 * 745.7; % Convert HP to W
V_r = 110; % Rated voltage, Vrms
V_B = V_r*sqrt(2);

I_B = P_B/V_B;

P = 4; % Number of poles
omega_b = 60 * 2 * pi;

T_B = P_B / (2/P * omega_b);

%rpm = 1710; % Rated mechanical speed, rpm
J = 1.46e-2; % Inertia constant, J*s^2 = kg*m^2

Z_B = V_B / I_B;

r_s = 2.02/Z_B; % Stator winding resistance
r_r_p = 4.12/Z_B;
X_ls = 2.79/Z_B;
X_ms = 66.8/Z_B;
X_lr_p = 2.12/Z_B;

H = 0.5 * (2/P) * J * omega_b / T_B;

v_in = 1;

%% Free acceleration
param = [
        omega_b
        r_s
        r_r_p
        X_ls
        X_lr_p
        X_ms
        H
        0
        v_in
        0
];

%%
[st,sy] = rk4(@(t,y) simulate(t,y,param,[],{}), [0,4],[0,0,0,0,omega_b,0],1e-4);

plot(st,sy(5,:))

i_qs = sy(1,:);
i_ds = sy(2,:);
i_qr_p = sy(3,:);
i_dr_p = sy(4,:);
omega_r = sy(5,:);
theta_r = sy(6,:);
T_e = X_ms * (i_qs .* i_dr_p - i_ds .* i_qr_p); %rotor ref. frame
figure;
plot(st,T_e)
y0 = sy(:,end);
%% 
[t,y] = rk4(@(t,y) simulate(t,y,param,[10,8],{@(t,y) 1,@(t,y) t>2}), [0,10],y0,1e-4);
figure;
plot(t,(60*2*pi-y(5,:))/60/2/pi)

%%

torques = 0:.1:4;
slips = zeros(size(torques));
slips2 = zeros(size(torques));
ends = false;

for T_i = 1:numel(torques)
    disp(T_i/numel(torques) * 100 + "%")
    param(8) = torques(T_i);
    if ~ends
        param(10) = 1;
        [st,sy] = rk4(@(t,y) odefun(t,y,param),[0,10],y0,1e-4);
        param(10) = 0;
        [t2,y2] = rk4(@(t,y) odefun(t,y,param),[0,10],y0,1e-4);
    end
    
    % Extract solutions
    i_qs = sy(1,:);
    i_ds = sy(2,:);
    i_qr_p = sy(3,:);
    i_dr_p = sy(4,:);
    omega_r = sy(5,:);
    theta_r = sy(6,:);

    T_e = X_ms * (i_dr_p .* i_qs - i_qr_p .* i_ds); %Torque eqn in arbitrary ref. frame
    T_e2 = X_ms * (y2(4,:).*y2(1,:) - y2(3,:).*y2(2,:));
    omega_r2 = y2(5,:);
    if ~ends
        slip = (omega_b - omega_r)/omega_b;
        slip2 = (omega_b-omega_r2)/omega_b;
        if slip(end) < 1
            slips(T_i) = slip(end);
        else
            slips(T_i) = NaN;
        end
        if slip2(end) < 1
            slips2(T_i) = slip2(end);
        else
            slips2(T_i) = NaN;
        end
    end
end
%%
figure;
hold on;
plot(slips,torques,'-','linewidth',2)
plot(slips2,torques,'-','linewidth',2)
format_plot("Slip","Torque (N$\cdot$m)")
legend(["Single-Phase","Two-Phase"],'interpreter','latex')

%%
function dydt = simulate(t,y,param,ps,vs)
    arguments
        t double
        y double
        param double
        ps (:,1) double = []
        vs cell = {@(t,y) 0}
    end
    
    for i = 1:numel(ps)
        param(ps(i)) = vs{i}(t,y);    
    end
    dydt = odefun(t,y,param);
end

function dydt = odefun(t,y,param)

    omega_b = param(1);
    r_s = param(2); % Stator winding resistance
    r_r_p = param(3);
    X_ls = param(4); % Stator leakage reactance at rated freq
    X_lr_p = param(5); % Field winding leakage reactance at rated freq (ref to stator side)
    X_ms = param(6);
    H = param(7);
    T_I = param(8); % mechanical torque    
    v_in = param(9);
    ss = param(10);
    % def y
    i_qs = y(1);
    i_ds = y(2);
    i_qr_p = y(3);
    i_dr_p = y(4);
    omega_r = y(5);
    theta_r = y(6);
    
    theta = 0;
    
    K_s = [
        cos(theta)   sin(theta)
        sin(theta)  -cos(theta)
    ];
   
    v_abs = v_in * [cos(t * omega_b);cos(t * omega_b-pi/2)];
            
    v_qds_r = (K_s * v_abs);
    
    v_qs = v_qds_r(1);
    v_ds = v_qds_r(2);
    
    
    %Shorted-circuit:
    if ss
        i_qs = 0;
    end
   
    X_eq = X_ls * X_lr_p + X_ls * X_ms + X_lr_p * X_ms;
    
%     i_qs = omega_b*(X_lr_p * l_qs + X_ms * l_qs - X_ms * l_qr_p)/X_eq;
%     i_ds = omega_b*(X_lr_p * l_ds + X_ms * l_ds - X_ms * l_dr_p)/X_eq;
%     i_qr_p = omega_b*(X_ls * l_qr_p - X_ms * l_qs + X_ms * l_qr_p)/X_eq;
%     i_dr_p = omega_b*(X_ls * l_dr_p - X_ms * l_ds + X_ms * l_dr_p)/X_eq;

    T_e = X_ms * (i_qs * i_dr_p - i_ds * i_qr_p); %rotor ref. frame

    lambda_qs = (X_ls * i_qs + X_ms * (i_qs+i_qr_p))/omega_b;
    lambda_ds = (X_ls * i_ds + X_ms * (i_ds+i_dr_p))/omega_b;
    lambda_qr_p = (X_lr_p * i_qr_p + X_ms * (i_qs + i_qr_p))/omega_b;
    lambda_dr_p = (X_lr_p * i_dr_p + X_ms * (i_ds + i_dr_p))/omega_b;

    % Stationary reference frame
    omega = 0;
    %Synchronous ref frame
    %omega = omega_b;
        
    if ~ss
        dydt = [
            omega_b* (X_lr_p*v_qs + X_ms*v_qs - X_lr_p*i_qs*r_s - X_ms*i_qs*r_s + X_ms*i_qr_p*r_r_p - X_ms*lambda_dr_p*omega_r)/X_eq
            omega_b* (X_lr_p*v_ds + X_ms*v_ds - X_lr_p*i_ds*r_s - X_ms*i_ds*r_s + X_ms*i_dr_p*r_r_p + X_ms*lambda_qr_p*omega_r)/X_eq
            omega_b*-(X_ms*v_qs + X_ls*i_qr_p*r_r_p - X_ms*i_qs*r_s + X_ms*i_qr_p*r_r_p - X_ls*lambda_dr_p*omega_r - X_ms*lambda_dr_p*omega_r)/X_eq
            omega_b*-(X_ms*v_ds + X_ls*i_dr_p*r_r_p - X_ms*i_ds*r_s + X_ms*i_dr_p*r_r_p + X_ls*lambda_qr_p*omega_r + X_ms*lambda_qr_p*omega_r)/X_eq
            omega_b / (2*H) * (T_e - T_I)
            omega_r
        ];
    else
        dydt = [
            0
            omega_b* (X_lr_p*v_ds + X_ms*v_ds - X_lr_p*i_ds*r_s - X_ms*i_ds*r_s + X_ms*i_dr_p*r_r_p + X_ms*lambda_qr_p*omega_r)/X_eq
            omega_b*-(X_ls*i_qr_p*r_r_p + X_ms*i_qr_p*r_r_p - X_ls*lambda_dr_p*omega_r - X_ms*lambda_dr_p*omega_r)/(X_eq+X_ms^2) 
            omega_b*-(X_ms*v_ds + X_ls*i_dr_p*r_r_p - X_ms*i_ds*r_s + X_ms*i_dr_p*r_r_p + X_ls*lambda_qr_p*omega_r + X_ms*lambda_qr_p*omega_r)/X_eq
            omega_b / (2*H) * (T_e - T_I)
            omega_r
        ];
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