clear;
close all;
data = readmatrix('data.csv');

rshunt = 0.514; % Measured shunt resistance
mapp = data(:,2)/1000;
V = data(:,4);
I = data(:,5)/rshunt;

g = 9.81;
R = 3.83 * 2.54 / 100;

T = mapp * g * R;


rpm = data(:,3)/8 * 60;
w = data(:,3)/8 * 2 * pi;

dt = data(:,6)/1000*60*2*pi * 180/pi;
pf = cosd(dt);

P = T .* w/746;
x_r = find(abs(mapp - 1) == min(abs(mapp-1)));
eff = (T.*w)./(V.*I.*pf);
slip = (1800-rpm)/1800;

%% Plots
figure;
hold on;
plot(rpm,P,'.-','markersize',20,'linewidth',2);
plot(rpm,V.*I.*pf/746,'.-','markersize',20,'linewidth',2)
plot(rpm,V.*I.*sqrt(1-pf.^2)/746,'.-','markersize',20,'linewidth',2)
plot([rpm(x_r),rpm(x_r)],[0,1],'k--','linewidth',1)
grid on;
format_plot("RPM","Power (hp)");
legend(["Mechanical","Electrical (Real)","Electrical (Imag)"])

figure;
hold on;
plot(P,pf,'linewidth',2)
plot(P,eff,'linewidth',2)
plot([.25,.25],[0,1],'k--','linewidth',1)
format_plot("Mechanical Power (hp)","Percentage (\%)");
legend(["Power Factor","Efficiency"])
grid on;

figure;
plot(rpm,T,'linewidth',2);
format_plot("RPM","Torque (N-m)");
grid on;



%% Fitting parameters
wl = 60*2*pi/2;
testP = @(b,x) 120.^2 ./(1+(b(1).*x.*(wl)).^2) .*(x.*(1-x))/(2.*b(2));
testT = @(b,x) 120.^2 ./(1+(b(1).*x.*(wl)).^2) .*x./(2.*b(2).*wl);

mdl = fitnlm(slip,T,testT,[1,1]);
S = linspace(0,max(slip),100);

format long;
disp("Estimate:");
disp("tau_r = " + mdl.Coefficients.Estimate(1));
disp("R'_r = " + mdl.Coefficients.Estimate(2));
disp("L'_r = " + mdl.Coefficients.Estimate(1)*mdl.Coefficients.Estimate(2));

%% Simulations

P_B = 0.25 * 745.7; % Convert HP to W
V_r = 120; % Rated voltage, Vrms

V_B = V_r*sqrt(2);
I_B = P_B/V_B;

Poles = 4; % Number of poles
omega_b = 60 * 2 * pi;

T_B = P_B / (2/Poles * omega_b);

J = 1e-2; % Inertia constant, J*s^2 = kg*m^2

Z_B = V_B / I_B;

r_s = 2.5/Z_B; % Stator winding resistance (measured)
r_r_p = (4.444-r_s*Z_B)/Z_B;
%r_r_p = mdl.Coefficients.Estimate(2)/Z_B;
%X_lr_p = omega_b*mdl.Coefficients.Estimate(2)*mdl.Coefficients.Estimate(1)/Z_B;
%X_ls = 4.61/Z_B; %Re-factor if needed
%X_ms = 66.42/Z_B; %Re-factor if needed

X_ls = 4.433/2/Z_B;
X_lr_p = 4.433/2/Z_B;
X_ms = (21.16-X_ls*Z_B)/Z_B;

H = 0.5 * (2/Poles) * J * omega_b / T_B;

v_in = 1;
%%
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

[st,sy] = rk4(@(t,y) simulate(t,y,param,[],{}), [0,10],[0,0,0,0,omega_b,0],1e-4);

y0 = sy(:,end);

torques = 0:0.1:2.5;
slips = zeros(size(torques));

for T_i = 1:numel(torques)
    disp(T_i/numel(torques) * 100 + "%")
    param(8) = torques(T_i);
    param(10) = 1;
    [st,sy] = rk4(@(t,y) simulate(t,y,param,[],{}),[0,10],y0,1e-4);
    
    % Extract solutions
    i_qs = sy(1,:);
    i_ds = sy(2,:);
    i_qr_p = sy(3,:);
    i_dr_p = sy(4,:);
    omega_r = sy(5,:);
    theta_r = sy(6,:);

    T_e = X_ms * (i_dr_p .* i_qs - i_qr_p .* i_ds); %Torque eqn in arbitrary ref. frame

    slipi = (omega_b - omega_r)/omega_b;

    if slipi(end) < 1
        slips(T_i) = slipi(end);
        y0 = sy(:,end);
    else
        slips(T_i) = NaN;
    end
end

%% Plot comparisons
figure;
hold on;
plot(slip,P,'.-','linewidth',2,'markersize',20);
format_plot("Slip","Mechanical Power (hp)")
plot(S,testP(mdl.Coefficients.Estimate,S)/746,'linewidth',2);
plot(slips, torques.*T_B.*(1-slips).*omega_b/746/2,'linewidth',2);
grid on;
legend(["Measured","Circuit Model","Numerical Model"])
%%
figure;
hold on;
plot(slip,T,'.-','linewidth',2,'markersize',20);
format_plot("Slip","Torque (N$\cdot$m)");
plot(S,testT(mdl.Coefficients.Estimate,S),'linewidth',2);
plot(slips,torques*T_B,'linewidth',2);
legend(["Measured","Circuit Model","Numerical Model"])
grid on;