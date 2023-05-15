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