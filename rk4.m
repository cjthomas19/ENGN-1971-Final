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