clear;
omega_b = 60*2*pi;
r_r_p = 1.077227818558939;
L_r_p = 0.048933010422115;
X_r_p = L_r_p * omega_b;

syms Xprm Rcore Xsl r_s real positive

Zs1 = (1/Rcore + 1/(1j*Xprm) + 1/(r_r_p+1j*X_r_p))^-1 + r_s + 1j*Xsl;
Zs0 = (1/Rcore + 1/(1j*Xprm) + 1/(r_r_p+2j*X_r_p))^-1 + r_s + 1j*Xsl;

Zs1 = ((r_s + 1j*Xsl + (r_r_p+1j*X_r_p))^-1 + 1/Rcore - 1j/Xprm)^-1;
Zs0 = ((r_s + 1j*Xsl + (r_r_p+2j*X_r_p))^-1 + 1/Rcore - 1j/Xprm)^-1;


eqns = [
    real(Zs1) == 4.444
    imag(Zs1) == 4.433
    real(Zs0) == 5.94
    imag(Zs0) == 21.16
];
soln = solve(eqns,[Rcore Xprm r_s Xsl])
