clear all;
clear;


load('Data.mat');


% ----------------------------------------------------------------
%                           Aircraft parameters
% ----------------------------------------------------------------

q_bar = 0.5 * r * U0^2;

X_u = -(q_bar * S ) / (m * U0) * (2 * CD_0 + CD_u);
X_w = ((q_bar * S) / (m * U0)) * (CL_0 * (1 - (2* CL_alpha / (pi * e * A))));


Z_u = -(q_bar * S / (m * U0)) * (2 * CL_0 + CL_u);
Z_w = -(q_bar * S / (m * U0)) * (CD_0 + CL_alpha); 
Z_w_dot = (q_bar * S * c_bar / (2 * m * U0^2)) * (CD_0 + CL_alpha_dot);
Z_q = (q_bar * S * c_bar / (2 * m * U0)) * (CL_q);

M_u = ((q_bar * S * c_bar) / (Iyy* U0)) * Cm_u;
M_w = ((q_bar * S * c_bar) / (Iyy* U0)) * Cm_alpha;
M_w_dot = ((q_bar * S * c_bar^2) / (2*Iyy* U0^2)) * Cm_alpha_dot;
M_q = ((q_bar * S * c_bar^2) / (2*Iyy* U0)) * Cm_q;


% ----------------------------------------------------------------
%                           Matrix A
% ----------------------------------------------------------------

A11 = X_u;
A12 = X_w;
A13 = 0;
A14 = -g*cos(deg2rad(theta_0));

A21 = Z_u;
A22 = Z_w;
A23 = U0;
A24 = -g*sin(deg2rad(theta_0));

A31 = M_u + Z_u*M_w_dot;
A32 = M_w + Z_w*M_w_dot;
A33 = M_q + U0*M_w_dot;
A34 = 0;

A41 = 0;
A42 = 0;
A43 = 1;
A44 = 0;

A = [A11 A12 A13 A14; A21 A22 A23 A24; A31 A32 A33 A34; A41 A42 A43 A44];


% ----------------------------------------------------------------
%                           Matrix B
% ----------------------------------------------------------------

X_delta_T = ((q_bar * S) / (m * U0)) * CD_delta_T;
X_delta_e = ((q_bar * S) / (m * U0)) * CD_delta_e;

Z_delta_T = ((q_bar * S) / (m * U0)) * CL_delta_T;
Z_delta_e = ((q_bar * S) / (m * U0)) * CL_delta_e;

M_delta_T = ((q_bar * S * c_bar) / (Iyy* U0)) * Cm_delta_T;
M_delta_e = ((q_bar * S * c_bar) / (Iyy* U0)) * Cm_delta_e;

B11 = X_delta_e;
B12 = X_delta_T;

B21 = Z_delta_e;
B22 = Z_delta_T;

B31 = M_delta_e + Z_delta_e*M_w_dot;
B32 = M_delta_T + Z_delta_T*M_w_dot;

B41 = 0;
B42 = 0;

B = [B11 B12; B21 B22; B31 B32; B41 B42];

% ----------------------------------------------------------------
%                    Caracteristic equation
% ----------------------------------------------------------------

syms lambda

eq_char = charpoly(A, lambda);
eq_char_simplified = vpa(eq_char, 4);

eigenvalues = eig(A);


% ----------------------------------------------------------------
%                  Modes of longitudinal stability
% ----------------------------------------------------------------


eq1 = (lambda - eigenvalues(1))*(lambda - eigenvalues(2));

eq2 = (lambda - eigenvalues(3))*(lambda - eigenvalues(4));

expanded_expr1 = vpa(expand(eq1),3);
expanded_expr2 = vpa(expand(eq2),3);

% Short period
omega_n_sp = sqrt(eigenvalues(1) * eigenvalues(2));
zeta_sp = (-real(eigenvalues(1))-real(eigenvalues(2)))/(2*omega_n_sp);

% Phugoid
omega_n_p = sqrt(eigenvalues(3) * eigenvalues(4));
zeta_p = (-real(eigenvalues(3))-real(eigenvalues(4)))/(2*omega_n_p);

% ----------------------------------------------------------------
%                Curves of longitudinal motion
% ----------------------------------------------------------------

%syms t;
t = linspace(0, 800, 5000);

[V, D] = eig(A); % V : matrix of eigenvectors / D:diagonal matrix of eigenvalues

V13 = real(V(1,3));
V14 = real(V(1,4));
V43 = real(V(4,3));
V44 = real(V(4,4));
V21 = real(V(2,1));
V22 = real(V(2,2));
V31 = real(V(3,1));
V32 = real(V(3,2));

du = V13 * exp(eigenvalues(3)*t) + V14 * exp(eigenvalues(4)*t);
d_theta = V43 * exp(eigenvalues(3)*t) + V44 * exp(eigenvalues(4)*t);

dw = V21 * exp(eigenvalues(1)*t) + V22 * exp(eigenvalues(2)*t);
dq = V31 * exp(eigenvalues(1)*t) + V32 * exp(eigenvalues(2)*t);

U = U0 + du;
alpha = atan(dw ./ U);

subplot(2,2,1);
plot(t, U); 
xlabel('Time (s)');
ylabel('Axial velocity U(m/s)');
title('U(t)');
grid on;


subplot(2,2,2);
plot(t, rad2deg(alpha)); 
xlabel('Time (s)');
ylabel('Angle of Attack (\circ)');
title('\alpha');
xlim([0, 5]);
grid on;

subplot(2,2,3);
plot(t, dq); 
xlabel('Time (s)');
ylabel('\delta q)');
title('Pitch Rate');
xlim([-1, 5]);
grid on;


subplot(2,2,4);
plot(t, d_theta); 
xlabel('Time (s)');
ylabel('\theta');
title('Pitch Angle');
grid on;



%----------------------------------------------------------------
%                       Transfer functions   
%----------------------------------------------------------------


syms s; 

I = eye(size(A));

D = adjoint(s*I-A)/det(s*I- A) * B;


for i = 1:numel(D)
    [num, den] = numden(vpa(D(i),3));
    
    if i == 1
        D_us = tf(sym2poly(num), sym2poly(den));
    end

    if i == 2
        D_ws = tf(sym2poly(num), sym2poly(den));
    end

    if i == 4
        D_thetas = tf(sym2poly(num), sym2poly(den));
    end

    if i == 3
        D_qs = tf(sym2poly(num), sym2poly(den));
    end
 
end

figure;

bode(D_us, D_ws,D_qs, D_thetas)
legend('U(s) : velocity axial ', 'w(s) : angle of attack', 'q(s) : Pitch rate', '\theta(s) : Pitch angle');
title('Bode diagram of the transfer function : Response to Elevator');

