clear all;
clear;


load('Data.mat');


% ----------------------------------------------------------------
%                           Matrix A
% ----------------------------------------------------------------


q_bar = 0.5 * r * U0^2;

Y_v = ((q_bar * S) / (m * U0)) * Cy_beta;
Y_p = ((q_bar * S * b ) / (2 * m * U0)) * Cy_p;
Y_r = ((q_bar * S * b ) / (2 * m * U0)) * Cy_r;

L_v = ((q_bar * S * b) / (Ixx * U0)) * Cl_beta;
L_p = ((q_bar * S * b^2) / (2*Ixx * U0)) * Cl_p;
L_r = ((q_bar * S * b^2) / (2*Ixx * U0)) * Cl_r;

N_v = ((q_bar * S * b) / (Izz * U0)) * Cn_beta;
N_p = ((q_bar * S * b^2) / (2*Izz * U0)) * Cn_p;
N_r = ((q_bar * S * b^2) / (2*Izz * U0)) * Cn_r;

A = [Y_v Y_p -(U0 - Y_r) -g*cos(deg2rad(theta_0));
     L_v L_p L_r 0;
     N_v N_p N_r 0;
     0 1 0 0];

% ----------------------------------------------------------------
%                           Matrix B
% ----------------------------------------------------------------

Y_dr = ((q_bar * S) /(m * U0)) * Cy_dr;
Y_da = ((q_bar * S) /(m * U0)) * Cy_da;

L_dr = ((q_bar * S * b) /(Ixx * U0)) * Cl_dr;
L_da = ((q_bar * S * b) /(Ixx * U0)) * Cl_da;

N_dr = ((q_bar * S * b) /(Iyy * U0)) * Cn_dr;
N_da = ((q_bar * S * b) /(Iyy * U0)) * Cn_da;

B = [Y_dr Y_da ; L_dr L_da ; N_dr N_da; 0 0];

% ----------------------------------------------------------------
%                    Caracteristic equation
% ----------------------------------------------------------------

syms lambda

eq_char = charpoly(A, lambda);
eq_char_simplified = vpa(eq_char, 4);

eigenvalues = eig(A);


% ----------------------------------------------------------------
%                 Different modes of lateral stability
% ----------------------------------------------------------------

% Spiral mode

lambda_spiral = eigenvalues(4);

% Roll mode

lambda_roll = eigenvalues(1);


% Dutch mode
omega_n_dutch = sqrt(eigenvalues(2)*eigenvalues(3));
zeta_dutch = (-real(eigenvalues(2))-real(eigenvalues(3)))/(2*omega_n_dutch);


% ----------------------------------------------------------------
%                   Curves of lateral motion
% ----------------------------------------------------------------

t = linspace(0, 20, 5000);

[V, D] = eig(A);

% need to write with the same forme as the one in the book (Lambda1 = roll / lambda2 = spiral / lambda34 = dutch)

V_spiral = V(:,4);
V_roll = V(:,1);
V_lambda3 = V(:,2);
V_lambda4 = V(:,3);

V = [V_roll, V_spiral, V_lambda3, V_lambda4];
lambda = [eigenvalues(1); eigenvalues(4); eigenvalues(2); eigenvalues(3)];


dv_roll = V(1,1) * exp(lambda(1)*t);
dp_roll = V(2,1) * exp(lambda(1)*t);
dr_roll = V(3,1) * exp(lambda(1)*t);
dphi_roll = V(4,1) * exp(lambda(1)*t);

dv_spiral = V(1,2) * exp(lambda(2)*t);
dp_spiral = V(2,2) * exp(lambda(2)*t);
dr_spiral = V(3,2) * exp(lambda(2)*t);
dphi_spiral = V(4,2) * exp(lambda(2)*t);

dv_dutch = V(1,3) * exp(lambda(3)*t) + V(1,4) * exp(lambda(4)*t);
dp_dutch = V(2,3) * exp(lambda(3)*t) + V(2,4) * exp(lambda(4)*t);
dr_dutch = V(3,3) * exp(lambda(3)*t) + V(3,4) * exp(lambda(4)*t);
dphi_dutch = V(4,3) * exp(lambda(3)*t) + V(4,4) * exp(lambda(4)*t);

figure;
subplot(2, 2, 1)
plot(t, dv_roll); 
xlabel('Time (s)');
ylabel('\Delta v');
title('Lateral Velocity');
grid on;

subplot(2, 2, 2)
plot(t, dp_roll); 
xlabel('Time (s)');
ylabel('\Delta p');
title('Roll rate');
grid on;

subplot(2, 2, 3)
plot(t, dr_roll); 
xlabel('Time (s)');
ylabel('\Delta r');
title('Yaw rate');
grid on;

subplot(2, 2, 4)
plot(t, dphi_roll); 
xlabel('\Delta \phi');
ylabel('caca');
title('Lateral angle');
grid on;

sgtitle('Lateral Stability - Roll Mode');

figure;
subplot(2, 2, 1)
plot(t, dv_spiral); 
xlabel('Time (s)');
ylabel('\Delta v');
title('Lateral Velocity');
grid on;

subplot(2, 2, 2)
plot(t, dp_spiral); 
xlabel('Time (s)');
ylabel('\Delta p');
title('Roll rate');
grid on;

subplot(2, 2, 3)
plot(t, dr_spiral); 
xlabel('Time (s)');
ylabel('\Delta r');
title('Yaw rate');
grid on;

subplot(2, 2, 4)
plot(t, dphi_spiral); 
xlabel('\Delta \phi');
ylabel('caca');
title('Lateral angle');
grid on;

sgtitle('Lateral Stability - Spirale Mode');

figure;
subplot(2, 2, 1)
plot(t, dv_dutch); 
xlabel('Time (s)');
ylabel('\Delta v');
title('Lateral Velocity');
grid on;

subplot(2, 2, 2)
plot(t, dp_dutch); 
xlabel('Time (s)');
ylabel('\Delta p');
title('Roll rate');
grid on;

subplot(2, 2, 3)
plot(t, dr_dutch); 
xlabel('Time (s)');
ylabel('\Delta r');
title('Yaw rate');
grid on;

subplot(2, 2, 4)
plot(t, dphi_dutch); 
xlabel('\Delta \phi');
ylabel('caca');
title('Lateral angle');
grid on;

sgtitle('Lateral Stability - Dutch Mode');

%----------------------------------------------------------------
%                       Transfer functions   
%----------------------------------------------------------------

syms s;

I = eye(size(A));

% Aileron response

R_aileron = inv(s*I-A) * B(:,2);


for i = 1:numel(R_aileron)
    [num, den] = numden(vpa(R_aileron(i),3));
    
    if i == 1
          d_vs_da = tf(sym2poly(num), sym2poly(den));
    end

    if i == 2
          d_ps_da = tf(sym2poly(num), sym2poly(den));
    end

    if i == 4
          d_rs_da = tf(sym2poly(num), sym2poly(den));
    end

    if i == 3
          d_phis_da = tf(sym2poly(num), sym2poly(den));
    end
 
end

figure;
bode(d_vs_da, d_ps_da, d_rs_da, d_phis_da);
legend("V(s) : Side velocity", " p(s) : Roll rate", " r(s) : Yaw rate", "\phi (s) :Side angle")
title('Bode diagram of the transfer function : (Response to Aileron)');

% Rudder response

R_rudder = inv(s*I-A) * B(:,1);

for i = 1:numel(R_aileron)
     [num, den] = numden(vpa(R_rudder(i),3));
     
     if i == 1
          d_vs_dr = tf(sym2poly(num), sym2poly(den));
     end
 
     if i == 2
          d_ps_dr = tf(sym2poly(num), sym2poly(den));
     end
 
     if i == 4
          d_rs_dr = tf(sym2poly(num), sym2poly(den));
     end
 
     if i == 3
          d_phis_dr = tf(sym2poly(num), sym2poly(den));
     end
  
end

figure;
bode(d_vs_dr, d_ps_dr, d_rs_dr, d_phis_dr);
legend("V(s) : Side velocity", " p(s) : Roll rate", " r(s) : Yaw rate", "\phi (s) :Side angle")
title('Bode diagram of the transfer function : (Response to Rudder)');