clc;
clear;
close all;

% Paramètres de simulation
animation = false;
   % peut être "reference", "cylinder", "spiral" ou "custom"

% Paramètres physiques
g = 10;
m = 1.63;
I_xx = 0.0151;
I_yy = 0.0092;
I_zz = 0.0093;
a1 = (I_yy - I_zz) / I_xx;
a2 = (I_zz - I_xx) / I_yy;
a3 = (I_xx - I_yy) / I_zz;
l = 1;
b1 = l / I_xx;
b2 = 1 / I_yy;
b3 = 1 / I_zz;
delta_t = 0.05;
sim_time = 550;

% Gains de contrôle
k1 = 10; k2 = 0.5;
k3 = 2; k4 = 10;
k5 = 12; k6 = 10;
k7 = 10; k8 = 0.5;
k9 = 10; k10 = 15;
k11 = 8; k12 = 10;
noise_level = 0;
trajectory = "spirale"; 
% Conditions initiales
x = [0,0,0,0,0,0,-5,0,-5,0,20,0]; % Initial state [phi,phi_dot,theta,theta_dot, psi,psi_dot, z,z_dot, x,x_dot, y , y_dot]
if trajectory == "reference"
% Définition des fonctions pour une trajectoire circulaire
% Paramètres de la trajectoire spirale
r_0 = 2;       % Rayon initial
alpha = 0.1;   % Taux d'augmentation du rayon
omega = 0.2;   % Vitesse angulaire
z_0 = 0;       % Altitude initiale
v_z = 0.1;     % Vitesse verticale

% Paramètres de la trajectoire droite
x_0 = 50; y_0 = 55; z_0 = 20; % Positions initiales
v_x = 0.5; v_y = 0.5; v_z = 0.2; % Vitesses constantes

% Définition des fonctions pour une trajectoire droite
xd_fun = @(t) x_0 + v_x * t;
yd_fun = @(t) y_0 + v_y * t;
zd_fun = @(t) z_0 + v_z * t;

xd_dotfun = @(t) v_x;
yd_dotfun = @(t) v_y;
zd_dotfun = @(t) v_z;

xd_2dotfun = @(t) 0;
yd_2dotfun = @(t) 0;
zd_2dotfun = @(t) 0;

psi_d_fun = @(t) 0+t*0;
    psi_d_dotfun = @(t) 0+t*0;
    psi_d_2dotfun = @(t) 0+t*0;
    
    phi_d_fun = @(t) 0+t*0;
    phi_d_dotfun = @(t) 0+t*0;
    phi_d_2dotfun = @(t) 0+t*0;
    
    theta_d_fun = @(t) 0+t*0;
    theta_d_dotfun = @(t) 0+t*0;
    theta_d_2dotfun = @(t) 0+t*0;
end
if trajectory == "circulaire"
    % Définition des fonctions pour une trajectoire circulaire
% Paramètres de la trajectoire circulaire
r = 10; % Rayon du cercle
omega = 0.2; % Vitesse angulaire
z_constant = 4; % Altitude constante
xd_fun = @(t) r * cos(omega * t);
yd_fun = @(t) r * sin(omega * t);
zd_fun = @(t) z_constant;

xd_dotfun = @(t) -r * omega * sin(omega * t);
yd_dotfun = @(t) r * omega * cos(omega * t);
zd_dotfun = @(t) 0;

xd_2dotfun = @(t) -r * omega^2 * cos(omega * t);
yd_2dotfun = @(t) -r * omega^2 * sin(omega * t);
zd_2dotfun = @(t) 0;

psi_d_fun = @(t) 0+t*0;
    psi_d_dotfun = @(t) 0+t*0;
    psi_d_2dotfun = @(t) 0+t*0;
    
    phi_d_fun = @(t) 0+t*0;
    phi_d_dotfun = @(t) 0+t*0;
    phi_d_2dotfun = @(t) 0+t*0;
    
    theta_d_fun = @(t) 0+t*0;
    theta_d_dotfun = @(t) 0+t*0;
    theta_d_2dotfun = @(t) 0+t*0;
end
if trajectory == "spirale"
% Paramètres de la trajectoire spirale
r_0 = 2;       % Rayon initial
alpha = 0.1;   % Taux d'augmentation du rayon
omega = 0.2;   % Vitesse angulaire
z_0 = 0;       % Altitude initiale
v_z = 0.1;     % Vitesse verticale

% Définition des fonctions pour une trajectoire spirale
r_fun = @(t) r_0 + alpha * t;
xd_fun = @(t) r_fun(t) * cos(omega * t);
yd_fun = @(t) r_fun(t) * sin(omega * t);
zd_fun = @(t) z_0 + v_z * t;

xd_dotfun = @(t) -r_fun(t) * omega * sin(omega * t) + alpha * cos(omega * t);
yd_dotfun = @(t) r_fun(t) * omega * cos(omega * t) + alpha * sin(omega * t);
zd_dotfun = @(t) v_z;

xd_2dotfun = @(t) -alpha * omega * sin(omega * t) - r_fun(t) * omega^2 * cos(omega * t);
yd_2dotfun = @(t) alpha * omega * cos(omega * t) - r_fun(t) * omega^2 * sin(omega * t);
zd_2dotfun = @(t) 0;

psi_d_fun = @(t) 0+t*0;
    psi_d_dotfun = @(t) 0+t*0;
    psi_d_2dotfun = @(t) 0+t*0;
    
    phi_d_fun = @(t) 0+t*0;
    phi_d_dotfun = @(t) 0+t*0;
    phi_d_2dotfun = @(t) 0+t*0;
    
    theta_d_fun = @(t) 0+t*0;
    theta_d_dotfun = @(t) 0+t*0;
    theta_d_2dotfun = @(t) 0+t*0;
end

   
    % Trajectoire désirée et simulée
desired_trajectory = [];
simulated_trajectory = [];

% Simulation
t_interval = 0:delta_t:sim_time;
for k = 1:length(t_interval)
    t = t_interval(k);

    % Trajectoire de référence
    xd = xd_fun(t);
    yd = yd_fun(t);
    zd = zd_fun(t);

    xd_dot = xd_dotfun(t);
    yd_dot = yd_dotfun(t);
    zd_dot = zd_dotfun(t);

    xd_2dot = xd_2dotfun(t);
    yd_2dot = yd_2dotfun(t);
    zd_2dot = zd_2dotfun(t);

    phi_d = phi_d_fun(t);
    phi_d_dot = phi_d_dotfun(t);
    phi_d_2dot = phi_d_2dotfun(t);

    theta_d = theta_d_fun(t);
    theta_d_dot = theta_d_dotfun(t);
    theta_d_2dot = theta_d_2dotfun(t);


    psi_d = psi_d_fun(t);
    psi_d_dot = psi_d_dotfun(t);
    psi_d_2dot = psi_d_2dotfun(t);
    
    desired_trajectory = [desired_trajectory; xd, yd, zd];

    %definintion des erreur
    e1 = phi_d - x(1);
    e1_dot = phi_d_dot -x(2);
    e2 = x(2) - phi_d_dot -k1*e1;


    e3 = theta_d - x(3);
    e3_dot = theta_d_dot - x(2);
    e4 = x(4) - theta_d_dot -k3*e3;

    e5 = psi_d - x(5);
    e5_dot = psi_d_dot - x(6);
    e6 = x(6) - psi_d_dot -k4*e4;

    e5 = psi_d - x(5);
    e5_dot = psi_d_dot - x(6);
    e6 = x(6) - psi_d_dot -k5*e5;

    e7 = zd - x(7);
    e7_dot = zd_dot - x(8);
    e8 = x(8) - zd_dot -k7*e7;

    e9 = xd - x(9);
    e9_dot = xd_dot - x(10);
    e10 = x(10) - xd_dot -k9*e9;

    e11 = yd - x(11);
    e11_dot = yd_dot - x(12);
    e12 = x(12) - yd_dot -k11*e11;

    % Lois de contrôle simplifiées
    U1 = (m/cos(x(1))*cos(x(3)))*(g+zd_2dot+k7*e7_dot+e7-k8*e8);
    U2 = (1/b1)*(x(4)*x(6)*a1 +phi_d_2dot+k1*e1_dot-e1-k2*e2);
    U3 = (1/b2)*(x(2)*x(6)*a2 +theta_d_2dot+k3*e3_dot-e3-k4*e4);
    U4 = (1/b3)*(x(2)*x(4)*a3 +psi_d_2dot+k5*e5_dot-e5-k6*e6);
    Ux = (m/U1)*(xd_2dot+k9*e9_dot+e9-k10*e10);
    Uy = (m/U1)*(yd_2dot+k11*e11_dot+e11-k12*e12);

    % Mise à jour des états (Euler)
    noise = noise_level * randn(1, 12);  % Generate Gaussian noise
    x(1:12) = x(1:12) + delta_t * [x(2), U2, x(4), U3, x(6), U1, x(8), (U1 - m * g), x(10), Ux, x(12), Uy] + noise;
    simulated_trajectory = [simulated_trajectory; x(9), x(11), x(7)];
end



figure;
plot3(desired_trajectory(:, 1), desired_trajectory(:, 2), desired_trajectory(:, 3), 'b--', 'LineWidth', 1.5); hold on;
plot3(simulated_trajectory(:, 1), simulated_trajectory(:, 2), simulated_trajectory(:, 3), 'r-', 'LineWidth', 1.5);
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Trajectoire : Désirée vs Simulée');
legend('Trajectoire désirée', 'Trajectoire simulée');
grid on;

