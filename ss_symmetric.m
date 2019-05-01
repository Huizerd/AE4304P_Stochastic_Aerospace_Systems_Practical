% Computation of the Symmetric State Space

clear; close all; clc

%% Aircraft Specifications & Flight Conditions

x_cg = 0.30;     % c
W    = 46042.;   % N
m    = 4695.;    % kg
S    = 24.2;     % m^2
c    = 2.022;    % m
b    = 13.36;    % m

V    = 126.4;    % m/s
h    = 6096.;    % m
rho  = 0.6528;   % kg/m^3
g    = 9.80665;  % m/s^2

mu_c = 147.;
mu_b = 22.;
K_Y  = sqrt(0.950);

%% Turbulence Parameters

disp(' ');
sigma_ug = input('   Enter hor. turb. intensity sigma_ug  (2) [m/s] : ');
sigma_wg = input('   Enter vert. turb. intensity sigma_wg (3) [m/s] : ');
L_g      = input('   Enter turb. scale length L_g         (150) [m] : ');

sigma_ugV = sigma_ug / V;
sigma_ag  = sigma_wg / V;

%% Longitudinal Stability & Control Derivatives

% Needed for computation of some derivatives, not given
C_mac    =  0.;
C_mh     =  0.;
l_h      =  0.;

C_X0     =  0.;
C_Xu     = -0.0638;
C_Xa     =  0.1382;
C_Xq     =  0.;
C_Xde    =  0.;
C_Xug    =  C_Xu;
C_Xudotg =  0.;
C_Xag    =  C_Xa;
C_Xadotg =  0.;

C_Z0     = -0.3650;
C_Zu     = -0.7300;
C_Za     = -5.4300;
C_Zadot  = -1.6200;
C_Zq     = -4.0700;
C_Zde    = -0.5798;
C_Zug    =  C_Zu;
C_Zudotg =  C_mac * 2.;
C_Zag    =  C_Za;
C_Zadotg =  C_Zadot - C_Zq;

C_mu     =  0.;
C_ma     = -0.5180;
C_madot  = -4.2000;
C_mq     = -7.3500;
C_mde    = -1.4440;
C_mug    =  C_mu;
C_mudotg =  C_mh * -2.*l_h/c;
C_mag    =  C_ma;
C_madotg =  C_madot - C_mq;

%% Symmetric State-Space Entries

xu     = V/c * C_Xu/(2*mu_c);
xa     = V/c * C_Xa/(2*mu_c);
xt     = V/c * C_Z0/(2*mu_c);
xq     = V/c * C_Xq/(2*mu_c);
xd     = V/c * C_Xde/(2*mu_c);
xug    = V/c * C_Xug/(2*mu_c);
xudotg = V/c * C_Xudotg/(2*mu_c);
xag    = V/c * C_Xag/(2*mu_c);
xadotg = V/c * C_Xadotg/(2*mu_c);

zu     = V/c *  C_Zu/(2*mu_c - C_Zadot);
za     = V/c *  C_Za/(2*mu_c - C_Zadot);
zt     = V/c * -C_X0/(2*mu_c - C_Zadot);
zq     = V/c * (2*mu_c + C_Zq)/(2*mu_c - C_Zadot);
zd     = V/c *  C_Zde/(2*mu_c - C_Zadot);
zug    = V/c *  C_Zug/(2*mu_c - C_Zadot);
zudotg = V/c *  C_Zudotg/(2*mu_c - C_Zadot);
zag    = V/c *  C_Zag/(2*mu_c - C_Zadot);
zadotg = V/c *  C_Zadotg/(2*mu_c - C_Zadot);

mu     = V/c *  (C_mu + C_Zu*C_madot/(2*mu_c - C_Zadot))/(2*mu_c*K_Y^2);
ma     = V/c *  (C_ma + C_Za*C_madot/(2*mu_c - C_Zadot))/(2*mu_c*K_Y^2);
mt     = V/c * -(C_X0*C_madot/(2*mu_c - C_Zadot))/(2*mu_c*K_Y^2);
mq     = V/c *  (C_mq + C_madot*(2*mu_c + C_Zq)/(2*mu_c - C_Zadot))/(2*mu_c*K_Y^2);
md     = V/c *  (C_mde + C_Zde*C_madot/(2*mu_c - C_Zadot))/(2*mu_c*K_Y^2);
mug    = V/c *  (C_mug + C_Zug*C_madot/(2*mu_c - C_Zadot))/(2*mu_c*K_Y^2);
mudotg = V/c *  (C_mudotg + C_Zudotg*C_madot/(2*mu_c - C_Zadot))/(2*mu_c*K_Y^2);
mag    = V/c *  (C_mag + C_Zag*C_madot/(2*mu_c - C_Zadot))/(2*mu_c*K_Y^2);
madotg = V/c *  (C_madotg + C_Zadotg*C_madot/(2*mu_c - C_Zadot))/(2*mu_c*K_Y^2);

%% Symmetric State-Space

A = [ xu  xa  xt    0                       xug         xag             0;
      zu  za  zt   zq  zug - zudotg*V/L_g*(c/V)         zag  zadotg*(c/V);
       0   0   0  V/c                         0           0             0;
      mu  ma  mt   mq  mug - mudotg*V/L_g*(c/V)         mag  madotg*(c/V);
       0   0   0    0                    -V/L_g           0             0;
       0   0   0    0                         0           0             1;
       0   0   0    0                         0  -(V/L_g)^2    -2*V/L_g ];
  
B = [ xd                                     0                                         0;
      zd  zudotg*(c/V)*sigma_ugV*sqrt(2*V/L_g)       zadotg*(c/V)*sigma_ag*sqrt(3*V/L_g);
       0                                     0                                         0;
      md  mudotg*(c/V)*sigma_ugV*sqrt(2*V/L_g)       madotg*(c/V)*sigma_ag*sqrt(3*V/L_g);
       0               sigma_ugV*sqrt(2*V/L_g)                                         0;
       0                                     0                    sigma_ag*sqrt(3*V/L_g);
       0                                     0  (1-2*sqrt(3))*sigma_ag*sqrt((V/L_g)^3) ];
   
%% Pitch Damper & Matrix At

K_t = -0.05;
K_q = -3.;
K   = [0 0 K_t K_q 0 0 0];  % feedback matrix
A_t = A - B(:,1)*K;         % new A matrix = (A - BK), due to feedback

%% C & D Matrices

% Computation of load factor for damped and undamped state matrix      
C_n_t = V/g * (A_t(3,:) - A_t(2,:));  % computation of load factor, defined as
C_n   = V/g * (A(3,:) - A(2,:));      % V/g * (thetadot - alphadot) (see Ex. 7.3)

C     = [eye(4) zeros(4,3); C_n];     % only a/c states as output
C_t   = [eye(4) zeros(4,3); C_n_t];                   

% Feedthrough matrix doesn't change with damping
D_n = V/g * (B(3,:) - B(2,:));
D   = [zeros(4,3); D_n];   

%% Clean Up / Save

% save(filename, 'A', 'At', 'B', 'C', 'D', 'sigma_ugV', 'sigma_ag', 'L_g', 'V', 'c', 'g')
%
% clearvars -except filename
%
% load(filename)

clearvars -except A A_t B C C_t D sigma_ug sigma_wg L_g V c g
