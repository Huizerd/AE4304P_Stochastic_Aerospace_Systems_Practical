% Time Domain Analysis of Aircraft in Turbulence
% Jesse Hagenaars, 07-02-18

clear; close all; clc

% rng('default')  % to make results reproducible

%% Load Data

ss_symmetric

% Detect type of turbulence for plotting/saving
if sigma_ug > 0 && sigma_wg == 0
    filet = '\textbf{Aircraft Response to Horizontal Turbulence}';
    filen = 'C:\Users\jesse\Google Drive\CnS\AE4304P_Stochastic_Aerospace_Systems_Practical\figures\timeresponse_hor.eps';
elseif sigma_ug == 0 && sigma_wg > 0
    filet = '\textbf{Aircraft Response to Vertical Turbulence}';
    filen = 'C:\Users\jesse\Google Drive\CnS\AE4304P_Stochastic_Aerospace_Systems_Practical\figures\timeresponse_vert.eps';
end
%% Define Input Vectors

% Time vector
T  = 200.;
dt = 0.01;
t  = 0:dt:T;
N  = length(t);

% Elevator and turbulence inputs
d_e = zeros(1,N);
w_1 = sigma_ug * randn(1,N) / sqrt(dt);      % sqrt(dt) because of lsim, amplitude linked to s.d.
w_3 = sigma_wg * randn(1,N) / sqrt(dt);      % sqrt(dt) because of lsim, amplitude linked to s.d.

u  = [d_e' w_1' w_3'];                       % both horizontal and vertical turbulence
                                             % B-matrix multiplies with 0 if either w1 or
                                             % w3 not needed

%% Simulation

y = lsim(A_t, B, C_t, D, u, t);              % make use of the damped A- and C-matrices

%% Plotting Results

set(0, 'DefaultAxesTickLabelInterpreter','latex')
set(0, 'DefaultLegendInterpreter','latex')
set(0, 'DefaultFigurePosition', [152.5 168 719 791.5])

subplot(5, 1, 1)
plot(t, y(:,1))
ylabel('$\hat{u}$ [-]', 'Interpreter', 'Latex'); title(filet, 'Interpreter', 'Latex')
grid on

subplot(5, 1, 2)
plot(t, y(:,2)*180/pi)
ylabel('$\alpha$ [deg]', 'Interpreter', 'Latex')
grid on

subplot(5, 1, 3)
plot(t, y(:,3)*180/pi)
ylabel('$\theta$ [deg]', 'Interpreter', 'Latex')
grid on

subplot(5 ,1, 4)
plot(t, y(:,4)*180/pi)
ylabel('$\frac{q\overline{c}}{V}$ [deg]', 'Interpreter', 'Latex')
grid on

subplot(5, 1, 5)
plot(t, y(:,5))
xlabel('$t$ [s]', 'Interpreter', 'Latex'); ylabel('$n_z$ [-]', 'Interpreter', 'Latex')
grid on

set(gcf, 'Renderer', 'Painters')
print('-painters', '-depsc', filen)