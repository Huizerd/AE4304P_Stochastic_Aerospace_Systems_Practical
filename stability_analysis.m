% Stability Analysis of Aircraft

clear; close all; clc

%% Load Data

ss_symmetric

%% Define Input Vectors

% Time vector
T  = 200.;
dt = 0.01;
t  = 0:dt:T;
N  = length(t);

% Elevator and turbulence inputs
d_e = [zeros(1,10) -0.05*ones(1,N-10)];  % input in radians
w_1 = zeros(1,N);                        % no turbulence
w_3 = zeros(1,N);                        % no turbulence

u  = [d_e' w_1' w_3'];

%% Simulation

y   = lsim(A, B, C, D, u, t);             % first with undamped matrix
y_t = lsim(A_t, B, C_t, D, u, t);         % then with pitch damper

%% Plotting Results

set(0, 'DefaultAxesTickLabelInterpreter','Latex')
set(0, 'DefaultLegendInterpreter','Latex')
set(0, 'DefaultFigurePosition', [152.5 168 719 791.5])

ax1 = subplot(5, 1, 1);
plot(t, y(:,1), '--', t, y_t(:,1))
ylabel('$\hat{u}$ [-]', 'Interpreter', 'Latex'); title('\textbf{Aircraft Response to Elevator Step}', 'Interpreter', 'Latex')
grid on
legend(ax1, 'Undamped', 'Pitch Damper')
legend('boxoff')
legend('Location', [0.80 0.892 0 0])

subplot(5, 1, 2)
plot(t, y(:,2)*180/pi, '--', t, y_t(:,2)*180/pi)
ylabel('$\alpha$ [deg]', 'Interpreter', 'Latex')
grid on

subplot(5, 1, 3)
plot(t, y(:,3)*180/pi, '--', t, y_t(:,3)*180/pi)
ylabel('$\theta$ [deg]', 'Interpreter', 'Latex')
grid on

subplot(5 ,1, 4)
plot(t, y(:,4)*180/pi, '--', t, y_t(:,4)*180/pi)
ylabel('$\frac{q\overline{c}}{V}$ [deg]', 'Interpreter', 'Latex')
grid on

subplot(5, 1, 5)
plot(t, y(:,5), '--', t, y_t(:,5))
xlabel('$t$ [s]', 'Interpreter', 'Latex'); ylabel('$n_z$ [-]', 'Interpreter', 'Latex')
grid on

set(gcf, 'Renderer', 'Painters')
savefig('C:\Users\jesse\Google Drive\CnS\AE4304P_Stochastic_Aerospace_Systems_Practical\figures\stability.fig')
print('-painters', '-depsc', 'C:\Users\jesse\Google Drive\CnS\AE4304P_Stochastic_Aerospace_Systems_Practical\figures\stability.eps')