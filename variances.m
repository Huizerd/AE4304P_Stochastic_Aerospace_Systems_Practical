% Variances of Aircraft States
% Jesse Hagenaars, 08-02-18

clear; close all; clc

% rng('default')  % to make results reproducible

%% Load Data

ss_symmetric

%% Variances using Analytical PSD from State-Space

% Frequency vector
N_ss = 300;
w_ss = logspace(-2, 2, N_ss);

% Select input row (horizontal or vertical turbulence)
if sigma_ug > 0 && sigma_wg == 0
    u_ss = 2;
    S_uu = sigma_ug^2;  % intensity W is equal to variance!
elseif sigma_ug == 0 && sigma_wg > 0
    u_ss = 3;
    S_uu = sigma_wg^2;  % intensity W is equal to variance!
end

% Compute PSD
S_xx = (bode(A_t, B, C_t, D, u_ss, w_ss)).^2 * S_uu;

% Numerical integration (crude) of PSD
var_ss = sum(diff(w_ss)' .* S_xx(1:end-1,:)) / pi;

%% Variances using Impulse Response Method

% Time vector
T_irm  = 200;
dt_irm = 0.01;
t_irm  = 0:dt_irm:T_irm;
N_irm  = length(t_irm);

% Zero input, initial condition from B-matrix
u_irm  = zeros(3,N_irm);
x0_irm = sqrt(S_uu) * B(:,u_ss);  % sqrt(S_uu) = sigma for turbulence case
                                  % not divide by sqrt(dt)? since no white?

% Simulation of impulse response
h_irm = lsim(A_t, B, C_t, D, u_irm, t_irm, x0_irm);

% Variances, integration as a fn. of time
var_irm = cumsum(dt_irm * h_irm.^2);

%% Variances from Time Signals using MATLAB var.m

% Time vector
T_ts  = 200.;
dt_ts = 0.01;
t_ts  = 0:dt_ts:T_ts;
N_ts  = length(t_ts);

% Elevator and turbulence inputs, multiple realizations
N_rl = 100;
d_e  = zeros(1,N_ts);
w_1  = sigma_ug * randn(N_rl,N_ts) / sqrt(dt_ts);  % sqrt(dt) because of lsim, amplitude linked to s.d.
w_3  = sigma_wg * randn(N_rl,N_ts) / sqrt(dt_ts);  % sqrt(dt) because of lsim, amplitude linked to s.d.

y_tsm   = zeros(N_ts, size(C_t,1), N_rl);
var_tsm = zeros(N_rl, size(C_t,1));

for i = 1 : N_rl
    
    u_ts  = [d_e' w_1(i,:)' w_3(i,:)'];

    % Simulation
    y_tsm(:,:,i) = lsim(A_t, B, C_t, D, u_ts, t_ts);
    
    % Compute variance
    var_tsm(i,:) = var(y_tsm(:,:,i));
    
end

% Average
var_ts = mean(var_tsm);

%% Plotting

set(0, 'defaultAxesTickLabelInterpreter','latex')
set(0, 'defaultLegendInterpreter','latex')

ax1 = subplot(5, 1, 1);
plot(t_irm, var_ss(1)*ones(1,N_irm), '--', t_ts, var_ts(1)*ones(1,N_ts), '-.', t_irm, var_irm(:,1))
ylabel('$\sigma^2_{\hat{u}}$ [rad$^2$]', 'Interpreter', 'Latex'); title('\textbf{Variances of Aircraft States due to Horizontal Turbulence}', 'Interpreter', 'Latex')
legend(ax1, 'Analytical', ['Using \texttt{var}, ' num2str(N_rl) ' Realizations'] , 'Impulse Response Method')
legend('boxoff')
grid on

subplot(5, 1, 2)
plot(t_irm, var_ss(2)*ones(1,N_irm), '--', t_ts, var_ts(2)*ones(1,N_ts), '-.', t_irm, var_irm(:,2))
ylabel('$\sigma^2_{\alpha}$ [rad$^2$]', 'Interpreter', 'Latex');
grid on

subplot(5, 1, 3)
plot(t_irm, var_ss(3)*ones(1,N_irm), '--', t_ts, var_ts(3)*ones(1,N_ts), '-.', t_irm, var_irm(:,3))
ylabel('$\sigma^2_{\theta}$ [rad$^2$]', 'Interpreter', 'Latex');
grid on

subplot(5, 1, 4)
plot(t_irm, var_ss(4)*ones(1,N_irm), '--', t_ts, var_ts(4)*ones(1,N_ts), '-.', t_irm, var_irm(:,4))
ylabel('$\sigma^2_{\frac{q\overline{c}}{V}}$ [rad$^2$]', 'Interpreter', 'Latex');
grid on

subplot(5, 1, 5)
plot(t_irm, var_ss(5)*ones(1,N_irm), '--', t_ts, var_ts(5)*ones(1,N_ts), '-.', t_irm, var_irm(:,5))
xlabel('$t$ [s]', 'Interpreter', 'Latex'); ylabel('$\sigma^2_{n_z}$ [rad$^2$]', 'Interpreter', 'Latex');
grid on
