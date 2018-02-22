% Spectral Analysis of Aircraft in Turbulence
% Jesse Hagenaars, 08-02-18

clear; close all; clc

% rng('default')  % to make results reproducible

%% Load Data

ss_symmetric

%% Analytical PSD from State-Space

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

%% Experimental PSD from FFT

% Time vector
T_fft  = 200.;
dt_fft = 0.01;
t_fft  = 0:dt_fft:T_fft;
N_fft  = length(t_fft);
fs_fft = 1 / dt_fft;

% Elevator and turbulence inputs, multiple realizations
N_rl = 100;
d_e  = zeros(1,N_fft);
w_1  = sigma_ug * randn(N_rl,N_fft) / sqrt(dt_fft);  % sqrt(dt) because of lsim, amplitude linked to s.d.
w_3  = sigma_wg * randn(N_rl,N_fft) / sqrt(dt_fft);  % sqrt(dt) because of lsim, amplitude linked to s.d.

y_fftm = zeros(N_fft, size(C_t,1), N_rl);
Y_fftm = zeros(N_fft, size(C_t,1), N_rl);
I_fftm = zeros(N_fft, size(C_t,1), N_rl);

for i = 1 : N_rl
    
    u_fft  = [d_e' w_1(i,:)' w_3(i,:)'];

    % Simulation
    y_fftm(:,:,i) = lsim(A_t, B, C_t, D, u_fft, t_fft);

    % FFT
    Y_fftm(:,:,i) = dt_fft * fft(y_fftm(:,:,i));  % DT transform of CT system!

    % Compute PSD (periodogram)
    I_fftm(:,:,i) = real((1/T_fft) * Y_fftm(:,:,i) .* conj(Y_fftm(:,:,i)));
    
end

% Average, take 1 realization for smoothing filter
I_fft  = mean(I_fftm, 3);
I_fft1 = I_fftm(:,:,1);

% Frequency vector for plotting
% Only plot 0 to pi
w_fft = 2*pi*fs_fft * (0:(N_fft/2) - 1) / N_fft;

%% Experimental PSD with Filter

% For smoothing, vector will have to be 2 elements shorter to allow shift
N_filt = N_fft - 2;
w_filt = 2*pi*fs_fft * (0:(N_filt/2) - 1) / N_filt;
I_filt = 0.25 * I_fft1(1:end-2,:) + 0.5 * I_fft1(2:end-1,:) + 0.25 * I_fft1(3:end,:);

%% Plotting

set(0, 'defaultAxesTickLabelInterpreter','latex')
set(0, 'defaultLegendInterpreter','latex')

ax1 = subplot(5, 1, 1);
loglog(w_ss, S_xx(:,1), '--', w_fft, I_fft(1:floor(N_fft/2),1), '-.', w_filt, I_filt(1:floor(N_filt/2),1))
axis(10.^[-2, 2, -15, 0])
ylabel('$S_{\hat{u}\hat{u}}$ [rad$^2$]', 'Interpreter', 'Latex'); title('\textbf{Power Spectral Densities of Aircraft States due to Vertical Turbulence}', 'Interpreter', 'Latex')
legend(ax1, 'Analytical', ['FFT, ' num2str(N_rl) ' Realizations'] , 'Filtered')
legend('boxoff')
grid on

subplot(5, 1, 2)
loglog(w_ss, S_xx(:,2), '--', w_fft, I_fft(1:floor(N_fft/2),2), '-.', w_filt, I_filt(1:floor(N_filt/2),2))
axis(10.^[-2, 2, -15, 0])
ylabel('$S_{\alpha\alpha}$ [rad$^2$]', 'Interpreter', 'Latex')
grid on

subplot(5, 1, 3)
loglog(w_ss, S_xx(:,3), '--', w_fft, I_fft(1:floor(N_fft/2),3), '-.', w_filt, I_filt(1:floor(N_filt/2),3))
axis(10.^[-2, 2, -15, 0])
ylabel('$S_{\theta\theta}$ [rad$^2$]', 'Interpreter', 'Latex')
grid on

subplot(5, 1, 4)
loglog(w_ss, S_xx(:,4), '--', w_fft, I_fft(1:floor(N_fft/2),4), '-.', w_filt, I_filt(1:floor(N_filt/2),4))
axis(10.^[-2, 2, -20, -5])
ylabel('$S_{\frac{q\overline{c}}{V}\frac{q\overline{c}}{V}}$ [rad$^2$]', 'Interpreter', 'Latex')
grid on

subplot(5, 1, 5)
loglog(w_ss, S_xx(:,5), '--', w_fft, I_fft(1:floor(N_fft/2),5), '-.', w_filt, I_filt(1:floor(N_filt/2),5))
axis(10.^[-2, 2, -10, 0])
xlabel('$\omega$ [rad/s]', 'Interpreter', 'Latex'); ylabel('$S_{n_zn_z}$ [rad$^2$]', 'Interpreter', 'Latex')
grid on