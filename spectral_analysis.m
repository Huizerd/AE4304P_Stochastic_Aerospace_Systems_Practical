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
    u_ss  = 2;
    S_uu  = sigma_ug^2;    % intensity W is equal to variance!
    filet = '\textbf{Power Spectral Densities of Aircraft States due to Horizontal Turbulence}';
    filen = 'C:\Users\jesse\Google Drive\CnS\AE4304P_Stochastic_Aerospace_Systems_Practical\figures\spectral_hor';
elseif sigma_ug == 0 && sigma_wg > 0
    u_ss  = 3;
    S_uu  = sigma_wg^2;     % intensity W is equal to variance!
    filet = '\textbf{Power Spectral Densities of Aircraft States due to Vertical Turbulence}';
    filen = 'C:\Users\jesse\Google Drive\CnS\AE4304P_Stochastic_Aerospace_Systems_Practical\figures\spectral_vert';
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

% Include automatic figure saving? For both turbulence cases?

set(0, 'DefaultAxesTickLabelInterpreter','latex')
set(0, 'DefaultLegendInterpreter','latex')
set(0, 'DefaultFigurePosition', [152.5 168 719 791.5])

colors = get(gca, 'ColorOrder');

ax1 = subplot(5, 1, 1);
h1 = loglog(w_filt, I_filt(1:floor(N_filt/2),1), 'Color', colors(3,:));
hold on
h2 = loglog(w_fft, I_fft(1:floor(N_fft/2),1), '-.', 'Color', colors(2,:));
h3 = loglog(w_ss, S_xx(:,1), '--', 'Color', colors(1,:));
hold off
axis(10.^[-2, 2, -15, 0])
ylabel('$S_{\hat{u}\hat{u}}$ $\Big[\frac{\mbox{rad$^2$}}{\mbox{rad/s}}\Big]$', 'Interpreter', 'Latex'); title(filet, 'Interpreter', 'Latex')
legend('Location', [0.782 0.887 0 0])
legend([h3 h2 h1], 'Analytical', ['FFT, ' num2str(N_rl) ' Realizations'] , 'Filtered');
legend('boxoff')
grid on

subplot(5, 1, 2)
loglog(w_filt, I_filt(1:floor(N_filt/2),2), 'Color', colors(3,:))
hold on
loglog(w_fft, I_fft(1:floor(N_fft/2),2), '-.', 'Color', colors(2,:))
loglog(w_ss, S_xx(:,2), '--', 'Color', colors(1,:))
hold off
axis(10.^[-2, 2, -15, 0])
ylabel('$S_{\alpha\alpha}$ $\Big[\frac{\mbox{rad$^2$}}{\mbox{rad/s}}\Big]$', 'Interpreter', 'Latex')
grid on

subplot(5, 1, 3)
loglog(w_filt, I_filt(1:floor(N_filt/2),3), 'Color', colors(3,:))
hold on
loglog(w_fft, I_fft(1:floor(N_fft/2),3), '-.', 'Color', colors(2,:))
loglog(w_ss, S_xx(:,3), '--', 'Color', colors(1,:))
hold off
axis(10.^[-2, 2, -15, 0])
ylabel('$S_{\theta\theta}$ $\Big[\frac{\mbox{rad$^2$}}{\mbox{rad/s}}\Big]$', 'Interpreter', 'Latex')
grid on

subplot(5, 1, 4)
loglog(w_filt, I_filt(1:floor(N_filt/2),4), 'Color', colors(3,:))
hold on
loglog(w_fft, I_fft(1:floor(N_fft/2),4), '-.', 'Color', colors(2,:))
loglog(w_ss, S_xx(:,4), '--', 'Color', colors(1,:))
hold off
axis(10.^[-2, 2, -20, -5])
ylabel('$S_{\frac{q\overline{c}}{V}\frac{q\overline{c}}{V}}$ $\Big[\frac{\mbox{rad$^2$}}{\mbox{rad/s}}\Big]$', 'Interpreter', 'Latex')
grid on

subplot(5, 1, 5)
loglog(w_filt, I_filt(1:floor(N_filt/2),5), 'Color', colors(3,:))
hold on
loglog(w_fft, I_fft(1:floor(N_fft/2),5), '-.', 'Color', colors(2,:))
loglog(w_ss, S_xx(:,5), '--', 'Color', colors(1,:))
hold off
axis(10.^[-2, 2, -10, 0])
xlabel('$\omega$ [rad/s]', 'Interpreter', 'Latex'); ylabel('$S_{n_zn_z}$ $\Big[\frac{\mbox{rad$^2$}}{\mbox{rad/s}}\Big]$', 'Interpreter', 'Latex')
grid on

set(gcf, 'Renderer', 'Painters')
savefig([filen '.fig'])
print('-painters', '-depsc', [filen '.eps'])