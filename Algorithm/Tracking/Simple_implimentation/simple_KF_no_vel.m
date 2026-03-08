clc; clear;

%% -- Read signal --
%{
filename = 'rx_stream-320.00-1975.00-2-nvme2.56-10db.bin';
tStart_ms = 0;
tEnd_ms   = 100;
[r, ~, ~, ~] = readSignal(filename, Fs, tStart_ms, tEnd_ms);
r = r(:);
%}

%% -- or Generate dummy signal --
Fs = 3.2e8;
Fsym = 2.304e8;
Tsim = 0.5;
Tcoh = 1e-3;          % 1 ms
Tpss = 400/Fsym;
Trep = 1e-3;
Nrep = round(Fs*Trep);
fd0_gen = 100000;
fddot = 2000;
N_gen = round(Fs*Tsim);
r = zeros(N_gen,1);
t = (0:N_gen-1).'/Fs;
snr_dB = 10;

[r_loc, loc] = localReplica(Fs, Fsym);
r_loc = r_loc';
Npss = length(r_loc);

for k = 0:floor((N_gen-1)/Nrep)
    i0 = k*Nrep + 1;
    i1 = min(i0 + Npss - 1, N_gen);
    r(i0:i1) = r(i0:i1) + r_loc(1:(i1-i0+1));
    Ps = mean(abs(r(i0:i1)).^2);
    sigma2 = Ps / (10^(snr_dB/10));
    w = sqrt(sigma2/2) * (randn(size(r(i0:i1))) + 1j*randn(size(r(i0:i1))));
    r(i0:i1) = r(i0:i1) + w;
end

fd = generateDoppler(fd0_gen, fddot, t, 0);
phi_sig = 2*pi*cumsum(fd)/Fs;
r = r.*exp(1j*phi_sig);

%% --- KF (1-state): fd only ---
T     = Tcoh;
N_coh = round(Fs*Tcoh);
K     = floor(length(r)/N_coh);
tblk  = (0:N_coh-1).'/Fs;

% State model
Fmat = 1;      % fd_{k+1} = fd_k
Hmat = 1;      % measure fd

% Init
x = fd(1);     % scalar
P = (25^2);    % scalar variance

% Process noise (random walk). Tune q_fd (Hz^2/s) roughly.
q_fd = 5;                % <-- tune
Q = q_fd * T;            % discrete-time RW variance

% Measurement noise
sigma_fd = 10;           % Hz
R = sigma_fd^2;

% Storage
doppler_est   = NaN(K,1);
fd_res_store  = NaN(K,1);
fd_meas_store = NaN(K,1);

t_k = ((0:K-1)*N_coh + (N_coh-1)/2).'/Fs;

% True Doppler per block (mean)
fd_true_blk = NaN(K,1);
for k = 1:K
    idx = (k-1)*N_coh + (1:N_coh);
    fd_true_blk(k) = mean(fd(idx));
end

z_prev = NaN;

%% --- Loop ---
for k = 1:K
    idx  = (k-1)*N_coh + (1:N_coh);
    xblk = r(idx);

    % Predict
    x = Fmat*x;
    P = Fmat*P*Fmat' + Q;

    % Wipe-off with predicted fd
    fd_pred = x;
    xw = xblk .* exp(-1j*2*pi*fd_pred*tblk);

    % Correlate
    rc = xcorr(xw, r_loc);
    [~, imax] = max(abs(rc));
    z = rc(imax);

    % Measurement from phase diff (needs previous z)
    if ~isnan(z_prev)
        dphi   = angle(z * conj(z_prev));
        fd_res = dphi/(2*pi*T);

        fd_res_store(k) = fd_res;

        fd_meas = fd_pred + fd_res;   % pseudo-measurement of fd
        fd_meas_store(k) = fd_meas;

        % Update (scalar KF)
        S  = Hmat*P*Hmat' + R;
        Kk = P*Hmat'/S;
        x  = x + Kk*(fd_meas - Hmat*x);
        P  = (1 - Kk*Hmat)*P;
    end

    z_prev = z;
    doppler_est(k) = x;
end

%% --- Plots ---
figure;
plot(t_k, doppler_est, 'r--'); hold on;
plot(t_k, fd_true_blk, 'k-');
legend('Estimated Doppler','True Doppler (block mean)');
xlabel('Time (s)'); ylabel('Doppler (Hz)'); grid on;

figure;
mask = ~isnan(fd_res_store);
plot(t_k(mask), fd_res_store(mask));
xlabel('Time (s)'); ylabel('Residual f_d (Hz)'); grid on;


%% Functions

function fd = generateDoppler(fd0, fddot, t, doptype)
    if doptype == 1
        idx = ones(1,length(t));
        fd = fd0.*idx; %+ fddot.*t;   % Hz
        fd = fd';
    else
        tpass = 600;
        fc = 2.9e9;
        c = 3*10^8;
        t = t(:);
        v = 7500;
        d = 400e3;
        tc=t - tpass/2;
        rho  = sqrt(d.^2 + (v.*tc).^2);
        rdot = (v.^2 .* tc) ./ rho;      % [m/s],
        fd   = -(fc/c) .* rdot;         % [Hz]
    end
end