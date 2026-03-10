clc; clear;

%% -- Read signal --
%{
filename = 'rx_stream-320.00-1975.00-2-nvme2.56-10db.bin';
tStart_ms = 0;
tEnd_ms   = 100;
[r, ~, ~, ~] = readSignal(filename, Fs, tStart_ms, tEnd_ms);
r = r(:);
K = floor(length(r)/N);
%}

%% -- or Generate dummy signal --
Fs = 3.2e8;
Fsym = 2.304e8;
Tsim = 0.5;
Tcoh = 10^-2;
Tpss = 400/Fsym;
Trep = 1e-3;
Nrep = round(Fs*Trep);
fd0_gen = 100000;
fddot = 2000;
N_gen = round(Fs*Tsim);
r = zeros(N_gen,1);
t = (0:N_gen-1).'/Fs;

[r_loc, loc] = localReplica(Fs, Fsym);
r_loc = r_loc';

Npss = length(r_loc);

for k = 0:floor((N_gen-1)/Nrep)
    i0 = k*Nrep + 1;
    i1 = min(i0 + Npss - 1, N_gen);
    r(i0:i1) = r(i0:i1) + r_loc(1:(i1-i0+1));
end

fd = generateDoppler(fd0_gen, fddot,t,0);
plot(fd, t)

phi = 2*pi*cumsum(fd)/Fs;
r = r.*exp(1j*phi);

%% --- Simulation parameters ---
Tcoh = 1e-2;      % 10 ms
Fs = 3.2e8;
T = Tcoh;          % dt del KF
phi = [1 T; 0 1];
C   = [1 0];
N_coh   = round(Fs*Tcoh);
tblk = (0:N_coh-1).'/Fs;
[r_local, ~] = localReplica(Fs, 2.304e8);
r_local = r_local(:);
K = floor(length(r)/N_coh);

%% --- EKF initialization ---
fd0 = fd(1);
fddot = (fd(2)-fd(1))/(t(2)-t(1)); 
x  = [fd0; fddot];                         % [fd; fd_dot]
P  = diag([ (10)^2, (10)^2 ]);       % 
phi = [1 T;
       0 1];
C = [1 0];
Q = (10^2)*[(T^3)/3 (T^2)/2; (T^2)/2 T];
R = 100^2; %Hz^2 
doppler_est = zeros(K,1);
z_prev = NaN;

%% --- KF loop with subaccumulation & wipeoff ---
for k = 1:K
    idx = (k-1)*N_coh + (1:N_coh);
    xblk = r(idx);
    % --- Wipe-off with predicted Doppler
    fd_pred = x(1);
    xw = xblk .* exp(-1j*2*pi*fd_pred*tblk);
    % --- Replica correlation
    rc = xcorr(xw, r_local);
    [~, imax] = max(abs(rc));
    z = rc(imax);     
   if ~isnan(z_prev)
        dphi = angle(z * conj(z_prev));    % [-pi, pi]
        fd_res = dphi/(2*pi*T);
        y = fd_pred + fd_res;    % medida ABSOLUTA
        Kk = P*C'/(C*P*C' + R);
        x  = x + Kk*(y - C*x);
        P  = (eye(2) - Kk*C)*P;
   end

    % --- Save estimates ---
    z_prev = z;
    doppler_est(k) = x(1);
    % --- Time Update ---
    x = phi*x;
    P = phi*P*phi' + Q;
end

%% --- Plot results ---

t_k = (0:K-1)*Tcoh;
fd_true_blk = fd(1 + (0:K-1)*N_coh);

figure;
plot(t_k, doppler_est, 'r--'); hold on;
plot(t_k, fd_true_blk, 'k-');
legend('Estimated Doppler','True Doppler');
xlabel('Time (s)'); ylabel('Doppler (Hz)'); grid on;

%% Functions
function fd = generateDoppler(fd0, fddot, t, doptype)
    if doptype == 1
        fd = fd0 + fddot.*t;   % Hz
    else
        tpass = t(end);
        fc = 2.9e9;
        c = 3*10^8;
        t = t(:);
        v = 7500;
        d = 1200e3;
        tc=t - tpass/2;
        rho  = sqrt(d.^2 + (v.*tc).^2);
        rdot = (v.^2 .* tc) ./ rho;      % [m/s], cambia de signo en t=0
        fd   = -(fc/c) .* rdot;         % [Hz]
    end
end