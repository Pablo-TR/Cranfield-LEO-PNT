clear all; close all
%   Signal generator
Fs = 3.2e8;
Fsym = 2.304e8;
Tsim = 0.05;
Tcoh = 10^-2;
Tpss = 400/Fsym;
Trep = 1e-3;
Nrep = round(Fs*Trep);
fd0 = 1000000;
fddot = 2000;

N = round(Fs*Tsim);
s = zeros(N,1);
t = (0:N-1).'/Fs;

[r_loc, loc] = localReplica(Fs, Fsym);
r_loc = r_loc';

Npss = length(r_loc);

for k = 0:floor((N-1)/Nrep)
    i0 = k*Nrep + 1;
    i1 = min(i0 + Npss - 1, N);
    s(i0:i1) = s(i0:i1) + r_loc(1:(i1-i0+1));
end

fd = generateDoppler(fd0, fddot,t);

phi = 2*pi*cumsum(fd)/Fs;
r = s.*exp(1j*phi);

N_block = round(Tcoh*Fs);
t = (0:N-1).' / Fs;

%% Test
[r_local, ~] = localReplica(Fs, Fsym);
fd1 = fd0;
fd2 = -20000;
r_wiped1 = r.*exp(-1j*2*pi*fd1*t);
r_wiped2 = r.*exp(-1j*2*pi*fd2*t);
mf = conj(flipud(r_local(:)));
corr1 = conv(r_wiped1(:), mf);
corr2 = conv(r_wiped2(:), mf);
metric1 = max(abs(corr1))
metric2 = max(abs(corr2))
%%

[doppler,timeArr] = dopplerSearchCustom(r, Fs, Tsim, Tcoh, t);

%%
function [fd] = generateDoppler(fd0, fddot, t)
    fd = fd0; %+ fddot.*t;
end