function out = leo_iekf(cfg)
%LEO_IEKF Short-name wrapper for the trusted 5-state IEKF solver.
%
% The implementation is intentionally delegated to the original solver file
% to preserve numerical behaviour.

out = real_strtg_sequential_EKF_fastrealistic(cfg);
end
