function out = leo_geom_iekf(cfg)
%LEO_GEOM_IEKF Short-name wrapper for the geometry/handover IEKF solver.

out = real_strtg_sequential_EKF_geom(cfg);
end
