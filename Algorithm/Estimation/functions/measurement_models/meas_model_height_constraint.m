function [zpred, H] = meas_model_height_constraint(x, iDU, nx)
H = zeros(1,nx);
H(iDU) = 1;
zpred = x(iDU);
end
