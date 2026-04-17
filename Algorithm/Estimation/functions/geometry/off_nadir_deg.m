function alpha = off_nadir_deg(rSat, rRec)
look  = rRec - rSat;
nadir = -rSat;
alpha = real(acosd(max(-1,min(1, dot(look,nadir)/(norm(look)*norm(nadir))))));
end
