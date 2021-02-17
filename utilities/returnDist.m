function distL = returnDist(J1,J2)
distL = sqrt(sum((J1-J2).^2));
end