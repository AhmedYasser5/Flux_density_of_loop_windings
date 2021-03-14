function H = loopWinding(a, Z, I, A, z, tol)
% a is the radius of the loops
% Z is a vector of Z coordinates of the loops
% I is the current
% A is a vector of points used for X and Y coordinates
% z is the plane at which the magnetic field intensity is computed
% tol is the tolerance for calculation near the loop
% H is a vector containing the magnetic field intensity of a number of
% points having the third dimension = z
% WARNINGS may appear because of the limit of intervals for integral
% function used here, it can be avoided using the variable tol >= 1e-3.
% This happens because Biot-Savart Law can't be used for the magnetic flux
% near the source
M = length(A);
H = zeros(M, M);
N = length(Z);
minZ = inf;
for i = 1: N
    minZ = min(minZ, abs(z - Z(i)));
end
minZ = minZ^2;
for iter1 = 1: M
    x = A(iter1);
    for iter2 = 1: M
        y = A(iter2);
        if 0 & (x^2 + y^2) * (1 - a / sqrt(x^2 + y^2))^2 + minZ < tol
            H(iter1, iter2) = inf;
            continue
        end
        % a_x direction
        tmp = 0;
        for j = 1: N
            tmp = tmp + (z-Z(j)) .* integral(@(t)(cos(t)./((x-a.*cos(t)).^2+(y-a.*sin(t)).^2+(z-Z(j)).^2).^1.5), 0, 2*pi);
        end
        tmp = tmp .* a;
        H(iter1, iter2) = tmp.^2;
        % a_y direction
        tmp = 0;
        for j = 1: N
            tmp = tmp + (z-Z(j)) .* integral(@(t)(sin(t)./((x-a.*cos(t)).^2+(y-a.*sin(t)).^2+(z-Z(j)).^2).^1.5), 0, 2*pi);
        end
        tmp = tmp .* a;
        H(iter1, iter2) = H(iter1, iter2) + tmp.^2;
        % a_z direction
        tmp = 0;
        for j = 1: N
            tmp = tmp - integral(@(t)(sin(t).*(y-a.*sin(t))./((x-a.*cos(t)).^2+(y-a.*sin(t)).^2+(z-Z(j)).^2).^1.5), 0, 2*pi);
            tmp = tmp - integral(@(t)(cos(t).*(y-a.*cos(t))./((x-a.*cos(t)).^2+(y-a.*sin(t)).^2+(z-Z(j)).^2).^1.5), 0, 2*pi);
        end
        tmp = tmp .* a;
        H(iter1, iter2) = H(iter1, iter2) + tmp.^2;
        H(iter1, iter2) = sqrt(H(iter1, iter2)).* I ./ (4 .* pi);
    end
end
surf(A, A, H);
colormap jet;
end
% TRY loopWinding(0.02, 0, 1, [-0.2: 0.0025: 0.2], 0.1, 0.001);
% TRY loopWinding(0.02, 0, 1, [-0.1: 0.0025: 0.1], 0, 0.001);
% TRY loopWinding(0.02, [-1:0.1:1.5], 1, [-0.1: 0.005: 0.1], 0.25, 0.0001);
% TRY loopWinding(0.04, [-1:0.1:1.5], 1, [-0.1: 0.005: 0.1], 0, eps);