function K = KfromCube(x,y,z,a,b,c)

K = indefInt(x+a,y+b,z+c) - ...
    indefInt(x-a,y+b,z+c) - ...
    indefInt(x+a,y-b,z+c) + ...
    indefInt(x-a,y-b,z+c) - ...
    indefInt(x+a,y+b,z-c) + ...
    indefInt(x-a,y+b,z-c) + ...
    indefInt(x+a,y-b,z-c) - ...
    indefInt(x-a,y-b,z-c);

end

function K = indefInt(x,y,z)
% This function returns the kernal to the integral Integrate[R,Norm[R]^3],
% numerically evaluated at the specified points.

% Frequently used quantity, saved to repeat duplicating calculation:
A = sqrt(x.^2+y.^2+z.^2);

% Analytic solution for the kernal of the integral Integrate[R,Norm[R]^3].
% I double checked the results of this formula with the output of
% Mathematica, so I think it's ok.
K = [y - x.*atan(y./x) + x.*atan(y.*z./(x.*A)) - z.*log(2*(y+A)) - ...
    y.*log(2*(z+A)), ...
    z - y.*atan(z./y) + y.*atan(z.*x./(y.*A)) - x.*log(2*(z+A)) - ...
    z.*log(2*(x+A)), ...
    x - z.*atan(x./z) + z.*atan(x.*y./(z.*A)) - y.*log(2*(x+A)) - ...
    x.*log(2*(y+A))];


% Alternative solutions that are numerically slower:
% Bx = -z*acoth(y./A) - y*acoth(z./A) + x.*atan(y.*z./(x.*A));
% By = -x.*acoth(z./A) + y.*atan(x.*z./(y.*A)) - z.*log(2*(x+A));
end