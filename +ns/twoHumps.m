function [y] = twoHumps(x,parms)
% Use this to fit a direction selective tuning function with
% one bump at the preferred, and a (potentially smaller) bump
% at the anti-preferred (=preferred +180).
%
% I tried this with the sum of two Von Mises functions:
% y = amp1*exp(kappa1*cos((x - preferred))) +amp2*exp(kappa2*cos((x - preferred-180)))+offset;
% but did not work well (solutions were oddly always pretty
% poor). Currently the function uses two wrapped/circular
% gaussians instead, which seems to work better.
%
% INPUT
% x - Stimulus angles in degrees.
% parms - parameter vector - [offset, preferred , width1, width2, amp1 ,amp2]
%
% OUTPUT
%  y - Value per x
% 
% See also ns.directionTuning, ns.Tuning

offset = parms(1); preferred = parms(2); sigma1 = parms(3); sigma2 = parms(4); amp1=parms(5); amp2=parms(6);

% Von Mises attempt
%deg2rad =pi/180;
%term1 = amp1*exp(sigma1*cos(deg2rad*(x-preferred)));
%term2 = amp2*exp(sigma2*cos(deg2rad*(x-preferred-180)));
%y = term1 + term2 + offset;

y = offset;
for k=-4:4
    y = y + amp1*exp(-(1/sigma1^2)*(x-preferred-360*k).^2);
    y = y + amp2*exp(-(1/sigma2^2)*(x-preferred-180-360*k).^2);
end
end