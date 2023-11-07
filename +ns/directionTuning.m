function out = directionTuning(direction,response,parms,npEstimate)
% Fit the response with the sum of two circular gaussians; one at the
% preferred, the other at the anti-preferred (+180) direction.
%
% This function is intended to be used with ns.Tuning.
%
% INPUT
% direction -  The direction use in each trial (in degrees) [nrTrials 1]
% response  -The response in each trial. [nrTrials 1]
% parms     - parms from ns.Tuning (not used here)
% npEstimate - the non-parametric estimates of the tuning curve as
%               determined by ns.Tuning, which are used to setup the
%               nonlinear parameter optimizatio in a reasonable starting
%               point.
% OUTPUT 
% A struct with the following fields:
% amplitude   # Response at the preferred direction. 
% preferred   # Preferred direction in degrees
% width       # Widht of the tuning at the preferred
% antiwidth   # Widht of the tuning at the anti preferred
% antiamplitude# Response at the anti-preferred direction.
% offset      # Offset response.
% ci          # Tuning parameters confidence intervals  [Offset Preferred Width AntiWidth Peak AntiPeak]
% error       # Error estimate [Offset Preferred Width AntiWidth Peak AntiPeak]
% splithalves # Quality of the fit based on split-halves cross validation
% r           # Correlation between fit and non-parametric estimate
% residualz   # Z-score of mean residuals (0=good)
% bootparms # Boostrap estimates of the parameters
%
% See also ns.twoHumps

npAntiAmplitude  = npEstimate.tc(npEstimate.tcx==npEstimate.preferred+180 | npEstimate.tcx==npEstimate.preferred-180);
fun = @ns.twoHumps;
%% Fit a parametric (twoHumps) function
bootfun = @(x,y) solve(fun,direction,response,npEstimate.min,npEstimate.preferred,npEstimate.peak-npEstimate.min,npAntiAmplitude);
bootOpts = statset;
bootOpts.UseParallel = ~isempty(gcp("nocreate")); % Use parallel pool only if the user has started it.
[bootEstimates]=bootstrp(parms.nrBoot,bootfun,npEstimates.tcx,signal,'Options',bootOpts);
% Determine preferred axis and circular standard deviation
% Multiply by two in case the preferred swaps
% between bootstrap sets (i.e. an orientation tuned
% neuron).
z = exp(2*pi/180*1i*bootEstimates(:,2));
R = mean(z);

preferredOrAntiPreferred = angle(R)/2; % Back to the original 360 space
% Because we always call the direction with the highest amplitude the
% preferred direction (see "SOLVE" above), we just have to count which occurs
% most frequently : preferredAxis or preferredAxis+180.
z = exp(pi/180*1i*bootEstimates(:,2));
inner = [real(z), imag(z)]*[cos(preferredOrAntiPreferred); sin(preferredOrAntiPreferred)];
if mean(inner>0)>0.5
    isPreferred = inner>0;
else
    isPreferred = inner<0;
end
bootEstimates(~isPreferred,[2 4 5]) = bootEstimates(~isPreferred,[2 5 4]) + [180 0 0];
bootEstimates(:,2)= mod(bootEstimates(:,2),360);% Convenience

estimate = mean(bootEstimates,1,'omitnan');
error = std(bootEstimates,0,1,"omitnan");

ci= prctile(bootEstimates,[parms.alpha/2 1-parms.alpha/2]);

% Correct circular vriables.
z = exp(pi/180*1i*bootEstimates(:,2));
R = mean(z);
estimate(2) = 180/pi*angle(R);
error(2) = 180/pi*sqrt(-2*log(abs(R))); % Circular standard deviation
% Use circular standard deviation as the
% confidence limits.
ci(:,2) = mod(estimate(2) + error(2)*[-1 1]',360);

fittedCurve = fun(uDirection,estimate);
gof = corr(tuningCurve,fittedCurve,'type','Spearman');

residuals = fun(direction,estimate)- spk;
meanR = mean(residuals,"all","omitnan");
stdR =  std(residuals,0,"all","omitnan");
residualZ = abs(meanR./stdR);

%% SplitHalves
r= nan(1,parms.nrSplitHalves);
for i=1:parms.nrSplitHalves
    [oneHalfTrials,otherHalfTrials] =resampleTrials(stimulusIx,false,0.5) ;
    % First half
    firstHalfEstimate = sbx.Tuning.solve(direction(oneHalfTrials),spk(oneHalfTrials),npMin,npPreferred,npAmplitude,npAntiAmplitude);
    otherHalfEstimate = sbx.Tuning.solve(direction(otherHalfTrials),spk(otherHalfTrials),npMin,npPreferred,npAmplitude,npAntiAmplitude);
    tc1 = sbx.Tuning.twoHumps(uDirection,firstHalfEstimate);
    tc2 = sbx.Tuning.twoHumps(uDirection,otherHalfEstimate);
    r(i) = corr(tc1,tc2);
end
splithalves  = mean(r);

%% Package into a struct 
out.offset = estimate(1);
out.preferred = estimate(2);
out.width = estimate(3);
out.antiWidth = estimate(4);
out.amplitude = estimate(5);
out.antiAmplitude = estimate(6);
out.ci = ci';
out.error = error';
out.splitHalves = splithalves;
out.r = gof;
out.bootParms = bootEstimates;
out.residualZ = residualZ;



end


function estimate = solve(fun,direction,spk,npMin,npPreferred,npAmplitude,npAntiAmplitude)
% Use the optimization toolbox to do a nonlinear least squares
% fit of the spiking data to the twoHumps function

% Setup the variables and bounds
offset = optimvar('offset',1,'LowerBound',0,'UpperBound',2*npMin);
preferred = optimvar('preferred',1,'LowerBound',npPreferred-45,'UpperBound',npPreferred+45);
sigma = optimvar('sigma',1,'LowerBound',0);
amplitude = optimvar('amplitude',1,'LowerBound',0,'UpperBound',2*npAmplitude);
antiAmplitude = optimvar('antiAmplitude',1,'LowerBound',0,'UpperBound',2*npAntiAmplitude);
antiSigma = optimvar('antiSigma',1,'LowerBound',0); % kappa= 30 is a 25 degree FWHM.

% Least square
objective = sum((spk  - fun(direction,[offset preferred sigma antiSigma amplitude antiAmplitude ])).^2);
problem = optimproblem("Objective",objective);

% Initial values
initialValue.offset = 0;
initialValue.sigma= 100;            % Start wide
initialValue.antiSigma = 100;
initialValue.amplitude= npAmplitude;  % Start with nonparameteric estiamtes
initialValue.antiAmplitude = npAntiAmplitude;
initialValue.preferred = npPreferred;

% Set options
opts =optimset;
opts.Display = 'off';
opts.Algorithm = 'trust-region-reflective';

% Solve and return as vector
solution =solve(problem,initialValue,'Options',opts);
estimate  = [solution.offset mod(solution.preferred,360) solution.sigma  solution.antiSigma solution.amplitude solution.antiAmplitude];
end


