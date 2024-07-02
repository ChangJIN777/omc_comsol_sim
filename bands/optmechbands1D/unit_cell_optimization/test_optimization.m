%% test_optimization

addpath('C:\Users\mhaas\Documents\GitHub\device_design\bands\optmechbands1D')

% vars.a = 640;
% vars.w = 870;
% vars.hx = 220;
% vars.hy = 570;

% x0 = [608, 891, 262, 597];
x0 = [700,700,300,500];
options = optimset('PlotFcns',@optimplotfval);

fun = @(x) optimize_unit(x);
x = fminsearch(fun,x0,options);

function F = optimize_unit(x)
vars.a = x(1);
vars.w = x(2);
vars.hx = x(3);
vars.hy = x(4);

disp('-----')
disp(vars)
datLoc = 'E:\michael\omc-comsol-master-files\opt_v3';

if min(vars.a-vars.hx,vars.w-vars.hy) < 100
    disp('Fabrication intolerant');
    F = 0;
else
    results = solveOptMechBands(vars,datLoc);

    oGap = results.opt_bds.gapSize;
    oFreq = results.opt_bds.midGap;
    
    % This is kind of a sneaky way to optimize bandgap around a target
    % frequency, just take the minimum of 
    target_freq = 194e12;
    sigma = 15e12;
    freq_penalty = exp(-((target_freq-oFreq)/sigma)^2);
    
%         mGap = results.mech_bds.gapSize;
%         mFreq = results.mech_bds.midGap;
%         mech_cont = min(mGap/mFreq,.4);
    opt_cont = oGap/oFreq*freq_penalty;
    F = -1*opt_cont;

 
end
disp(F)
disp('-----')
end