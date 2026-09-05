%% optimization_client
%The goal of this script is to take the cross bandstrcture code that we
%wrote and automate it. Few notes:
% 1. We're going to run this with no symmetry condition
% 2. We'll run it for kpts = 2


addpath([fullfile('.','test'),filesep])

% a, wc, hc

x0 = [450,378,95];
options = optimset('PlotFcns',@optimplotfval);

fun = @(x) optimize_cross(x);
x = fminsearch(fun,x0,options);

function F = optimize_cross(x)
clear P

vars.a = x(1);
vars.hc = x(2);
vars.wc = x(3);

disp('-----')
disp(vars)

P.th = 160e-9;
P.r1 = 30e-9;
P.r2 = 30e-9;
target_freq = 13e9;
P.datLoc = [fullfile('.','test'),filesep];

if min(vars.a-vars.hc,vars.wc) < 60
    disp('Fabrication intolerant');
    F = 0;
else
    results = solveBands(vars,P);
    
    if ~isempty(results.full.gapSize)
    mGaps = results.full.gapSize;
    mFreqs = results.full.midGap;

    gap_ind = find(abs(target_freq-mFreqs)<mGaps);
    if isempty(gap_ind)
        disp('Target frequency not in any gap, finding closest bandgap')
        gap_ind = find(min(abs(target_freq-mFreqs) - mGaps/2));  
    end
    
    disp(gap_ind)
    
    mGap = mGaps(gap_ind);
    mFreq = mFreqs(gap_ind);
    
    sigma = 1e10;
    freq_penalty = exp(-((target_freq-mFreq)/sigma)^2);
    disp(['Using gap number ', num2str(gap_ind)])
    disp(['Midgap Frequency = ', num2str(mFreq/1e9), ' GHz'])
    disp(['Gap size = ', num2str(mGap/mFreq*100), ' %'])
    disp(['Frequency penalty = ', num2str(freq_penalty)])

    mech_cont = mGap/mFreq*freq_penalty;
    F = -1*mech_cont;
    else 
        disp('No complete gaps')
        F = 0
    end
end
disp('-----')
end