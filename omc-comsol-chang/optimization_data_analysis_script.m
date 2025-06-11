clc, clear;
data_path = 'C:\Users\User\Documents\comsol\omc_comsol_sim\omc-comsol-chang\test\snowflake\05302025\';
matfiles = dir(fullfile([data_path,'*.mat']));
for i=1:length(matfiles)
    filename = matfiles(i).name;
    fullpathname = fullfile([data_path,filename]);
    data = importdata(fullpathname);
    bandGapSize = data.full.gapSize;
    midGap = data.full.midGap;
    fitness_list = calFitnessSnowflake(midGap,bandGapSize);
    [fitness,fitness_idx] = min(fitness_list);
    fprintf([filename,'\n']);
    fprintf('fitness: %f \n', fitness);
    fprintf('midGap: %f GHz \n',midGap(fitness_idx)/(1e9));
    fprintf('gapSize: %f GHz \n',bandGapSize(fitness_idx)/(1e9));
end

%% function definitions
function fitness = calFitnessSnowflake(fullMidBand,fullGapSize)
    targetFreq = 20e9; % the target frequency of the bandgap 
    freqTolerance = 1e9; % the tolerance for frequency offsets 
    % the fitness function associated with the optimization code
    fitness = -fullGapSize.*exp(-((fullMidBand-targetFreq)./freqTolerance).^2);
end