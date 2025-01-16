function F = readOpticalBandStruct(P,dataLoc)
    data = importdata(dataLoc);
    kNorm = str2double(data.textdata(:,1));
    nbands = length(kNorm)/(1+P.kpts);
    freq_list = data.data(:,1);
    F = ones((P.kpts+1),nbands);
    % sorting the data 
    for i=1:nbands
        for j=1:(P.kpts+1)
            F(j,i)=freq_list((j-1)*(P.kpts+1)+i)*1e12;
        end
    end
end