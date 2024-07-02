% script to analyze optimization runs for different geometry parameters
clear optimds
close all
%% load files
% specify as cell array with each row containing directory, optimization
% data structure, and starting points

%% tri omcs
% dataLoad = {'D:\Files\OMC-SiV\45oOMC\opti6p_nd8g0wOQLapMax\',...
%                 'optimds20191027T230854.mat',(1:2);...
%             'D:\Files\OMC-SiV\45oOMC\opti6p_nd8g0wOQLapMax\',...
%                 'optimds20191028T191525.mat',(3:7);...
%             'D:\Files\OMC-SiV\45oOMC\opti6p_nd9g0wOQLapMax\',...
%                 'optimds20191030T004246.mat',[1:8,10];...
%             'D:\Files\OMC-SiV\45oOMC\opti6p_nd9g0wOQLapMax_hi_a_low_w\',...
%                 'optimds20191102T232606.mat',(1:11)};
% 
% datLoc = 'D:\Files\OMC-SiV\45oOMC\';

% dataLoad = {'D:\Files\OMC-SiV\35oOMC\opti6p_nd8g0wOQLapMax\',...
%                 'optimds20200103T223303.mat',(1:4);...
%             'D:\Files\OMC-SiV\35oOMC\opti6p_nd8g0wOQLapMax2\',...
%                 'optimds20200105T163835.mat',(1:4);...
%             'D:\Files\OMC-SiV\35oOMC\opti6p_nd8g0wOQLapMax3\',...
%                 'optimds20200107T002626.mat',(1:5)};

dataLoad = {'D:\Files\OMC-SiV\35oOMC\opti6p_nd8g0wOQLapMaxT2\',...
                 'optimds20200108T171917.mat',(1:6)};
datLoc = 'D:\Files\OMC-SiV\35oOMC\';

%% rect omcs
% dataLoad = {'D:\Files\OMC-SiV\RectOMC\opti_g0wOQLMax\',...
%                 'optimds20191022T012703.mat',(1:5);...
%             'D:\Files\OMC-SiV\RectOMC\opti_g0wOQ(1e6)LapMean\',...
%                 'optimds20191023T040019.mat',(1:7);...
%             'D:\Files\OMC-SiV\RectOMC\opti_g0(2e5)wOQ(1e6)LapMean\',...
%                 'optimds20191024T215205.mat',(1:8)};
% 
% datLoc = 'D:\Files\OMC-SiV\RectOMC\';

%%
[npath,nR] = size(dataLoad);
savePlotEachDir = 0;

% plot one result per fig
% results: optwvl, Q, wM, gOM, LSiV, F
% each fig: plot for different geom params 
% params: w, a, hx/a, hy/w, dmax, obl
% include th, w/th for rect

resHdrs = {'\lambda_{opt}','Q','\omega_M','g_{OM}','\lambda_{SiV,max}','F'};
resIdxs = [10,11,12,15,16,22];
resFnames = {'optWvl','Q','wM','gOM','LSiVMax','F'};

hdl = zeros(1,6);
spTot = 0;
itrTot = 0;
cbarTitle = zeros(1,6);
resAll = [];
% resNew = [];

for ni = 1:npath
    clear optimds; load([dataLoad{ni,1},dataLoad{ni,2}]);
    spAll = dataLoad{ni,3};
    for sp = spAll
        spTot = spTot+1;
        res0 = optimds.sp(sp).res;
        [itrSP,ncols] = size(res0);
        itrTot = itrTot + itrSP;
        resNew = [];
        
        
        fdtdlist = dir([dataLoad{ni,1},'sp',num2str(sp),'_*Vmode*.mat']);
%         disp(length(fdtdlist))
        if ~isempty(fdtdlist)
            % add FDTD results into optimds if not already done
            hdrsF = [optimds.hdrs(1:21),'lambdaFDTD','Vmode','Qtime',...
                             'Qx1','Qx2','Qy1','Qy2','Qz1','Qz2',...
                             'Qsc','Qwg','Qt','Trans',optimds.hdrs(22:23)];
            optimds.hdrsF = hdrsF;
            for nItr = 1:itrSP
                clear ds
%                 disp(nItr)
                dlist = dir([dataLoad{ni,1},'sp',num2str(sp),'_itr',num2str(nItr),'_*Vmode*.mat']);
                load([dataLoad{ni,1},dlist(1).name]);
                
                resNew(nItr,1:36) = [res0(nItr,1:21),ds.lambda,ds.Vmode,...
                                     ds.Qtime,ds.Qx1,ds.Qx2,ds.Qy1,ds.Qy2,...
                                     ds.Qz1,ds.Qz2,ds.Qsc,ds.Qwvg,ds.Qt,...
                                     ds.Trans,res0(nItr,22:23)];
            end
            optimds.sp(sp).resNew = resNew;
%             save([dataLoad{ni,1},'new_',dataLoad{ni,2}],'optimds');
            res = resNew;
        else
            res = res0;
        end
        
        % add new columns for directory and start point indices
        res(:,end+1) = sp;
        res(:,end+1) = ni;
        
        
        % recalc fitness function
        % tri: normmin(g*wO,1e5/196e12)*normmin(Q,Qcutoff)*LSiVAperMean/wM*optWvl/1550nm
        Fnew = min(abs(res(:,15).*res(:,10)./299792458),2.0e5/1.96e14)./(2.0e5/1.96e14)...
                .*min(res(:,11),1e6)./1e6...
                .*real(res(:,20))./1e6...
                .*res(:,10)./1550e-9...
                .*(res(:,34)>0.25).*(res(:,25)>0).*(res(:,26)>0).*(res(:,27)>0)...
                .*(res(:,29)>0).*(res(:,30)>0);
        
        % rect:
        % normmin(g*wO,2e5/196e12)*normmin(Q,3e6)*LSiVAperMean/wM*optWvl/1550nm
%         Fnew = min(abs(res(:,15).*res(:,10)./299792458),2.0e5/1.96e14)./(2.0e5/1.96e14)...
%                 .*min(res(:,11),3e6)./3e6...
%                 .*real(res(:,21))./1e6./(res(:,12)./1e9)...
%                 .*res(:,10)./1550e-9;
        
        res(:,end+1) = Fnew;
        
        % compile all results
        resAll = [resAll;res];
    end
end

% sort compiled results
[~,idxSort] = sort(resAll(:,end),1,'descend');
resAllSort = resAll(idxSort,:);