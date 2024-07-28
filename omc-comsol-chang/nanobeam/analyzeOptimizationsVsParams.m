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

dataLoad = {'D:\Files\OMC-SiV\35oOMC\opti6p_nd8g0wOQLapMax\',...
                'optimds20200103T223303.mat',(1:4);...
            'D:\Files\OMC-SiV\35oOMC\opti6p_nd8g0wOQLapMax2\',...
                'optimds20200105T163835.mat',(1:4);...
            'D:\Files\OMC-SiV\35oOMC\opti6p_nd8g0wOQLapMax3\',...
                'optimds20200107T002626.mat',(1:5)};

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

for ni = 1:npath
    clear optimds; load([dataLoad{ni,1},dataLoad{ni,2}]);
    spAll = dataLoad{ni,3};
    for sp = spAll
        spTot = spTot+1;
        res = optimds.sp(sp).res;
        [itrSP,~] = size(res);
        itrTot = itrTot + itrSP;
        for fi = 1:6
            if ~hdl(fi)
                hdl(fi) = figure('units','pixels','outerposition',[0 0 800 600]);
            else
                figure(hdl(fi))
            end
            
            for subi = 1:4
                subplot(2,2,subi); hold on; box on;
                if subi == 1
                    X = res(:,2)*1e9; Y = res(:,5)*1e9;
                    xlbl = 'a/nm'; ylbl = 'w/nm';
                elseif subi == 2
                    X = res(:,2)*1e9; Y = res(:,1)*1e9;
                    xlbl = 'a/nm'; ylbl = 'th/nm';
                elseif subi == 3
                    X = res(:,6)./res(:,2); Y = res(:,7)./res(:,5);
                    xlbl = 'hx/a'; ylbl = 'hy/w';
                elseif subi == 4
                    X = res(:,8); Y = res(:,9);
                    xlbl = 'maxdef'; ylbl = 'oblong';
                end
                
                plot3(X,Y,res(:,resIdxs(fi)),'-','Color',[170 170 170]./256)
                scatter3(X,Y,res(:,resIdxs(fi)),30,res(:,resIdxs(fi)),'filled')
                xlabel(xlbl); ylabel(ylbl)
            end
        end
    end
    
    for fi = 1:6
        figure(hdl(fi));
%         disp(ni)
%         disp(cbarTitle)
        if ~savePlotEachDir && ni == npath
            cbarTitle = zeros(1,6);
        end
        if ~cbarTitle(fi)
            
            subplot(2,2,4)
            hold on; box on;    
            colorbar('Position',[0.925 0.15 0.03 0.75])
            
            ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
            text(0.5, 1,{['\bf ',resHdrs{fi},' vs geometry params, total ',...
                          num2str(itrTot),' itrs across ',...
                          num2str(spTot),' start points']},...
                'HorizontalAlignment' ,'center','VerticalAlignment', 'top')
            

            set(gcf,'PaperPositionMode','auto')
        end
        cbarTitle(fi) = 1;
        
        if savePlotEachDir
            datLocSave = dataLoad{ni,1};
            pathFig = [datLocSave,'optim',resFnames{fi},'VsGeomParams'];
            saveas(gcf,[pathFig,'.fig']);
            saveas(gcf,[pathFig,'.png']);
            close
        else
            datLocSave = datLoc;
            pathFig = [datLocSave,'optim',resFnames{fi},'VsGeomParams'];
            if ni==npath
                saveas(gcf,[pathFig,'.fig']);
                saveas(gcf,[pathFig,'.png']);
                close
            end
        end
    end
    if savePlotEachDir
        % reset counters
        hdl = zeros(1,6);
        spTot = 0;
        itrTot = 0;
        cbarTitle = zeros(1,6);
    end
end