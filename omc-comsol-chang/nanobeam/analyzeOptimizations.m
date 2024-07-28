% script to analyze optimization runs for different starting points and
% iterations

clear optimds
close all
%% load files
%% opti6p_nd8g0wOQLapMax
% datLocLoad = 'D:\Files\OMC-SiV\45oOMC\opti6p_nd8g0wOQLapMax\';
% optimdsLoad = {'optimds20191027T230854.mat';...
%                'optimds20191028T191525.mat'};
% spLoad = {(1:2);...
%           (3:7)};
% geomStr = 'Tri 45o';
      
%% opti6p_nd9g0wOQLapMax
% datLocLoad = 'D:\Files\OMC-SiV\45oOMC\opti6p_nd9g0wOQLapMax\';
% optimdsLoad = {'optimds20191030T004246.mat'};
% spLoad = {(1:10)};
% geomStr = 'Tri 45o';
      
%% opti6p_nd9g0wOQLapMax
% datLocLoad = 'D:\Files\OMC-SiV\45oOMC\opti6p_nd9g0wOQLapMax_hi_a_low_w\';
% optimdsLoad = {'optimds20191102T232606.mat'};
% spLoad = {(1:11)};
% geomStr = 'Tri 45o';

%% RectOMC\opti_g0wOQLMax
% datLocLoad = 'D:\Files\OMC-SiV\RectOMC\opti_g0wOQLMax\';
% optimdsLoad = {'optimds20191022T012703.mat'};
% spLoad = {(1:5)};
% geomStr = 'Rect';

%% RectOMC\opti_g0wOQ(1e6)LapMean
% datLocLoad = 'D:\Files\OMC-SiV\RectOMC\opti_g0wOQ(1e6)LapMean\';
% optimdsLoad = {'optimds20191023T040019.mat'};
% spLoad = {(1:7)};
% geomStr = 'Rect';

%% RectOMC\opti_g0(2e5)wOQ(1e6)LapMean
% datLocLoad = 'D:\Files\OMC-SiV\RectOMC\opti_g0(2e5)wOQ(1e6)LapMean\';
% optimdsLoad = {'optimds20191024T215205.mat'};
% spLoad = {(1:8)};
% geomStr = 'Rect';

%% 35oOMC\opti6p_nd8g0wOQLapMax
datLocLoad = 'D:\Files\OMC-SiV\35oOMC\opti6p_nd8g0wOQLapMaxT2\';
optimdsLoad = {'optimds20200108T171917.mat'};
spLoad = {(3:6)};
geomStr = 'Tri 35o';
    
%%
datLocSave = datLocLoad;

% custom limit function - to handle cases where plot min and max limits are
% identical: add 10% either side of limit
eps = @(x) 0.1*mean(x)*(abs(min(x)-max(x))/mean(x)<1e-6);

%% plot results vs itr
for di = 1:length(optimdsLoad)
    clear optimds; load([datLocLoad,optimdsLoad{di}]);
    %%
    for sp = spLoad{di}
        figure('units','pixels','outerposition',[0 0 1920 1080])
        X = optimds.sp(sp).res(:,23); % iteration index
        
        % fitness
        subplot(4,4,1:4); box on; hold on; Y = optimds.sp(sp).res(:,22);
        plot(X,Y,'o-','Linewidth',2);
        xlabel('iteration no.'); ylabel('F'); xlim([min(X)-eps(X) max(X)+eps(X)]); ylim([min(Y)-eps(Y) max(Y)+eps(Y)]);
        [FMaxSP,FMaxSPIdx] = max(Y);
        
        % unit cell lattice constant
        subplot(4,4,5); box on; hold on; Y = optimds.sp(sp).res(:,2)*1e9;
        plot(X,Y,'o-','Linewidth',0.5,'MarkerSize',4);
        ylabel('a/nm'); xlim([min(X)-eps(X) max(X)+eps(X)]); ylim([min(Y)-eps(Y) max(Y)+eps(Y)]);
        
        % beam width
        subplot(4,4,9); box on; hold on; Y = optimds.sp(sp).res(:,5)*1e9;
        plot(X,Y,'o-','Linewidth',0.5,'MarkerSize',4);
        ylabel('w/nm'); xlim([min(X)-eps(X) max(X)+eps(X)]); ylim([min(Y)-eps(Y) max(Y)+eps(Y)]);
        
        % beam thickness
        subplot(4,4,13); box on; hold on; Y = optimds.sp(sp).res(:,1)*1e9;
        plot(X,Y,'o-','Linewidth',0.5,'MarkerSize',4);
        ylabel('th/nm'); xlim([min(X)-eps(X) max(X)+eps(X)]); ylim([min(Y)-eps(Y) max(Y)+eps(Y)]);
        
        % max defect
        subplot(4,4,6); box on; hold on; Y = optimds.sp(sp).res(:,8);
        plot(X,Y,'o-','Linewidth',0.5,'MarkerSize',4);
        ylabel('maxdef'); xlim([min(X)-eps(X) max(X)+eps(X)]); ylim([min(Y)-eps(Y) max(Y)+eps(Y)]);
        
        % unit cell hole x-height
        subplot(4,4,10); box on; hold on; Y = optimds.sp(sp).res(:,6)*1e9;
        plot(X,Y,'o-','Linewidth',0.5,'MarkerSize',4);
        ylabel('hx/nm'); xlim([min(X)-eps(X) max(X)+eps(X)]); ylim([min(Y)-eps(Y) max(Y)+eps(Y)]);
        
        % unit cell hole y-width
        subplot(4,4,14); box on; hold on; Y = optimds.sp(sp).res(:,7)*1e9;
        plot(X,Y,'o-','Linewidth',0.5,'MarkerSize',4);
        ylabel('hy/nm'); xlim([min(X)-eps(X) max(X)+eps(X)]); ylim([min(Y)-eps(Y) max(Y)+eps(Y)]);
        
        % taper oblong
        subplot(4,4,7); box on; hold on; Y = optimds.sp(sp).res(:,9);
        plot(X,Y,'o-','Linewidth',0.5,'MarkerSize',4);
        ylabel('oblong'); xlim([min(X)-eps(X) max(X)+eps(X)]); ylim([min(Y)-eps(Y) max(Y)+eps(Y)]);
        
        % optical wavelength
        subplot(4,4,11); box on; hold on; Y = optimds.sp(sp).res(:,10)*1e9;
        plot(X,Y,'o-','Linewidth',0.5,'MarkerSize',4);
        ylabel('\lambda_{opt}/nm'); xlim([min(X)-eps(X) max(X)+eps(X)]); ylim([min(Y)-eps(Y) max(Y)+eps(Y)]);
        
        % optical Q
        subplot(4,4,15); box on; hold on; Y = optimds.sp(sp).res(:,11)/1e6;
        plot(X,Y,'o-','Linewidth',0.5,'MarkerSize',4);
        ylabel('Q/10^6'); xlim([min(X)-eps(X) max(X)+eps(X)]); ylim([min(Y)-eps(Y) max(Y)+eps(Y)]);
        
        % mechanical frequency
        subplot(4,4,8); box on; hold on; Y = optimds.sp(sp).res(:,12)/1e9;
        plot(X,Y,'o-','Linewidth',0.5,'MarkerSize',4);
        ylabel('\omega_{mech}/GHz'); xlim([min(X)-eps(X) max(X)+eps(X)]); ylim([min(Y)-eps(Y) max(Y)+eps(Y)]);
        
        % optomechanical coupling
        subplot(4,4,12); box on; hold on; Y = abs(optimds.sp(sp).res(:,15))/1e3;
        Y2 = abs(optimds.sp(sp).res(:,13))/1e3;
        Y3 = abs(optimds.sp(sp).res(:,14))/1e3; 
        plot(X,Y,'o-','Linewidth',0.5,'MarkerSize',4);
        plot(X,Y2,'ro:','Linewidth',0.5,'MarkerSize',2);
        plot(X,Y3,'ko:','Linewidth',0.5,'MarkerSize',2);
        ylabel('g/kHz'); xlim([min(X)-eps(X) max(X)+eps(X)]); ylim([min([Y;Y2;Y3]) max([Y;Y2;Y3])]);
        legend({'tot','MB','PE'},'location','southoutside','orientation','horizontal')
        
        % orbital coupling
        subplot(4,4,16); box on; hold on; Y = abs(optimds.sp(sp).res(:,16))/1e6;
        Y2 = abs(optimds.sp(sp).res(:,20))/1e6;
        Y3 = abs(optimds.sp(sp).res(:,21))/1e6; 
        plot(X,Y,'o-','Linewidth',0.5,'MarkerSize',4);
        plot(X,Y2,'ro:','Linewidth',0.5,'MarkerSize',2);
        plot(X,Y3,'ko:','Linewidth',0.5,'MarkerSize',2);
        ylabel('\lambda_{SiV}/MHz'); xlim([min(X)-eps(X) max(X)+eps(X)]); ylim([min([Y;Y2;Y3]) max([Y;Y2;Y3])]);
        legend({'global max','aperture max','aperture mean'},'location','southoutside','orientation','horizontal')
        
        % global title
        nholes = optimds.sp(sp).res(1,3);
        ndef = optimds.sp(sp).res(1,4);
        optimdsName = strrep(optimdsLoad{di},'.mat','');
        optimdsName = strrep(optimdsName,'optimds','ID ');
        fitStr = strrep(optimds.fitStr,'omds.cpl.','');
        fitStr = strrep(fitStr,'299792458','c');
        fitStr = strrep(fitStr,'SiV_111.full','SiV111');
        fitStr = strrep(fitStr,'./','/');
        fitStr = strrep(fitStr,'.*','*');
        ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
        text(0.5, 1,{['\bf OMC optimization ',optimdsName,': ',...
                geomStr,', n_{hole} = ',num2str(nholes),', ', ...
                'n_{def} = ',num2str(ndef),', ',...
                'start point ',num2str(sp),', ',num2str(max(X)),' itrs, '...
                'max F = ',num2str(FMaxSP,'%.3e'),' at itr ',num2str(FMaxSPIdx)];...
                ['F = ',fitStr]},...
                'HorizontalAlignment' ,'center','VerticalAlignment', 'top')

        % save figure
        set(gcf,'PaperPositionMode','auto')
        pathFig = [datLocSave,'\optimResultsVsItr_SP',num2str(sp)];
        saveas(gcf,[pathFig,'.fig']);
        saveas(gcf,[pathFig,'.png']);
        close
    end
end