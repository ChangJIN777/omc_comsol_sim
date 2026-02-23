% Draw 1/8 of the nanobeam cavity (using all available symmetry planes).
% The structure P is assumed to have the following fields, which define the
% geometry of the nanobeam cavity:
%
% P.a = lattice constant;
% P.wid = beam width; P.th = beam thickness;
% P.hx = hole height;
% P.hy = hole width;
% P.nholes = # of holes in 1/2 beam length;
% P.ndef = # of holes in 1/2 the defect region;
% P.maxdef = max. fractional increase in spacing/hole size in the defect
% P.oblong = hy is scaled by alpha^(1+oblong), hx by alpha^(1-oblong), e.g.
%   oblong == 1 results in constant height but changing width in defect;
%
% M.J. Burek, 08/16

function PlotDefectCells(P)
%% Geometry parameters

%a = P.ahole;
hx = P.geom(:,1);
hy = P.geom(:,2);
xpos = P.geom(:,3);
ypos = P.geom(:,4);
asym = P.asym;

if isfield(P, 'ctrholeadd') && (P.ctrholeadd>0)
    nctr = P.ctrholeadd;
else
    nctr = 0;
end

% no. of hole for each type (waveguide mirror, defect, mirror)
nwgmR = P.wvgmir;
if isfield(P,'phonMir') && P.nphonmir>0
    nphmR = P.nphonmir;
else
    nphmR = 0;
end

ndefR = P.ndef + nctr;
nmirR = P.nholes - P.ndef;
ntotR = nwgmR + nphmR + ndefR + nmirR;

if isfield(P,'PL') && P.asymCav
    PL = P.PL;
    if isfield(PL, 'ctrholeadd') && (PL.ctrholeadd>0)
        nctrL = PL.ctrholeadd;
    else
        nctrL = 0;
    end
    if isfield(PL,'phonMir') && PL.nphonmir>0
        nphmL = PL.nphonmir;
    else
        nphmL = 0;
    end
    nwgmL = PL.wvgmir;
    ndefL = PL.ndef + nctrL;
    nmirL = PL.nholes - PL.ndef;
else
    nphmL = nphmR;
    nwgmL = nwgmR;
    ndefL = ndefR;
    nmirL = nmirR;
end

ntotL = nwgmL + nphmL + ndefL + nmirL - P.holeatctr;

%% Plot defect region
figure; %set(gcf,'position',[9 1108-500 913 500])
ax = axes('position',[0.1 0 0.85 0.45]);
hold(ax,'on')

% assemble center positons for each mirror segment
xdat = 0.5*(xpos(2:end)'+xpos(1:end-1)');
a_hole = abs(xpos(2:end)'-xpos(1:end-1)');
hx_hole = hx';
hy_hole = hy';

hdl = plot(ax,xdat,ones(length(xdat)),'--r','LineWidth',1.5);
hdl(2) = plot(ax,xdat,a_hole/P.a,'-ks','LineWidth',1, ...
    'MarkerSize',6,'Markerfacecolor','k');

% assemble center positons for each hole position
hdl(3) = plot(ax,xpos,hx_hole/P.a,'-s','LineWidth',1, ...
    'color',[0 102 51]/255,'MarkerSize',6,'Markerfacecolor',[0 102 51]/255);
hdl(4) = plot(ax,xpos,hy_hole/P.a,'-bs','LineWidth',1, ...
    'MarkerSize',6,'Markerfacecolor','b');

axis(ax,'tight')
axl = axis(ax);
ylim(ax,[0, min([2,0.05+axl(4)])]);
xlabel(ax,'mirror segment #','FontName','arial','FontSize',12);
ylabel(ax,'(a, hx, hy)/a_{nom}','FontName','arial','FontSize',12);

legend([hdl(2) hdl(3) hdl(4)],'a','h_x','h_y','location','southoutside','orientation','horizontal');
set(ax,'FontName','arial','FontSize',10,'xtick',[],'xticklabel','');

hold(ax,'off');
box(ax,'on');


%% Plot geometry in the xy-plane

bx = axes('position',[0.1 0.4 0.85 0.45]);
% bx = axes('position',[0.1 0.1 0.85 0.85]);
hold(bx,'on')

%waveguide overlay
if isfield(P, 'tapSup') && isfield(P.tapSup, 'vertices')
    vertices = P.tapSup.vertices;
    gXL = [xpos(1)-P.a/2; vertices(:,1); xpos(end)+P.a/2];
    gYL = [-P.w/2; -1*vertices(:,2); -P.w/2];
    gXU = [xpos(1)-P.a/2; vertices(:,1); xpos(end)+P.a/2];
    gYU = [P.w/2; vertices(:,2); P.w/2];
else
    if isfield(P,'phonMir') && P.nphonmir>0
        gXL = [xpos(1)-P.phonMir.a/2 xpos(nphmL+nwgmL)+P.phonMir.hx/2 ...
               xpos(nwgmL+nphmL+1)-P.a/2 xpos(ntotL+ndefR+nmirR)+P.a/2 ...
               xpos(ntotL+ndefR+nmirR+1)-P.phonMir.hx/2 ...
               xpos(end)+P.phonMir.a/2];
        gYL = -1/2*[P.phonMir.w P.phonMir.w ...
                    P.w P.w ...
                    P.phonMir.w P.phonMir.w ];
        gXU = gXL;
        gYU = -gYL; 
    else
        gXL = [xpos(1)-P.a/2 xpos(end)+P.a/2];
        gYL = 1/2*[-P.w -P.w];
        gXU = [xpos(1)-P.a/2 xpos(end)+P.a/2];
        gYU = 1/2*[P.w P.w]; 
    end
end

plot(bx,gXU,gYU,'k','linewidth',1)
plot(bx,gXL,gYL,'k','linewidth',1)

%airhole overlay
for j=1:length(hy);
    xdat = linspace(-hx(j)/2,hx(j)/2,100)';
    ydatu(:,j) = sqrt(1-linspace(-hx(j)/2,hx(j)/2,100).^2/(hx(j)/2)^2)*hy(j)/2 + ypos(j);
    ydatd(:,j) = -sqrt(1-linspace(-hx(j)/2,hx(j)/2,100).^2/(hx(j)/2)^2)*hy(j)/2 + ypos(j);
    plot(bx,[xdat + xpos(j),xdat + xpos(j)], ...
        [ydatu(:,j),ydatd(:,j)],'k','linewidth',1);
    
    % lines in defect region
    if (j >= nphmL+nwgmL+nmirL+1 && j <= ntotL) || ...
       (j >= ntotL+1 && j <= ntotL+ndefR) 
        plot(bx,[xpos(j),xpos(j)],2*P.w*[-1,1], ...
            'k','linestyle',':','linewidth',.5);
    
    % lines at start of mirror region
    elseif j == nphmL+nwgmL+1 || j == nphmL+nwgmL+nmirL || ...
           j == ntotL+ndefR+1 || j == ntotL+ndefR + nmirR
        plot(bx,[xpos(j),xpos(j)],2*P.w*[-1,1], ...
            'b','linestyle','--','linewidth',.5);
    
    % lines at end of mirror region
%     elseif j == P.wvgmir+1 || j == (length(hy)-P.wvgmir)
%         plot(bx,[xpos(j),xpos(j)],2*P.w*[-1,1], ...
%             'b','linestyle','--','linewidth',.5);
    end
end
plot(bx,[xpos(ntotL+1),xpos(ntotL+1)],2*P.w*[-1,1], ...
    'r','linestyle','--','linewidth',.5);
set(bx,'xtick',[],'ytick',[]);
hold(bx,'off');
box(bx,'on')
set(bx,'xcolor','w','ycolor','w');
daspect([1 1 1])
xlim(bx,[min(gXU) max(gXU)])
bxl = axis(bx);
% ylim(bx,[bxl(1)/5 bxl(2)/5]);
ylim(bx,[-0.1*(bxl(2)-bxl(1)) 0.1*(bxl(2)-bxl(1))]);
linkaxes([bx,ax],'x')

% title
if strcmp(P.xsect,'tri')
    ttl1 = ['\theta=',num2str(P.theta,'%.0f'),'^o, '];
elseif strcmp(P.xsect,'rect')
    ttl1 = ['th=',num2str(P.th*1e9,'%.0f'),'nm, '];
elseif strcmp(P.xsect,'isoFit')
    ttl1 = ['th=',num2str(P.th*1e9,'%.0f'),'nm, '];
end
ttl1 = [ttl1,...
        'a=',num2str(P.a*1e9,'%.0f'),'nm, ',...
        'w=',num2str(P.w*1e9,'%.0f'),'nm, ',...
        'hx=',num2str(P.hx*1e9,'%.0f'),'nm, ',...
        'hy=',num2str(P.hy*1e9,'%.0f'),'nm, ',...
        'obl=',num2str(P.oblong,'%.4f'),', ',...
        'd_{max}=',num2str(P.maxdef,'%.5f')];

if isfield(P, 'ctrholeadd')
    ctrholeStr = [', n_{ctr} = ',num2str(P.ctrholeadd+1,'%.0f')];
else
    ctrholeStr = '';
end

if isfield(P,'taperTo') && strcmp(P.taperTo,'custom')
    ctrholeStr = [ctrholeStr,', aC=',num2str(P.a_ctr*1e9,'%.0f'),'nm, ',...
                      'hxC=',num2str(P.hx_ctr*1e9,'%.0f'),'nm, ',...
                      'hyC=',num2str(P.hy_ctr*1e9,'%.0f'),'nm'];
end

if isfield(P,'wgmTaper') && P.wvgmir > 0
    if strcmp(P.wgmTaper.func,'linear')
        fval = 'lin';
    elseif strcmp(P.wgmTaper.func,'quadratic')
        fval = 'quad';
    elseif strcmp(P.wgmTaper.func,'cubic')
        fval = 'cub';
    end
    wgmStr = [', ',fval, ' end taper, '];
    if strcmp(P.wgmTaper.endtype,'custom')
        wgmStr = [wgmStr,'a_{end}=',num2str(P.wgmTaper.a_end*1e9,'%.0f'),'nm, ',...
                         'hx_{end}=',num2str(P.wgmTaper.hx_end*1e9,'%.0f'),'nm, ',...
                         'hy_{end}=',num2str(P.wgmTaper.hy_end*1e9,'%.0f'),'nm'];
    elseif strcmp(P.wgmTaper.endtype,'maxdef')
        wgmStr = [wgmStr,'d_{max,end}=',num2str(P.wgmTaper.maxdef_end*1e9,'%.5f'),', ',...
                         'obl_{end}=',num2str(P.wgmTaper.oblong_end*1e9,'%.4f')];
    end
else
    wgmStr = '';
end

if isfield(P,'PL') && P.asymCav
    RLbl = 'R: ';
else
    RLbl = '';
end

ttl2 = [RLbl,'n_{h}=',num2str(P.nholes,'%.0f'),', ',...
        'n_{d}=',num2str(P.ndef,'%.0f'),', ',...
        'n_{wgm}=',num2str(P.wvgmir,'%.0f'),...
        ctrholeStr,wgmStr];
TAll = {ttl1;ttl2};

if isfield(P,'PL') && P.asymCav
    PL = P.PL;
    if isfield(PL, 'ctrholeadd')
        ctrholeLStr = [', n_{ctr} = ',num2str(PL.ctrholeadd+1,'%.0f')];
    else
        ctrholeLStr = '';
    end

    if isfield(PL,'wgmTaper') && PL.wvgmir > 0
        if strcmp(PL.wgmTaper.func,'linear')
            fval = 'lin';
        elseif strcmp(PL.wgmTaper.func,'quadratic')
            fval = 'quad';
        elseif strcmp(PL.wgmTaper.func,'cubic')
            fval = 'cub';
        end
        wgmLStr = [', ',fval, ' wgm taper, '];
        if strcmp(PL.wgmTaper.endtype,'custom')
            wgmLStr = [wgmLStr,'a_{end}=',num2str(PL.wgmTaper.a_end*1e9,'%.0f'),'nm, ',...
                             'hx_{end}=',num2str(PL.wgmTaper.hx_end*1e9,'%.0f'),'nm, ',...
                             'hy_{end}=',num2str(PL.wgmTaper.hy_end*1e9,'%.0f'),'nm'];
        elseif strcmp(PL.wgmTaper.endtype,'maxdef')
            wgmLStr = [wgmLStr,'d_{max,end}=',num2str(PL.wgmTaper.maxdef_end*1e9,'%.5f'),', ',...
                             'obl_{end}=',num2str(PL.wgmTaper.oblong_end*1e9,'%.4f')];
        end
    else
        wgmLStr = '';
    end
    ttl3 = ['L: n_{h}=',num2str(PL.nholes,'%.0f'),', ',...
        'n_{d}=',num2str(PL.ndef,'%.0f'),', ',...
        'n_{wgm}=',num2str(PL.wvgmir,'%.0f'),...
        ctrholeLStr,wgmLStr];
    TAll = [TAll;ttl3];
end

% Tflds1 = {'theta','a','w','hx','hy'};
% Tflds2 = {'oblong','maxdef','nholes','ndef','wvgmir'};
% if isfield(P,'ctrholeadd') && P.ctrholeadd > 0
%     Tflds2{end+1} = 'ctrholeadd';
% end
% Tline1 = GenFileName(P,Tflds1,'title');
% Tline2 = GenFileName(P,Tflds2,'title');
% TAll = {Tline1;Tline2};
% 
% Tflds3 = '';
% if isfield(P,'wgmTaper') && P.wvgmir > 0
%     Tflds3 = {'wgmTaper.func'};
%     if strcmp(P.wgmTaper.endtype,'custom')
%         Tflds3 = [Tflds3,'wgmTaper.a_end','wgmTaper.hx_end','wgmTaper.hy_end'];
%     elseif strcmp(P.wgmTaper.endtype,'maxdef')
%         Tflds3 = [Tflds3,'wgmTaper.maxdef_end','wgmTaper.oblong_end'];
%     end
% end
% if ~isempty(Tflds3)
%     Tline3 = GenFileName(P,Tflds3,'title');
%     TAll = [TAll;Tline3];
% end
% 
% 
% Tflds4 = '';
% if isfield(P,'taperTo') && strcmp(P.taperTo,'custom')
%     Tflds4 = {'taperFunc','a_ctr','hx_ctr','hy_ctr'};
% end
% if ~isempty(Tflds4)
%     Tline4 = GenFileName(P,Tflds4,'title');
%     TAll = [TAll;Tline4];
% end

% title(bx,TAll,'fontname','arial','fontsize',10)

ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text(0.5, 1,TAll,...
        'HorizontalAlignment' ,'center','VerticalAlignment', 'top',...
        'FontSize',10,'FontWeight','bold')


%% Plot geometry in the yz-plane
if asym > 0
    cx = axes('position',[0.825 0.8 0.125 0.125]);
    hold(cx,'on')
    
    %waveguide overlay
    plot(cx,1/2*[-P.w P.w], 1/2*[P.th P.th],'k','linewidth',1)
    plot(cx,1/2*[-P.w -2*asym], 1/2*[P.th, -P.th],'k','linewidth',1)
    plot(cx,1/2*[-2*asym P.w], 1/2*[-P.th, P.th],'k','linewidth',1)
        
    plot(cx,[0,0],1.05*P.th*[-1,1], ...
        'r','linestyle','--','linewidth',.5);
    plot(cx,[-asym,-asym],1.05*P.th*[-1,1], ...
        'k','linestyle','--','linewidth',.5);
   
    set(cx,'xtick',[],'ytick',[]);
    hold(cx,'off');
    box(cx,'on')
    set(cx,'xcolor','w','ycolor','w');
    daspect([1 1 1])
    xlim(cx,0.65*[-P.w P.w]);
    ylim(cx,[-P.th P.th]);  
    text(0,-(P.th/2),[' ~',num2str(asym*1e9,'%.0f'),'nm'],'fontsize',10)
end





