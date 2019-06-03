%%%%%%% plot structure factor at different q point
clear
Nxy = 16;  
qsampling = 1;
%CL = 14:16; CT1 = 1:10; CT2 = 1:10;
%filenm1 = 'SFactor_Si_bulk_all.dat';
%filenm1 = 'SFactor_Si_bulk_all_16x16x16_16qpts.dat'
filenm1 = 'SFactor_Si_8x8x8_int_Ge_64_atoms_grid_all.dat'
% filenm1 = 'SFactor_Si3nm_all.dat';
% filenm2 = 'SFactor_SiO3nm_all.dat';
% filenm3 = 'SFactor_SiO3nm_interf.dat';
% filenm4 = 'SFactor_SiO3nm_Si.dat';
% %filenm5 = 'SFactor_SiO3nmV0.2_sp1_interf.dat';
% filenm5 = 'SFactor_SiO3nmNh10_interf2.dat';
N=1;
% for N = 1:5
    eval(['filenm = filenm',num2str(N),';']);
fid = fopen(filenm,'r');
[status,Nline] = system(['wc -l ', filenm,'  | awk ''{print $1}'' ']);
Nline = str2num(Nline)
Nread = Nline/3;
if Nread ~= round(Nline/3)
    error('Nline in structure file mistake');
end

format = repmat('%f',[1,qsampling*Nxy+2]);
aa = textscan(fid,format,Nread,'headerlines',1);
SL = [aa{:}];

aa = textscan(fid,format,Nread,'headerlines',1);
ST1 = [aa{:}];

aa = textscan(fid,format,Nread,'headerlines',1);
ST2 = [aa{:}];

SL(1:3,:) = []; ST1(1:3,:) = []; ST2(1:3,:) = [];

%%%%%%%%%%%%%%% make a convolution to smooth the curves %%%%%%%%%%%%%%%
dom = SL(2,1) - SL(1,1)
dwin = 0.125 %0.25 % THz
win = round(dwin/dom)
g = gausswin_my(win)
g = g/sum(g)
for i = 1:Nxy+1
    SL(:,i+1) = conv(SL(:,i+1),g,'same');
    ST1(:,i+1) = conv(ST1(:,i+1),g,'same');
    ST2(:,i+1) = conv(ST2(:,i+1),g,'same');
end


ID = find(SL(:,1) < 20);
SL = SL(ID,:);
ID = find(ST1(:,1) < 20);
ST1 = ST1(ID,:);
ID = find(ST2(:,1) < 20);
ST2 = ST1(ID,:);
x = (1:qsampling*Nxy)/(qsampling*Nxy);
yL = SL(:,1); yT1 = ST1(:,1); yT2 = ST2(:,1);


SL(:,1:2) = [];  SL = SL/max(SL(:));
ST1(:,1:2) = []; ST1 = ST1/max(ST1(:));
ST2(:,1:2) = []; ST2 = ST2/max(ST2(:));
eval(['SL',num2str(N),'=SL;']);
eval(['ST',num2str(N),'1=ST1;']);
eval(['ST',num2str(N),'2=ST2;']);
%end

[XL,YL] = meshgrid(x,yL);

SF1 = SL1 + 0.5*ST1+0.5*ST2;
%SF1 = SL
%SF1 = SL1 + 0.5*ST11+0.5*ST12;
% SF2 = SL2 + 0.5*ST21+0.5*ST22;
% SF3 = SL3 + 0.5*ST31+0.5*ST32;
% SF4 = SL4 + 0.5*ST41+0.5*ST42;
% SF5 = SL5 + 0.5*ST51+0.5*ST52;

figure % n,m,p and position can not be used at the same time in subplot, 
      % that is: subplot(n,m,p) and set(gca,'Position',[x,y,w,h]) can not
      % be used at the same time.
% subplot('Position',[0.115 0.2 0.16 0.75])
pcolor(XL,YL,SF1); shading interp
%set(gca, 'clim', [0.14 1]);
%colormap([0 0 0; jet]);
%colorbar;
ylim([0 17]); xlim([0.0625 1]); box on;
ylabel('Frequency (THz)','FontSize',30);
set(gca,'FontSize',20,'XTick',[0.0625 1],'XTickLabel',{'[000]','[100]'});
set(gca,'YTick',[0.3 2 4 6 8 10 12 14 16],'YTickLabel',{'0' '2' '4' '6' '8' '10' '12' '14' '16'});
xlabel('q (pi/a)','FontSize',30);
%text(0.08,13.3,'(a)','FontSize',24,'BackGroundColor','w','FontWeight','Bold');

% subplot('Position',[0.293 0.2 0.16 0.75])
% pcolor(XL,YL,SF2); shading interp
% ylim([0.3 14]); box on;
% set(gca,'FontSize',20,'YTick',[]);
% set(gca,'XTick',[0.1 0.97],'XTickLabel',{'0','1'});
% text(0.08,13.3,'(b)','FontSize',24,'BackGroundColor','w','FontWeight','Bold');
% 
% subplot('Position',[0.471 0.2 0.16 0.75])
% pcolor(XL,YL,SF3); shading interp
% ylim([0.3 14]); box on;
% xlabel('q (\pi/a)','FontSize',30);
% set(gca,'FontSize',20,'YTick',[]);
% set(gca,'XTick',[0.1 0.97],'XTickLabel',{'0','1'});
% text(0.08,13.3,'(c)','FontSize',24,'BackGroundColor','w','FontWeight','Bold');
% 
% subplot('Position',[0.649 0.2 0.16 0.75])
% pcolor(XL,YL,SF4); shading interp
% ylim([0.3 14]); box on;
% set(gca,'FontSize',20,'YTick',[]);
% set(gca,'XTick',[0.1 0.97],'XTickLabel',{'0','1'});
% text(0.08,13.3,'(d)','FontSize',24,'BackGroundColor','w','FontWeight','Bold');
% 
% subplot('Position',[0.827 0.2 0.16 0.75])
% pcolor(XL,YL,SF5); shading interp
% ylim([0.3 14]); box on;
% set(gca,'FontSize',20,'YTick',[]);
% set(gca,'XTick',[0.1 0.97],'XTickLabel',{'0','1'});
% text(0.08,13.3,'(e)','FontSize',24,'BackGroundColor','w','FontWeight','Bold');
