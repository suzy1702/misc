%%%%%%% plot structure factor at different q point
clear
Nxy = 8;  CL = 14:16; CT1 = 1:10; CT2 = 1:10;
filenm = 'SFactor_Si_bulk_all.dat';
%filenm = 'SFactor_Si_bulk_interf.dat';
%filenm = 'SFactor_SiO3nmV0.2_sp1_interf.dat';

fid = fopen(filenm,'r');
[status,Nline] = system(['wc -l ', filenm,'  | awk ''{print $1}'' ']);
Nline = str2num(Nline);
Nread = Nline/3;
if Nread ~= round(Nline/3)
    error('Nline in structure file mistake');
end

format = repmat('%f',[1,Nxy+2]);
aa = textscan(fid,format,Nread,'headerlines',1);
SL = [aa{:}];

aa = textscan(fid,format,Nread,'headerlines',1);
ST1 = [aa{:}];

aa = textscan(fid,format,Nread,'headerlines',1);
ST2 = [aa{:}];

%SL(1:2,:) = [];
%ST1(1:2,:) = [];
%ST2(1:2,:) = [];

%%%%%%%%%%%%%%% make a convolution to smooth the curves %%%%%%%%%%%%%%%
dom = SL(2,1) - SL(1,1);
dwin = 14.01; % THz
win = round(dwin/dom);
g = gausswin_my(win);
g = g/sum(g);
for i = 1:Nxy+1
    SL(:,i+1) = conv(SL(:,i+1),g,'same');
    ST1(:,i+1) = conv(ST1(:,i+1),g,'same');
    ST2(:,i+1) = conv(ST2(:,i+1),g,'same');
end




x = (1:Nxy)/Nxy;
yL = SL(:,1); yT1 = ST1(:,1); yT2 = ST2(:,1);

SL(:,1:2) = [];  %SL = SL/max(SL(:));
ST1(:,1:2) = []; %ST1 = ST1/max(ST1(:));
ST2(:,1:2) = []; %ST2 = ST2/max(ST2(:));

[XL,YL] = meshgrid(x,yL);
[XT1,YT1] = meshgrid(x,yT1);
[XT2,YT2] = meshgrid(x,yT2);
figure
subplot(1,4,1)
pcolor(XL,YL,SL); shading interp
subplot(1,4,2)
pcolor(XT1,YT1,ST1); shading interp
subplot(1,4,3)
pcolor(XT2,YT2,ST2); shading interp

%ylim([0 14])
