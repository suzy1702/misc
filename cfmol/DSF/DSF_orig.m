clear

a = 5.4310;  Nxy = 8;  Lxy = a*Nxy; % a is the lattice parameter and Nxy is the number of unit cells in the x and y directions.
CS = 'interf'; % CS = 'all','Si','Siall','interf','SiO';
timestep = 0.00016; %0.00016;0.0008 % in unit of ps, MD timestep
%Nstep = 21000;
filenm = 'Si_bulk'; %'SiO3nmV0.2_sp1';

fid = fopen(['position_',filenm],'r'); % open the lammps position file to read the number of atoms, number of atom types and the original atom positions.
aa = textscan(fid,'%d%s',1,'headerlines',1);
Nat = aa{1};
aa = textscan(fid,'%d%s%s',1,'headerlines',2);
Ntype = aa{1};
aa = textscan(fid,'%f%f%f%f%f','headerlines',10+Ntype);
position = [aa{3} aa{4} aa{5}];
fclose(fid); 
if Nat ~= size(position,1)
    error('atom number mismatch');
end


filename = strcat(filenm,'.pos.dat'); % the position file contains atom ID TYPE x y z
q = zeros(3,Nxy+1); 
q(1,:) = [(0:Nxy)/Nxy]; % set the k points (or q points)
qL = q; qT1 = q; qT2 = q;
qL = 2*pi*qL/a; qT1(2,:) = 2; qT1 = 2*pi*qT1/a;
qT2(3,:) = 2; qT2 = 2*pi*qT2/a;
qall = [qL, qT1, qT2];



fid = fopen(filename,'r');
aa = textscan(fid,'%f',1,'headerlines',1);
step1 = aa{1};
aa = textscan(fid,'%f',1,'headerlines',2);
Natom = aa{1};
frewind(fid);
aa = textscan(fid,'%f',1,'headerlines',Natom+10);
frewind(fid);
step2 = aa{1};
dt = (step2 - step1)*timestep

if ~exist('Nstep','var')
    [status,Nline] = system(['wc -l ', filename,'  | awk ''{print $1}'' ']); % use the linux or Mac command to get the number of lines for the trajectory file
    Nline = str2num(Nline);
    Nstep = floor(double(Nline)/(Natom+9)); % calculate the number of snapshots
    Nstep = floor(Nstep/2)*2  % make Nstep an even value
end

QdotR = zeros(Nstep,3*(Nxy+1));
ID = 1:Nat;
for i = 1:Nstep
    if i == 1
        aa = textscan(fid,'%d%d%f%f%f',Natom,'headerlines',9);
        pos = [aa{3} aa{4} aa{5}];
        pos = pos(ID,:);
        pos_ori = pos;

        rho = exp(-1i*pos*qall);  % space Fourier transform
        QdotR(i,:) = sum(rho,1);
    else
        aa = textscan(fid,'%d%d%f%f%f',Natom,'headerlines',10);
        pos = [aa{3} aa{4} aa{5}];
        pos = pos(ID,:);
        dpos = pos - pos_ori;
        IDb = find(dpos > Lxy/2);  % elimilate the cross boundary effect in the periodic directions! 
        pos(IDb) = pos(IDb) - Lxy; % this works only for Lx = Ly;
        IDb = find(dpos < -Lxy/2);
        pos(IDb) = pos(IDb) + Lxy;

        rho = exp(-1i*pos*qall);
        QdotR(i,:) = sum(rho,1);
    end

end
fclose(fid);


%%%%%%%%%%%%% do time Fourier transform %%%%%%%%%%
n = (0:Nstep/2-1);
fre = n/(dt*Nstep);
df = 1/(dt*Nstep);
Sfactor = zeros(Nstep/2,3*(Nxy+1));

for i = 1:3*(Nxy+1)
    Sfft = abs(fft(QdotR(:,i),Nstep));
    Sfactor(:,i) = Sfft(1:Nstep/2);
end
%%%%%%%%%%%%% end of time Fourier transform %%%%%%%%%%%

SL = Sfactor(:,1:Nxy+1); SL = [fre',SL];
ST1 = Sfactor(:,Nxy+2:2*(Nxy+1)); ST1 = [fre',ST1];
ST2 = Sfactor(:,2*(Nxy+1)+1:end); ST2 = [fre',ST2];

fid = fopen(['SFactor_',filenm,'_',CS,'.dat'],'w');

format = repmat('%15.6f ',[1,Nxy+2]);
format = [format,'\n'];
fprintf(fid,'Longitudinal mode at q=(0:%d)/%d\n',Nxy,Nxy);
fprintf(fid, format, SL');

fprintf(fid,'Transverse mode 1 at qL + (0,2,0)\n');
fprintf(fid, format, ST1');

fprintf(fid,'Transverse mode 2 at qL + (0,0,2)\n');
fprintf(fid, format, ST2');
fclose(fid);
