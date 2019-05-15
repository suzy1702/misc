%% Prep data
load 'SED.mat'
timestep = 0.0005e-12; % pico seconds
duration = 10000; % number of time steps 

tscale = 1/(duration*timestep);
tscale = tscale*10^-12; %scale to terahertz

halve = zscore(sed(1:100,1:120)); %take a section of the data
[nrow, ncol] = size(halve);

trim = 5; %trim off the first couple columns when plotting

%% Peak tracing algorithm

totpeaksH = {};
totpeaksL = {};

cut = trim:ncol-trim;
for j = cut
    % for every column
    peak = halve(:,j);
    peakLocal = [];
    peakHeight = [];
    peakMedian = median(peak);
    
    for i=5:(nrow-1)
        %Check if a point is larger than its neighbors and a threshold
            if(peak(i) > peak(i+1) && peak(i) > peak(i-1) && peak(i) > 0)
                peakHeight = [peakHeight , peak(i)];
                peakLocal = [peakLocal, i];
            end
    end
    totpeaksH{end+1} = peakHeight;
    totpeaksL{end+1} = peakLocal;
end

%% Plot points
scale = max(max(halve));
xvec = cut;
xvec = repelem(xvec,5);
mar = ones(5*50) * (scale+1);

X = linspace(1,nrow,nrow)*tscale;
Y = linspace(1,ncol,ncol);
surf(Y,X,halve/scale,'EdgeColor','None');
view(2)
hold on

for i=1:(length(totpeaksL))
    %For every peak in the cell, plot the peaks above the 3d plot
    tmp = cell2mat(totpeaksL(i))*tscale;
    mar = (scale+1).*ones(length(tmp));
    xvec = (cut(i)).*ones(length(tmp));
    plot3(xvec,tmp,mar,'r*','MarkerSize',10,'LineWidth',1);
end


title('Band tracing for Bulk Silicon','FontSize',24);
ylabel('Frequency (Terahertz)','FontSize',24);
xlabel('k Points','FontSize',24);
set(gca,'FontSize',20)
legend('SED','Peak')
names = { 'G'; 'X' ; 'K' ; 'G' ; 'L'};
set(gca,'xtick',[0,40,80,120,160],'xticklabel',names)
