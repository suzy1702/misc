% Band Tracing of SED data
% Experimental does not work!!

%% New data
h5disp('si.sed.h5');
z = h5read('si.sed.h5','/sed');
halve = log(z(:,1:end/2));
halve = halve';

timestep = 0.0005e-12; % pico seconds
duration = 10000; % number of time steps

tscale = 1/(duration*timestep);
tscale = tscale*10^-12; %scale to terahertz

[nrow, ncol] = size(halve);

trim = 5; %trim off the first couple columns when plotting

surf(halve,'EdgeColor','None');
view(2)

%% Peak tracing algorithm

totpeaksH = {};
totpeaksL = {};
cut = trim:ncol-trim;
for j = cut
    % Standardize the data for every column
    peak = zscore(halve(:,j));
    peakLocal = [];
    peakHeight = [];
    peakMedian = median(peak);
    
    for i=1:(nrow-1)
        %For every point check if its a peak and mark its height and location
            if(peak(i) > peak(i+1) && peak(i) > peak(i-1) && peak(i) > 0)
                peakHeight = [peakHeight , peak(i)];
                peakLocal = [peakLocal, i];
            end
    end
    totpeaksH{end+1} = peakHeight;
    totpeaksL{end+1} = peakLocal;
end

%% trace sketch
scale = max(max(halve));
xvec = cut;
xvec = repelem(xvec,5);
mar = ones(5*50) * (scale+1);

figure;
Y = linspace(1,131072,131072)*tscale;
X = linspace(1,50,50);

surf(X,Y,halve/scale,'EdgeColor','None');
hold on

plotTrace(totpeaksL,cut,scale,tscale);
view(2)
% title('Band tracing for Bulk Silicon','FontSize',24);
% ylabel('Frequency (Terahertz)','FontSize',24);
% xlabel('k Points','FontSize',24);
% set(gca,'FontSize',20)
% legend('SED','Peak')
% names = { 'G'; 'X' ; 'K' ; 'G' ; 'L'};
% set(gca,'xtick',[0,40,80,120,160],'xticklabel',names)


