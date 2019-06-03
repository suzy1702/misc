% Phonon Lifetime (Riley Hadjis 5/13/2019)
% Each section can be run seperatly and is titled
% Note that the section Dynamic peak fitting can be overided to use
%   a static fitting method 
%   (provide the correct slices of data in peakStartStop)

%% Prep data
load 'SED.mat'
[nrow, ncol] = size(sed);
timestep = 0.0005e-12; % pico seconds
duration = 10000; % Pulled from lammps file
kdist = 70; %kdist to find life times at

tscale = 1/(duration*timestep);
tscale = tscale*10^-12; %scale to terahertz

origSlice = sed(1:end/2,kdist);

% %Plot SED
% surf(sed,'EdgeColor','None');
% view(2)

%% Dynamic Peak finding and prep for fitting
slice = zscore(origSlice);
peakLocal = [];
valleyLocal = [];
for i=2:(length(slice)-1)
    % See if a point is a peak or a valley by checking its neighbors
    if(slice(i-1) < slice(i) && slice(i) > slice(i+1) && slice(i) > 0)
        peakLocal = [peakLocal , i];
    end
    len = length(peakLocal);
    if(slice(i-1) > slice(i) && slice(i) < slice(i+1) && len ~= 0 && length(valleyLocal)~= len)
        valleyLocal = [valleyLocal, i];
    end
end

% Find a symmertic distance around each peak by finding its closest valley
% This will be used to section the data to fit peaks

dist = valleyLocal -  peakLocal;
l = length(peakLocal);
peakStartStop = zeros(2,l);
peakStartStop(:,1) = [peakLocal(2) - dist(1) ; peakLocal(1) + dist(1)];

for i=2:l
    if(dist(i) > dist(i-1))
        shift = peakLocal(i) - (peakLocal(i-1)+dist(i-1));
        peakStartStop(:,i) = [peakLocal(i) - shift ; peakLocal(i) + shift];
    else
        peakStartStop(:,i) = [peakLocal(i) - dist(i) ; peakLocal(i) + dist(i)];
    end   
end

%% Fit peaks and plot
X = linspace(1,250,250)*tscale;
plot(X,log(origSlice),'Color','g')
hold on

for i=1:length(peakLocal)
    peak = peakStartStop(1,i):peakStartStop(2,i);
    [y PARAMS RESNORM RESIDUAL] = lorentzfit(peak,origSlice(peak)');
    disp(PARAMS)
    % MATLAB will output the following display
    % PARAMS = [P1 P2 P3 C]
    % y = P1./((X - P2).^2 + P3) + C
    fit = y(1:length(peak));
    plot(peak*tscale,log(fit),'Color','k')
end

title('Fitted Peaks for a Single k-Point of Bulk Silicon')
xlabel('Frequency (Terahertz)');
ylabel('log(SED)');
legend('SED','Fitted peak')
set(gca,'FontSize',20)
