function [] = plotTrace(peakDat,cut,scale,tscale)
%findTrace given dynamic peak finding data, plot individual traces
level = 0.8;
bandN = 1;

l = length(peakDat);
size(peakDat)
yep = cell2mat(peakDat(1));
yep = yep(bandN)*tscale;
trace1 = [yep];
for i=2:l
    % find near points and interpolate if there is a hole
    yep = cell2mat(peakDat(i));
    yep = sort(yep);
    yep = yep(bandN)*tscale;
%     avg = (yep + trace1(i-1))/2;
    avg = yep;
    trace1 = [trace1 avg];
    
    if(abs(trace1(i)-trace1(i-1)) > level)
        disp('hole found')
        trace1(i) = (trace1(i-1)/cut(i-1))*cut(i);
    end
end

avg = [trace1(1)];
for i = 2:l-1
    val = trace1(i)+trace1(i-1)+trace1(i+1);
    avg = [avg val/3];
end
avg = [avg trace1(end)];
trace1 = avg;

tmp = (scale+1).*ones(length(trace1));
plot3(cut,trace1,tmp,'LineWidth',5)
end

