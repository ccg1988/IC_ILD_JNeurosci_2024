function [xi,Fitb]=fit_scatter(scatter_x, scatter_y, x_range)
Fita=polyfit(scatter_x, scatter_y, 1);%linear fit
if strcmp(x_range, 'free')
xi=min(scatter_x):0.005:max(scatter_x);%the range of SC values
else
xi=0.1:0.005:0.55;
end
Fitb=polyval(Fita,xi);