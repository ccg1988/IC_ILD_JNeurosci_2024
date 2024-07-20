function [p,r]=fit_pr(scatter_x, scatter_y)
fit_line=   fitlm(scatter_x, scatter_y); 
p=fit_line.Coefficients.pValue(2);
r=sqrt(fit_line.Rsquared.Ordinary);