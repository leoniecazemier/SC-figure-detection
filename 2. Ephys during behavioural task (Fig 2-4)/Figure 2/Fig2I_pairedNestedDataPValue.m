function [p_value, F, R] = pairedNestedDataPValue(y, condition, mouse, session, unit) 
% y: The response variable (dependent variable) 
% group: The grouping variable indicating the different levels of nesting 
% subject: The subject identifier variable indicating the nesting within groups 

data = table(y, condition, unit, session, mouse); 
lme_model = fitlme(data, 'y ~ condition + (1|mouse) + (1|session) + (1|unit)', 'FitMethod','REML','StartMethod','random'); 
[p_value,F,R] = coefTest(lme_model); end



