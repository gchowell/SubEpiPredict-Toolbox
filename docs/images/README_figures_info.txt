Figures made by make_readme_figures.m on 2026-11-22 07:49
MATLAB 24.1.0.2689473 (R2024a) Update 6 | fit revision m1 | random seed 1
Data: coronavirus deaths USA (daily), area (column) 52, data to 05-11-2020, time step DT = 1
Fit: method1 = 3, dist1 = 3, flag1 = 1, up to 2 sub-epidemic(s), onset_fixed = 0, smoothfactor1 = 7, calibrationperiod1 = 90, numstartpoints = 20, topmodelsx = 4, B = 40
Forecast: horizon 30, weight_type1 = -1

model_fit.png          figure 305  4 ranked models: fit with bootstrap simulations, sub-epidemic components, residuals
parameters.png         figure 301  bootstrap parameter distributions of ranked model 1
rankings.png           figure 99   AICc, relative likelihood and evidence ratio of the ranked models
ranking_profiles.png   figure 100  sub-epidemic profiles of the ranked models
ranked_forecasts.png   figure 10   forecasts of the ranked models, 30 time units ahead
ensembles.png          figure 500  ensemble forecasts, 30 time units ahead, weight_type1 = -1
