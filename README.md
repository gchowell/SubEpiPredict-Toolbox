# ensemble_n-subepidemic_framework
A MATLAB toolbox for fitting and forecasting epidemic trajectories using the ensemble n-subepidemic framework

<p> It carries out the following tasks: </p> 
<ul>
    <li>fitting models to time series data,</li>
    <li>estimation of parameters with quantified uncertainty,</li>
    <li>plotting the model fits, sub-epidemic profiles, and residuals of the top-ranked models,</li>
    <li>plotting the empirical distributions of the model parameters associated with each sub-epidemic</li>
    <li>plotting the calibration performance metrics of the top-ranked models, </li>    
    <li>plotting the AICc values, relative likelihood, and evidence ratio of the top-ranked models,</li>
    <li>plotting forecasts derived from the top-ranked and ensemble models,</li>
    <li>plotting the forecasting performance metrics of the top-ranked models and the ensemble models, </li>    
    <li>plotting the effective reproduction number derived from the top-ranked models.</li>
    
</ul>

<p> Additional features include:</p>

<ul>
    <li>fitting models using different parameter estimation approaches (least-squares, maximum likelihood estimation),</li>
    <li>fitting models using assuming different error structures (normal, Poisson, negagive binomial),</li>
    <li>user can select the underlying function for the sub-epidemic building block (generalized-logistic model (GLM), Richards model, generalized Richards model (GRM), Gompertz model),</li>
    <li>user can select whether the sub-epidemics start synchronously at time 0 or asynchronously as defined by parameter C_thr.</li>
    
</ul>
    
# Installation requirements

The n-subepidemic framework toolbox requires a MATLAB installation.

# Fitting the model to your data

To use the toolbox to fit the ensemble n-subepidemic framework to your data, you just need to:

<ul>
    <li>download the code </li>
    <li>create 'input' folder in your working directory where your data is located </li>
    <li>create 'output' folder in your working directory where the output files will be stored</li>   
    <li>open a MATLAB session </li>
    <li>define the model parameter values and time series parameters by editing <code>options.m</code> </li>
    <li>run the function <code>Run_Fit_subepidemicFramework.m</code> </li>
</ul>
  
# Plotting the fits of the top-ranked models and parameter estimates

After fitting the model to your data, you can use the toolbox to plot the model fits and parameter estimates as follows:

<ul>
    <li>run the function <code>plotFit_subepidemicFramework.m</code> </li>
</ul>
    
The function also outputs a file with the calibration performance metrics of the top-ranked models.
    
# Plotting the top-ranked subepidemic model profiles and the corresponding AIC values

After fitting the model to your data, you can use the toolbox to plot the subepidemic profiles and AICc values as follows:

<ul>
    <li>run the function <code>plotRankings_subepidemicFramework.m</code></li>
</ul>
    
# Generating and plotting forecasts of the top-ranked and ensemble subepidemic models

After fitting the model to your data, you can use the toolbox to plot forecasts derived from the top-ranked and ensemble subepidemic models as follows:

<ul>
    <li>define the forecasting parameters by editing <code>options_forecast.m</code></li>
    <li>run the function <code>plotForecast_subepidemicFramework.m</code></li>
</ul>

The function also outputs files with the fit and forecasts of the top-ranked and ensemble models as well as the forecasting performance metrics for the top-ranked and ensemble models.

# Generating and plotting reproduction number forecasts from the top-ranked models

After generating forecasts from top-ranked models, you can use the toolbox to generate and plot forecasts of the effective reproduction number from the top-ranked models as follows:

<ul>
    <li>define the generation interval parameters by editing the function <code>options_rt.m</code></li>
    <li>run the function <code>plotReproductionNumber.m</code></li>
</ul>

# Publications

<ul>

 <li>Chowell, G., Dahal, S., Tariq, A., Roosa, K., Hyman, J. M., & Luo, R. (2022). An ensemble n-sub-epidemic modeling framework for short-term forecasting epidemic trajectories: Application to the COVID-19 pandemic in the USA. PLOS Computational Biology, 18(10), e1010602. </li>

 <li>Chowell, G., Tariq, A., & Hyman, J. M. (2019). A novel sub-epidemic modeling framework for short-term forecasting epidemic waves. BMC medicine, 17(1), 1-18.  </li>
 
<li>Bleichrodt, A., Dahal, S., Maloney, K., Casanova, L., Luo, R., & Chowell, G. (2022). Real-time forecasting the trajectory of monkeypox outbreaks at the national and global levels, July–October 2022. BMC Medicine (In Press). </li>

</ul>

# Disclaimer

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, version 3 of the License.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License for more details.  
