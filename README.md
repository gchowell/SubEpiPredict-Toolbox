# SubEpiPredict Toolbox

**SubEpiPredict** is a user-friendly **MATLAB toolbox** designed for fitting and forecasting epidemic trajectories using the **ensemble n-subepidemic modeling framework**. This approach is particularly effective for capturing complex epidemic patterns, including multiple waves and overlapping sub-epidemics.

📄 **Tutorial Paper**  
Chowell et al. (2024), *SubEpiPredict: A tutorial-based primer and toolbox for fitting and forecasting growth trajectories using the ensemble n-sub-epidemic modeling framework*, Infectious Disease Modelling.  
👉 [Read the paper](https://www.sciencedirect.com/science/article/pii/S2468042724000125)

🎥 **Video Tutorial**  
👉 [Watch on YouTube](https://www.youtube.com/watch?v=lj_-2Kre1qw)

---

## Features

The toolbox offers the following capabilities:

- **Fitting models to time series data**
- **Estimating parameters with quantified uncertainty**
- **Plotting model fits, sub-epidemic profiles, and residuals of top-ranked models**
- **Visualizing empirical distributions of model parameters for each sub-epidemic**
- **Assessing calibration performance metrics of top-ranked models**
- **Analyzing AICc values, relative likelihoods, and evidence ratios**
- **Generating forecasts from top-ranked and ensemble models**
- **Evaluating forecasting performance metrics**
- **Estimating and plotting the effective reproduction number (Rt) from top-ranked models**

**Additional features include:**

- **Support for different parameter estimation approaches (least squares, maximum likelihood estimation)**
- **Flexibility in error structures (normal, Poisson, negative binomial)**
- **Choice of sub-epidemic building block functions:**
  - **Generalized Logistic Model (GLM)**
  - **Richards Model**
  - **Generalized Richards Model (GRM)**
- **Option to model sub-epidemics starting synchronously at time 0 or asynchronously using parameter `C_thr`**

---
    
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

## Publications

- Chowell, G., Dahal, S., Bleichrodt, A., Tariq, A., Hyman, J. M., & Luo, R. (2024). *SubEpiPredict: A tutorial-based primer and toolbox for fitting and forecasting growth trajectories using the ensemble n-sub-epidemic modeling framework*. Infectious Disease Modelling, 9(2), 411-436.  
  [https://www.sciencedirect.com/science/article/pii/S2468042724000125](https://www.sciencedirect.com/science/article/pii/S2468042724000125)
- Chowell, G., Dahal, S., Tariq, A., Roosa, K., Hyman, J. M., & Luo, R. (2022). *An ensemble n-sub-epidemic modeling framework for short-term forecasting epidemic trajectories: Application to the COVID-19 pandemic in the USA*. PLOS Comput Biol, 18(10), e1010602.
- Chowell, G., Tariq, A., & Hyman, J. M. (2019). *A novel sub-epidemic modeling framework for short-term forecasting epidemic waves*. BMC Medicine, 17(1), 164.
- Bleichrodt, A., Dahal, S., Maloney, K., Casanova, L., Luo, R., & Chowell, G. (2022). *Real-time forecasting the trajectory of monkeypox outbreaks at the national and global levels, July–October 2022*. BMC Medicine, 21(1), 1–20.

---

## License

This project is licensed under the terms of the **GNU General Public License v3.0**.  
See the [LICENSE](LICENSE) file for more information.

---

## Contact

For questions or feedback, please contact:  
**Gerardo Chowell**  
[https://github.com/gchowell](https://github.com/gchowell)
