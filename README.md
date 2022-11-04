# ensemble_n-subepidemic_framework
A Matlab toolbox for fitting and forecasting epidemic trajectories using the ensemble n-subepidemic framework

It carries out the following tasks:

    fitting models to time series data,
    estimation of parameters with quantified uncertainty,
    plotting the fits of the top-ranked models,
    plotting the AICc values of the top-ranked models,
    generates forecasts of the top-ranked models
    gnerates ensemble forecasts based on the top-ranked models.

Additional features include:

    fitting models using different parameter estimation approaches (least-squares, maximum likelihood estimation),
    fitting models using assuming different error structures (normal, Poisson, negagive binomial),
    user can select the underlying function for the sub-epidemic building block (generalized-logistic model, Richards model, Gompertz model),
    user can select whether the sub-epidemics start synchronously at time 0 or asynchronously at different times as defined by parameter C_thr.
    
  
