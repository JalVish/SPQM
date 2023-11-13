# SPQM functions quantile map projected climate simulations.
Quantile map precipitation
Description: Quantile maps projected precipitation from a climate model simulation at a daily scale
Usage
SPQM_Precipitation_Ft(Observations, Simulations, Threshold)
Arguments
Observations: A dataframe containing reference precipitation data at a daily scale. The dataframe must contain the date and data
Simulations: A dataframe containing projected precipitation data to be quantile mapped. The dataframe must contain the date and data
Threshold: Threshold below which the precipitation is set to zero (considered as a dry day)
References
Rajulapati, C. R., & Papalexiou, S. M. (2023). Precipitation bias correction: A novel semi-parametric quantile mapping method. Earth and Space Science, 10, e2023EA002823. https://doi.org/10.1029/2023EA002823

Quantile map temperature
Description: Quantile maps projected minimum and maximum temperature from a climate model simulation at a daily scale
Usage
SPQM_Temperature_Ft (Observations, Simulations)
Arguments
Observations: A dataframe containing reference minimum and maximum temperature data at a daily scale. The dataframe must contain the date, minimum temperature, and maximum temperature
Simulations: A dataframe containing projected minimum and maximum temperature data to be quantile mapped. The dataframe must contain the date, minimum temperature, and maximum temperature 
References
Rajulapati, C. R., Abdelmoaty, H.M., Nerantzaki, S.D., & Papalexiou, S. M. (2022). Changes in the risk of extreme temperatures in megacities worldwide, Climate Risk Management,  36, 100433 https://doi.org/10.1016/j.crm.2022.100433
Rajulapati, C. R., Gaddam, R.K., Nerantzaki, S.D., Papalexiou, S. M., Cannon, A.J., &  Clark, M.P. (2022). Exacerbated heat in large Canadian cities, Urban Climate,  42, 101097 https://doi.org/10.1016/j.uclim.2022.101097


