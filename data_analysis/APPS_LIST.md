# List of applications and their use:

- <plot_signals.cpp signals plots> plots signals obtained in each SiPM. More the input files, more the plots.

- <first_times.cpp initial times analysis> plots and fits arrival times for each SiPM, for each SiPM side, of the difference dx_side - sx_side and plot the position_x versus time points and their fit.

- <charges.cpp charge analysis> plots and fits the total charges in dx and sx side of SiPM for each generated position. We evaluate and fit also some functions of the charges, for example charge_dx/charge_sx. The app do it via EventPerEvent method and via classical estimation.

- <exclusion_band.cpp time versus charge> uses random generated points to produce a scatter plot of arrival time vs charge of the signal for each SiPM side. We expect to see anti-correlation.