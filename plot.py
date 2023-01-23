from plot_goesr import make_goes16_ena_daily_movie, make_goes16_ena_daily_plots, download_daily_data

year, month, day = 2018, 2, 12 
download_daily_data(year, month, day)
make_goes16_ena_daily_plots(year, month, day)
make_goes16_ena_daily_movie(year, month, day)
