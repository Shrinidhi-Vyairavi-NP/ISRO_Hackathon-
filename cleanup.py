import lightkurve as lk 

star_name = "Pi Mensae" #name of star to analyze 
search_result = lk.search_lightcurve(star_name, mission = 'TESS', author = 'SPOC') #search for lightcurve data from TESS mission
print(search_result) #print the search results 

lc = search_result[0].download() #download the first lightcurve data - search_result is a list
print(lc) #print the lightcurve data 

ax = lc.plot() #plot the lightcurve data 
ax.figure.savefig("trial_raw_lightcurve.png", dpi = 150) #save the plot as a png with 150 res

clean_lc = lc.flatten() #flatten the lightcurve data to remove those weird outliers like extreme peaks
ax_clean = clean_lc.plot() #plot the cleaned lightcurve data 
ax_clean.figure.savefig("trial_clean_lightcurve.png", dpi = 150) #save plot 

clean_lc = clean_lc.remove_outliers(sigma = 5) #sigma clipping extremes aka noise the stats logic we learnt 
ax2_clean = clean_lc.plot()
ax2_clean.figure.savefig("trial_clipped_lightcurve.png", dpi = 150)

periodogram = clean_lc.to_periodogram(method = "bls") #box least squared method 
print(periodogram.period.min(), periodogram.period.max())

best_period = periodogram.period_at_max_power
print(f"Best-fit period found: {best_period}")

ax3 = periodogram.plot()
ax3.figure.savefig("trial_bls_lightcurve.png", dpi = 150)