import lightkurve as lk 
from scipy.signal import find_peaks
import numpy as np


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


power_values = periodogram.power.value
period_values = periodogram.period.value

peak_indices, _ = find_peaks(power_values, distance = 50)
peak_periods = period_values[peak_indices]
peak_powers = power_values[peak_indices]

sorted_order = peak_powers.argsort()[::-1]

top_n = 5
top_indices = sorted_order[:top_n]
top_periods = peak_periods[top_indices]
top_powers = peak_powers[top_indices]

print(top_periods)
print(top_powers)

for period, power in zip(top_periods, top_powers):
    print(f"Period: {period:.4f} days, Power: {power:.2f}")

print(len(period_values))


num_bins = 15
bin_edges = np.linspace(period_values.min(), period_values.max(), num_bins + 1)

bin_indices = np.digitize(period_values, bin_edges)

bin_winners_period = []
bin_winners_power = []

for bin_num in range(1, num_bins + 1):
    mask = bin_indices == bin_num
    
    if mask.sum() == 0:
        continue
    
    periods_in_bin = period_values[mask]
    powers_in_bin = power_values[mask]
    
    best_local_index = powers_in_bin.argmax()
    
    bin_winners_period.append(periods_in_bin[best_local_index])
    bin_winners_power.append(powers_in_bin[best_local_index])

for p, pw in zip(bin_winners_period, bin_winners_power):
    print(f"Bin winner — Period: {p:.4f} days, Power: {pw:.2f}")


bin_winners_period = np.array(bin_winners_period)
bin_winners_power = np.array(bin_winners_power)

sorted_order = bin_winners_power.argsort()[::-1]
top_5_periods = bin_winners_period[sorted_order[:5]]
top_5_powers = bin_winners_power[sorted_order[:5]]

for p, pw in zip(top_5_periods, top_5_powers):
    print(f"Top candidate — Period: {p:.4f} days, Power: {pw:.2f}")


folded_curves = []

for period in top_5_periods:
    folded = clean_lc.fold(period=period)
    folded_curves.append(folded)

for i, folded in enumerate(folded_curves):
    ax = folded.plot()
    ax.set_title(f"Folded at period {top_5_periods[i]:.4f} days")
    ax.figure.savefig(f"folded_candidate_{i}.png", dpi=150)

''' First binning hard coded
for i, period in enumerate(top_5_periods):
    folded = clean_lc.fold(period=period)
    binned = folded.bin(time_bin_size=0.02)
    
    ax = binned.plot()
    ax.set_title(f"Binned fold at period {period:.4f} days")
    ax.figure.savefig(f"binned_candidate_{i}.png", dpi=150)'''

''' Second bining : period/100 method
for i, period in enumerate(top_5_periods):
    folded = clean_lc.fold(period=period)
    
    bin_width = period / 100
    binned = folded.bin(time_bin_size=bin_width)
    
    ax = binned.plot()
    ax.set_title(f"Binned fold at period {period:.4f} days")
    ax.figure.savefig(f"binned_v2_candidate_{i}.png", dpi=150)'''


zoomed = clean_lc[(clean_lc.time.value > 1325) & (clean_lc.time.value < 1332)]
ax = zoomed.plot()
ax.figure.savefig("zoomed_window.png", dpi=150)

top_5_durations = []
for period in top_5_periods:
    idx = (np.abs(periodogram.period.value - period)).argmin()
    duration = periodogram.duration.value[idx]
    top_5_durations.append(duration)

for p, d in zip(top_5_periods, top_5_durations):
    print(f"Period: {p:.4f} days, Duration: {d:.4f} days")

for i, (period, duration) in enumerate(zip(top_5_periods, top_5_durations)):
    folded = clean_lc.fold(period=period)
    
    bin_width = duration / 5
    binned = folded.bin(time_bin_size=bin_width)
    
    ax = binned.plot()
    ax.set_title(f"Binned fold at period {period:.4f} days (bin={bin_width:.4f}d)")
    ax.figure.savefig(f"binned_v3_candidate_{i}.png", dpi=150)

def get_transit_times(t_start, t_end, period, phase_offset):
    n_transits = int((t_end - t_start) / period)
    transit_times = [t_start + phase_offset + i * period for i in range(n_transits)]
    return transit_times

def get_phase_offset(clean_lc, period):
    folded = clean_lc.fold(period=period)
    binned = folded.bin(time_bin_size=period/100)
    
    min_index = binned.flux.argmin()
    phase_offset = binned.time.value[min_index]
    
    return phase_offset

for period in top_5_periods:
    phase_offset = get_phase_offset(clean_lc, period)
    
    t_start = clean_lc.time.value.min()
    t_end = clean_lc.time.value.max()
    
    transit_times = get_transit_times(t_start, t_end, period, phase_offset)
    
    print(f"Period {period:.4f} days -> {len(transit_times)} predicted transit times")
    print(transit_times)

def measure_transit_depths(clean_lc, transit_times, half_window):
    depths = []
    
    for t_center in transit_times:
        mask = (clean_lc.time.value > t_center - half_window) & (clean_lc.time.value < t_center + half_window)
        
        if mask.sum() == 0:
            depths.append(None)
            continue
        
        local_flux = clean_lc.flux.value[mask]
        depth = 1 - local_flux.mean()
        depths.append(depth)
    
    return depths

top_5_depths = []

for period, duration in zip(top_5_periods, top_5_durations):
    phase_offset = get_phase_offset(clean_lc, period)
    t_start = clean_lc.time.value.min()
    t_end = clean_lc.time.value.max()
    
    transit_times = get_transit_times(t_start, t_end, period, phase_offset)
    
    half_window = duration / 2
    depths = measure_transit_depths(clean_lc, transit_times, half_window)
    
    top_5_depths.append(depths)
    
    print(f"\nPeriod {period:.4f} days:")
    print(f"Depths per occurrence: {depths}")

def compute_consistency_score(depths):
    clean_depths = [d for d in depths if d is not None]
    
    if len(clean_depths) < 2:
        return 0
    
    clean_depths = np.array(clean_depths)
    mean_depth = clean_depths.mean()
    std_depth = clean_depths.std()
    
    if std_depth == 0:
        return 0
    
    score = mean_depth / std_depth
    return score

for period, depths in zip(top_5_periods, top_5_depths):
    score = compute_consistency_score(depths)
    print(f"Period: {period:.4f} days, Consistency score: {score:.2f}")


