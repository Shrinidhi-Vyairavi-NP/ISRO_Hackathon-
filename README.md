# ISRO_Hackathon-
Detect Exoplanet transit signals from noisy astronomical light curves 

Progress so far : 
A transit occurs when an exoplanet passes between a star and its observer - This causes a dip in the lightcurve (star's brightness is hindered) 
Learn more : https://science.nasa.gov/exoplanets/whats-a-transit/

AIM : To detect this exoplanet by finding the dip in the lightcurve, identifying the dip and type of exoplanet 

LIGHTCURVE : brightness vs time curve 

METHOD USED SO FAR : Test run with known exoplanets and train the AI Model to identify 
Test star : Pi Mensae 
Test exoplanet : Pi Mensae C Orbital Period : 6.3 days 
Learn more : https://science.nasa.gov/exoplanet-catalog/pi-mensae-c/

Using built in python libraries like lightkurve we clean the noisy data and identify the dips by using the Box Least Squared Algorithm. 
Train a classifier to then take other parameters like dip depth, length, etc. to actually identify the exoplanet [NOT REACHED THERE YET] 

Notable/Important Methods : 
.plot()
.flatten()
.remove_outliers(sigma = 5)
.to_periodogram(method = "bls")
