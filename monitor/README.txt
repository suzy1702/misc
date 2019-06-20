If you dont have it, download and install lm-sensors. Then you might have to edid get_temps to include the correct path to the sensors executable.

By default, cpu temps are recorded every 60 seconds. If you want to change it, edit TIMELAG variable in get_temps. 

Then, run ./get_temps

It will continue to write temps to results.txt until you stop it, i.e. CTRL+C

Then it will call the python script to read and plot the temperatures to show if any cpus have gotten close to the max temp and critical temps.


