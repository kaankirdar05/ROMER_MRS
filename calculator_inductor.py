#Calculates the inductance from the given voltage across the inductor (According to the simple RL high-pass filter)
import numpy as np

Vout = 0.21/2  # Voltage across the inductor
Vin = 4/2      # Input voltage amplitude
#Vl = np.sqrt(Vin**2 - Vout**2)
R = 2200    #Do not forget to add the internal resistance of the signal generator (usually 50 ohms)
f = 500 * 1000 #Frequency of the signal generator
w = 2 * np.pi * f #Angular frequency

L = (Vout * R) / (w * np.sqrt(Vin**2 - Vout**2))
#Calculate the inductance
print(L * 1e6)  # μH    #Print the inductance in microhenries