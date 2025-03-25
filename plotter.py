import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Load the data
data = pd.read_csv("data.csv")
dft = pd.read_csv("dft.csv")

# 1. Plot the time-domain signal
plt.figure(figsize=(10, 5))  # Adjust figure size as needed
plt.plot(data["time"], data["real"])
plt.title("Sine Wave Signal")
plt.xlabel("Time")
plt.ylabel("Amplitude")
plt.grid(True)  # Add grid lines for better readability
plt.show()

# 2. Plot the DFT (Magnitude)
plt.figure(figsize=(10, 5))
#   Calculate magnitude of DFT
dft['magnitude'] = dft['real'] + dft['imaginary']
plt.plot(dft["index"], dft["magnitude"])
plt.title("DFT Magnitude")
plt.xlabel("Frequency Index")
plt.ylabel("Magnitude")
plt.grid(True)
plt.show()