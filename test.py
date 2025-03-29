import matplotlib.pyplot as plt
import pandas as pd

# Plot 1: Time-domain signal
time_data = pd.read_csv("sine_wave_data.csv")
plt.figure()
plt.plot(time_data['time'], time_data['real'])
plt.title('Time-Domain Signal')
plt.xlabel('Time')
plt.ylabel('Amplitude')
plt.show()

# Plot 2: DFT Real and Imaginary vs. Sample Index
dft_data = pd.read_csv("sine_wave_dft.csv")
plt.figure()
plt.plot(dft_data['index'], dft_data['real'], label='Real')
plt.plot(dft_data['index'], dft_data['imaginary'], label='Imaginary')
plt.title('DFT vs. Sample Index')
plt.xlabel('Sample Index')
plt.ylabel('DFT Value')
plt.legend()
plt.show()

# Plot 3: DFT Magnitude vs. Frequency
mag_data = pd.read_csv("sine_wave_magnitude.csv")
plt.figure()
plt.plot(mag_data['frequency'], mag_data['magnitude'])
plt.title('DFT Magnitude vs. Frequency')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Magnitude')

# Annotate significant peaks
peak_threshold = 0.1 * mag_data['magnitude'].max()
peaks = mag_data[mag_data['magnitude'] > peak_threshold]
for _, peak in peaks.iterrows():
    plt.annotate(f"{peak['frequency']:.3f} Hz",
                 xy=(peak['frequency'], peak['magnitude']),
                 xytext=(0, 10), textcoords='offset points', ha='center')

plt.show()
