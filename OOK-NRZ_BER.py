import math
import numpy as np
import matplotlib.pyplot as plt

# Calculate the Bit Error Rate (BER) for an On-Off Keying (OOK) None Return to Zero (NRZ) optical communication system with a given SNR.
def calculate_ber_ook(snr_db):
    snr_linear = 10 ** (snr_db / 10.0)  # Convert SNR from dB to linear scale
    ber = 0.5 * math.erfc(math.sqrt(snr_linear))
    return ber

# Range of SNR values in dB
snr_db_range = np.arange(0, 13, 0.5)

# Calculate BER for each SNR value
ber_ook = [calculate_ber_ook(snr_db) for snr_db in snr_db_range]

# Plot the BER vs SNR using a semilog scale
plt.semilogy(snr_db_range, ber_ook, marker='o')
plt.xlabel('SNR (dB)')
plt.ylabel('Bit Error Rate (BER)')
plt.title('OOK Optical Communication System')

plt.grid()
plt.show()