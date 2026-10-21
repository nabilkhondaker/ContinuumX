"""Simple frequency estimation from tip history."""
import numpy as np

def dominant_frequency(signal, dt):
    if len(signal) < 8:
        return 0.0
    spec = np.abs(np.fft.rfft(signal - np.mean(signal)))
    freqs = np.fft.rfftfreq(len(signal), dt)
    return float(freqs[np.argmax(spec[1:]) + 1])
