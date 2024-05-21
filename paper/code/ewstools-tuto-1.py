import time
import numpy as np
import matplotlib.pyplot as plt
import ewstools
from ewstools.models import simulate_ricker

def ewstools_setup():
    # Initialize time series and spectrum computation
    series = simulate_ricker(tmax=1000, F=[0,2.7])
    ts = ewstools.TimeSeries(data=series, transition=860)
    ts.detrend(method='Lowess', span=0.2)
    ts.state[['state','smoothing']].plot()
    ts.compute_spectrum(rolling_window=0.5, ham_length=40)
    return ts

np.random.seed(0)           # Set random seed for reproducibility
rw = 0.5                    # Rolling window for variance computation
ts = ewstools_setup()       # Setup time series
n = 100                     # Number of iterations for timing
t_elapsed = np.zeros(n)     # Array to store elapsed times
t_minruntime = np.zeros(8) # Array to store minimum runtime

for i in range(n):
    t0 = time.time()
    ts.compute_var(rolling_window=rw)
    t_elapsed[i] = time.time() - t0
t_minruntime[0] = min(t_elapsed)

for i in range(n):
    t0 = time.time()
    ts.compute_cv()
    t_elapsed[i] = time.time() - t0
t_minruntime[1] = min(t_elapsed)

for i in range(n):
    t0 = time.time()
    ts.compute_skew(rolling_window=rw)
    t_elapsed[i] = time.time() - t0
t_minruntime[2] = min(t_elapsed)

for i in range(n):
    t0 = time.time()
    ts.compute_kurt()
    t_elapsed[i] = time.time() - t0
t_minruntime[3] = min(t_elapsed)

for i in range(n):
    t0 = time.time()
    ts.compute_auto(rolling_window=rw, lag=1)
    t_elapsed[i] = time.time() - t0
t_minruntime[4] = min(t_elapsed)

for i in range(n):
    t0 = time.time()
    ts.compute_smax()
    t_elapsed[i] = time.time() - t0
t_minruntime[5] = min(t_elapsed)

for i in range(n):
    t0 = time.time()
    ts.compute_ktau()
    t_elapsed[i] = time.time() - t0
t_minruntime[6] = min(t_elapsed)

for i in range(n):
    t0 = time.time()
    ewstools.core.block_bootstrap(ts.state.residuals, 1, bs_type='Stationary', block_size=10)
    t_elapsed[i] = time.time() - t0
t_minruntime[7] = min(t_elapsed)

np.savetxt('../data/ewstools_ricker.csv', ts.state[['residuals']].values)
np.savetxt('../data/ewstools_perfo.csv', t_minruntime, delimiter=',')