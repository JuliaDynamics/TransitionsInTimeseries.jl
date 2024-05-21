import time
import numpy as np
import matplotlib.pyplot as plt
import ewstools
from ewstools.models import simulate_ricker

# Set seed for reproducibility
np.random.seed(0)

# Initialize time series and spectrum computation
series = simulate_ricker(tmax=1000, F=[0,2.7])
ts = ewstools.TimeSeries(data=series, transition=860)
ts.detrend(method='Lowess', span=0.2)
ts.state[['state','smoothing']].plot()
ts.compute_spectrum(rolling_window=0.5, ham_length=40)

# Initialize parameters for timing functions
rw = 0.5
n = 100
t_elapsed = np.zeros(10)

# Time functions (in a not very elegant way)
t0 = time.time()
for i in range(n):
    ts.compute_var(rolling_window=rw)
t_elapsed[0] = time.time() - t0

t0 = time.time()
for i in range(n):
    ts.compute_cv()
t_elapsed[1] = time.time() - t0

t0 = time.time()
for i in range(n):
    ts.compute_skew(rolling_window=rw)
t_elapsed[2] = time.time() - t0

t0 = time.time()
for i in range(n):
    ts.compute_kurt()
t_elapsed[3] = time.time() - t0

t0 = time.time()
for i in range(n):
    ts.compute_auto(rolling_window=rw, lag=1)
t_elapsed[4] = time.time() - t0

t0 = time.time()
for i in range(n):
    ts.compute_smax()
t_elapsed[5] = time.time() - t0

t0 = time.time()
for i in range(n):
    ts.compute_ktau()
t_elapsed[6] = time.time() - t0

t0 = time.time()
surro = ewstools.core.block_bootstrap(ts.state.residuals, n, bs_type='Stationary', block_size=10)
t_elapsed[7] = time.time() - t0

np.savetxt('ricker.csv', ts.state[['residuals']].values)
np.savetxt('t_elapsed.csv', t_elapsed)