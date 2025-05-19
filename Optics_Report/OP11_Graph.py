# WARNING ========================

# SPAGHETTI CODE: TOO MUCH REPETITION

# IMPORTS ========================

import numpy as np
from scipy.odr import ODR, Model, RealData
import matplotlib.pyplot as plt

#[0.027,0.025,0.021,0.029,0.025,0.021] #separation errors withoud dividing by n

# VARIABLES ======================

y_data = np.array([3.71,3.46,3.24,2.99,2.74,2.50,2.29,2.02,1.77])#N
x_data = np.array([10,15,20,25,30,35,40,45,50]) #a
y_err = np.array([0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01])
x_err = np.array([0.000001,0.000001,0.000001,0.0000001,0.0000001,0.0000001,0.0000001,0.0000001,0.0000001])

# FUNCTIONS ======================

def linear_model(params, x):
    m, c = params
    return m * x + c

# MAIN CODE ======================

# Wrap the model in a scipy ODR Model
model = Model(linear_model)

# Combine the data and their uncertainties
data = RealData(x_data, y_data, sx=x_err, sy=y_err)

# Create an ODR object
odr = ODR(data, model, beta0=[0.0, 0.0002])  # Initial guess for m and c

output = odr.run()
output.pprint()  # Prints the fit parameters and statistics

# Extract parameters
m, c = output.beta
m_err, c_err = output.sd_beta

# Compute R²
y_pred = m * x_data + c
SS_res = np.sum((y_data - y_pred) ** 2)
SS_tot = np.sum((y_data - np.mean(y_data)) ** 2)
r2 = 1 - (SS_res / SS_tot)

print(f"Slope (m): {m} ± {m_err}")
print(f"Intercept (c): {c} ± {c_err}")

# Generate fitted line
x_fit = np.linspace(min(x_data), max(x_data), 100)
y_fit = linear_model([m, c], x_fit)

# Plot data with error bars
plt.errorbar(x_data, y_data, xerr=x_err, yerr=y_err, fmt='.', label="Data", capsize=3)

# Plot the fitted line
plt.plot(x_fit, y_fit,linestyle='-.', label=f"Fit: {m:.3f} ± {m_err:.4f} + {c:.2f} ± {c_err:.2f}")

x_text = x_data[len(x_data) // 2]-10
y_text = m * x_text + c + 0.15
plt.text(x_text, y_text, f"y = ({m:.3f} ± {m_err:.4f})x + {c:.2f} ± {c_err:.2f}\n$R^2 = {r2:.4f}$", 
         fontsize=11, color='black', bbox=dict(facecolor='white', alpha=0.8))

plt.legend(['Best fit', 'Data'])
plt.xlabel('n (mm^-1)')
plt.ylabel('d_measured (mm)')
plt.title('Distance from mirror to optical axis (d_measured) vs Number of lines per mm (n)')
plt.grid(color='k', linestyle=':', linewidth=0.25)
plt.show()
 
