import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
#
sys.path.append('../')
from PyUDV.read_MFprof import (read_MFprof, amplitude_from_MFprof_reading,
                               velocity_from_MFprof_reading)

path_data = 'src/data_sample.mfprof'

# #### Loading data
Data, Parameters, Info, Units = read_MFprof(path_data)
Amplitude_data = amplitude_from_MFprof_reading(Data, Parameters)
Velocity_data = velocity_from_MFprof_reading(Data, Parameters)
time = Data['profileTime']
z_coordinates = Data['DistanceAlongBeam']*1e2
indmax_time = np.argwhere(Data['transducer'] == 0)[1][0] - 1
indmax_z = -1
ind_bottom = 572

# #### plotting velocity
divnorm = colors.TwoSlopeNorm(vcenter=0, vmin=-0.2, vmax=0.2)
fig, ax = plt.subplots(1, 1, constrained_layout=True)
c = ax.pcolormesh(time[:indmax_time], z_coordinates[:indmax_z],
                  Velocity_data[:indmax_z, :indmax_time], cmap='PuOr', norm=divnorm,
                  rasterized=True, shading='auto')
ax.axhline(y=z_coordinates[ind_bottom], color='k', lw='0.5', ls='--')
ax.invert_yaxis()
fig.colorbar(c, label='Velocity [m/s]')
ax.set_xlabel('Time [s]')
ax.set_ylabel('DistanceAlongBeam [cm]')
# fig.draw_without_rendering()
plt.savefig('plots/Spatio_temporal_velocity.pdf', dpi=600)
plt.show()

# #### plotting amplitude
divnorm = colors.TwoSlopeNorm(vcenter=0, vmin=-0.1, vmax=0.1)
fig, ax = plt.subplots(1, 1, constrained_layout=True)
c = ax.pcolormesh(time[:indmax_time], z_coordinates[:indmax_z],
                  Amplitude_data[:indmax_z, :indmax_time],
                  rasterized=True, shading='auto', cmap='PuOr', norm=divnorm)
ax.axhline(y=z_coordinates[ind_bottom], color='k', lw='0.5', ls='--')
ax.invert_yaxis()
fig.colorbar(c, label='Amplitude [V]')
ax.set_xlabel('Time [s]')
ax.set_ylabel('DistanceAlongBeam [cm]')
plt.savefig('plots/Spatio_temporal_amplitude.pdf', dpi=600)
plt.show()
