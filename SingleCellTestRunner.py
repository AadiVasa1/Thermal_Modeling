import SingleCellTest as sct
import matplotlib.pyplot as plt, numpy as np

sct.solver_setter(200)
print(sct.num_nodes)

sct.test_data_grabber_top_and_middle("Cold_Plate_Test_Data_Filtered.xlsx")

k_estimate, c_estimate = sct.thermal_prop_guess(10, 1000)
time_v, x_v, temp_c = sct.cell_simulation(k_estimate,c_estimate)
location = sct.temp_sensor_location()
print(f"Estimated Thermal Conductivity (W/mK): {k_estimate}")
print(f"Estimated Specific Heat Capacity (J/kgC): {c_estimate}")
print(f"Average Temperature of cell at end of run: {np.mean(temp_c[-1, :])}")

plt.figure()
plt.plot(x_v*1000, temp_c[-1,:])
plt.xlabel("Distance from Cell Top (mm)")
plt.ylabel("Temperature (C)")
plt.show(block=False)

plt.figure()
plt.plot(sct.time_measured, sct.temp_measured, label = "Actual Data, Top")
plt.plot(sct.time_measured, sct.temp_measured2, label = "Actual Data, Middle")
plt.plot(time_v, temp_c[:,0], label = "Sim, Top")
plt.plot(time_v, temp_c[:, sct.num_nodes//2], label = "Sim, Middle")
plt.xlabel("Time (s)")
plt.ylabel("Temperature (C)")
plt.legend()
plt.show()


#modeling steady state
# runtime = 2500
# time_v, x_v, temp_c = cell_simulation(5.088,786.868)

# plt.figure()
# plt.plot(x_v, temp_c[-1,:])
# plt.xlabel("Distance from Cell Top (mm)")
# plt.ylabel("Temperature (C)")
# plt.show(block=False)

# plt.figure()
# plt.plot(time_measured, temp_measured, label = "Actual Data, Top")
# plt.plot(time_measured, temp_measured2, label = "Actual Data, Middle")
# plt.plot(time_v, temp_c[:,0], label = "Sim, Top")
# plt.plot(time_v, temp_c[:, num_nodes//2], label = "Sim, Middle")
# plt.xlabel("Time (s)")
# plt.ylabel("Temperature (C)")
# plt.legend()
# plt.show()




# test_data_grabber_single()

# k_estimate, c_estimate = thermal_prop_guess(15, 1000)
# time_v, x_v, temp_c = cell_simulation(k_estimate,c_estimate)
# location = temp_sensor_location()


# t,x,tt = cell_simulation(8, 769)

# plt.plot(time_measured, temp_measured, label = "Actual Data, Top")
# plt.plot(time_v, temp_c[:,location], label = "Sim, Top")
# plt.plot(t,tt[:,location])
# plt.legend()
# plt.show()
# print(f"Estimated Thermal Conductivity (W/mK): {k_estimate}")
# print(f"Estimated Specific Heat Capacity (J/kgC): {c_estimate}")

# plt.figure()
# plt.plot(x_v*1000, temp_c[-1,:])
# plt.plot(x*1000,tt[-1,:])
# plt.xlabel("Distance from Cell Top (mm)")
# plt.ylabel("Temperature (C)")
# plt.show()

