from SingleCellTest import SingleCellSimulation
import matplotlib.pyplot as plt, numpy as np

s = SingleCellSimulation()
s.cell_properties_setter(mass = .070, diameter = 21.55/1000, height = 70.15/1000, internal_resistance = 15/1000)
s.test_conditions_setter(current = 70/3, time = 400, temperature_sensor_location = 0)
s.test_data_grabber_top_and_middle("Cold_Plate_Test_Data_Filtered.xlsx")

s.solver_setter(50)
print(s.num_nodes)

k_estimate, c_estimate = s.thermal_prop_guess(10, 1000)
x_v, time_v, temp_c = s.cell_simulation(k_estimate,c_estimate)
location = s.temp_sensor_location
print(f"Estimated Thermal Conductivity (W/mK): {k_estimate}")
print(f"Estimated Specific Heat Capacity (J/kgC): {c_estimate}")
print(f"Average Temperature of cell at end of run: {np.mean(temp_c[-1, :])}")

plt.figure()
plt.plot(x_v*1000, temp_c[-1,:])
plt.xlabel("Distance from Cell Top (mm)")
plt.ylabel("Temperature (C)")
plt.show(block=False)

plt.figure()
plt.plot(s.time_measured, s.temp_measured, label = "Actual Data, Top")
plt.plot(s.time_measured, s.temp_measured2, label = "Actual Data, Middle")
plt.plot(time_v, temp_c[:,0], label = "Sim, Top")
plt.plot(time_v, temp_c[:, s.num_nodes//2], label = "Sim, Middle")
plt.xlabel("Time (s)")
plt.ylabel("Temperature (C)")
plt.legend()
plt.show()


#modeling steady state
# runtime = 2500
# x_v, time_v, temp_c = cell_simulation(5.088,786.868)

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




# s.test_data_grabber_single("Cold_Plate_Test_Data_Filtered.xlsx")

# s.solver_setter(50)
# print(s.num_nodes)

# location = s.temp_sensor_location

# k_estimate, c_estimate = s.thermal_prop_guess(15, 1000)
# x_v, time_v, temp_c = s.cell_simulation(k_estimate,c_estimate)

# x, t,tt = s.cell_simulation(8, 769)

# plt.plot(s.time_measured, s.temp_measured, label = "Actual Data, Top")
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

