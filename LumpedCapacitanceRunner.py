import LumpedCapacitance as lc
import matplotlib.pyplot as plt

s = lc.LumpedCapacitanceModel(70, 128, 4, .02155, .07015, .07, .015, 790, 33)

t, temp = s.temperature_over_time(avg_convection_coefficient = 35.82, runtime = 3600)

print(s.required_convection_coefficient(60))
print(s.biot_number(.02155, 35, 1))

plt.plot(t, temp)
plt.show()