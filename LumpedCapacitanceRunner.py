#By Aadi Vasa
#aadivtx@gmail.com
import LumpedCapacitance as lc
import matplotlib.pyplot as plt

s = lc.LumpedCapacitanceModel(70, 140, 4, .0186, .0652, .047, .017, 790, 33, 38)

t, temp = s.temperature_over_time(avg_convection_coefficient = 56, runtime = 3600)

print(s.required_convection_coefficient(60))
print(s.biot_number(.02155, 35, 1))
print(s.required_volumetric_flow_rate())

plt.plot(t, temp)
plt.show()