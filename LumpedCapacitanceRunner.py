#By Aadi Vasa
#aadivtx@gmail.com
import LumpedCapacitance as lc
import matplotlib.pyplot as plt

s = lc.LumpedCapacitanceModel(60, 26, 5, .0186, .0652, .047, .017, 790, 33, 60)

t, temp = s.temperature_over_time(avg_convection_coefficient = 99, runtime = 1800)

print(f"h: {s.required_convection_coefficient(60)}")
# print(s.biot_number(.0186, 101, 5.1))

area_inlet = .0652*.0208*5*2

print(f"v: {s.required_volumetric_flow_rate()/area_inlet}")
print(f"q: {s.required_volumetric_flow_rate()}")

plt.plot(t, temp)
plt.show()