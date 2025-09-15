import CellStaggeredSpacing as cs
import matplotlib.pyplot as plt
import numpy as np

c = cs.CellStaggeredSpacing(26, 3, 2, 18.6e-3, 65.2e-3, uneven_stagger=True)

# ST=[] #Transverse Pitch
# SL=[] #Longitudinal Pitch
# x = 0.02
# while x < 0.025:
#     ST.append(x)
#     x = x + 0.0001
# x = 0.0174
# while x < .025:
#     SL.append(x)
#     x = x + 0.0001

# fig = plt.figure(figsize = (10,10))
# ax = plt.axes(projection='3d')
# X,Y = np.meshgrid(ST,SL)

# h=[]
# for j in SL:  # Outer loop for rows
#     row = []  # Create an empty row
#     for i in ST:  # Inner loop for columns
#         # v, p = c.known_h_and_spacing(73, i, j)
#         # row.append(v)  # Append an initial value (e.g., 0) to the row

#         hh, p = c.known_v_and_spacing(.9, i, j)
#         # q = c.heat_transfer_rate(2.189, h, i)

#         row.append(p)
#     h.append(row)  # Append the filled row to the main array
# #print(h)
# Z=np.array(h)

# surf = ax.plot_surface(X, Y, Z,  cmap = plt.cm.cividis)

# # ax.scatter(.028,.022,22692, color = "red", s=50)

# ax.set_xlabel('Transverse Pitch (m)', labelpad=20)
# ax.set_ylabel('Longitudinal Pitch (m)', labelpad=20)
# ax.set_zlabel('Convection Coefficient over Pressure Drop', labelpad=20)
# ax.set_title('Cell Spacing vs Convection Coefficient')

# plt.show()

q_req = .017 * ((62/5)**2)

# v, p = c.known_v_and_spacing(2.1, 31.6e-3, 11.9e-3)
# v_max = c._v_max(.0186, .0316, .0119, 2.1)
# print(p)
# print(v_max)

t, l, inlet_v, h, p, _ = c.optimize_v_and_spacing([.1, 10], [.0323, .046], [.0118, .0119], q_req, 60, 33, None, .0001, .0001, .1)
print(f"max velocity: {c._v_max(.0186, t, l, inlet_v)}")
print()

q, outlet_temp = c.heat_transfer_rate(inlet_v, h, t, outlet_temp_report = True)

print(f"total q dissipated: {q}")

last_cell_dissipation, q_cell_array, outlet_temp_array = c.last_cell_test(t, inlet_v, h, q_temp_graph_by_row = True)

print(f"last cell q dissipated: {last_cell_dissipation}")
print(f"q needed: {q_req}")
print(f"outlet temperature: {outlet_temp}")

print(c._s_d(.0316, .0119))

cell_row = []
for i in range(26):
    cell_row.append(i+1)

plt.plot(cell_row, q_cell_array)
plt.title('Q dissipated per Cell by Row')
plt.xlabel('Cell Row Number')
plt.ylabel('Q (W)')
plt.grid(True)
plt.show()

plt.plot(cell_row, outlet_temp_array)
plt.title('Air Temperature per Cell by Row')
plt.xlabel('Cell Row Number')
plt.ylabel('T (C)')
plt.grid(True)
plt.show()

# h, p = c.known_v_and_spacing(1.3, .0225, .01725)
# print()
# print(f"h: {h}")
# print(f"p: {p}")
# q = c.heat_transfer_rate(1.3, h, .0225)
# print(f"q dissipated: {q}")
# print(f"q needed: {q_req}")

# print()
# print(f"h: {h}")
# print(f"p: {p}")

# inlet_v = 2.75
# print(f"v: {inlet_v}")
# print()

# t, l, _ = c.known_v_optimize_spacing(inlet_v, [.02, .025], [.0174, .025], q_req, .0001, .0001)
# h, p = c.known_v_and_spacing(inlet_v, t, l)
# print()
# print(f"h: {h}")
# print(f"p: {p}")
# print(f"max velocity: {c._v_max(.0186, .02, .0174, inlet_v)}")
# print()
#
# q = c.heat_transfer_rate(inlet_v, h, t)
# print(f"q dissipated: {q}")
# print(f"q needed: {q_req}")



# c.known_h_optimize_spacing(134,[.02, .03], [.0174, .03], .001, .001)
# v, p = c.known_h_and_spacing(134, .020, .0174)
# print(f"velocity: {v}")
# print(f"pressure drop: {p}")
# print(c.known_v_and_spacing(.2, .02, .0174))
# print(c._v_max(.0186, .02, .0174, v))
# print()
# q = c.heat_transfer_rate(v, 134, .025)
# print(f"q dissipated: {q}")
# print(f"q needed: {.017 * ((70/4)**2) * 24 * 4}")

