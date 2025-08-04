import CellStaggeredSpacing as cs
import matplotlib.pyplot as plt
import numpy as np

c = cs.CellStaggeredSpacing(28, 18.6e-3)
# print(c.known_v_and_spacing(6,.0313,.0343))
# print(c.known_v_and_spacing(6,.0313,.0343, True))

# ST=[] #Transverse Pitch
# SL=[] #Longitudinal Pitch
# x = 0.022
# while x < 0.04:
#     ST.append(x)
#     SL.append(x)
#     x=x+0.001

# fig = plt.figure(figsize = (10,10))
# ax = plt.axes(projection='3d')
# X,Y = np.meshgrid(ST,SL)

# h=[]
# for j in SL:  # Outer loop for rows
#     row = []  # Create an empty row
#     for i in ST:  # Inner loop for columns
#         # v, p = c.known_h_and_spacing(50, i, j)
#         # row.append(1/(v*p))  # Append an initial value (e.g., 0) to the row

#         hh, p = c.known_v_and_spacing(2, i, j)
#         row.append(hh/i/p)
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

# for i in range(15):
#     c.known_v_optimize_spacing(i+1,[.022, .07], [.022, .07])
#     print(i)
#     print()

# for i in range(99):
#     print(i*5)
#     c.known_h_optimize_spacing(5*(i+1), [.022, .07], [.022, .07])
#     print()
# print(c.known_v_and_spacing(20,.028, .022))
# print(c.known_v_and_spacing(20,.028, .022))
# print()
# print(c.known_h_and_spacing(410, .028, .022))
# print(c.known_h_and_spacing(410, .028, .022)
c.known_h_optimize_spacing(56,[.02, .04], [.0174, .04])
v, p = c.known_h_and_spacing(56, .020, .0174)
print(f"velocity: {v}")
print(f"pressure drop: {p}")
print(c.known_v_and_spacing(.2, .02, .0174))
print()
print()
c.known_v_optimize_spacing(.2, [.02, .0225], [.0174, .0225], 56)
# print(c.known_h_and_spacing(50, .022, .020)) #(0.11450290661719202, 0.0005691083508591201)
# print(c.known_v_and_spacing(.11817839031716551, .022,.022))
