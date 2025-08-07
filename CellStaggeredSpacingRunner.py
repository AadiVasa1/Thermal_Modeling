import CellStaggeredSpacing as cs
import matplotlib.pyplot as plt
import numpy as np

c = cs.CellStaggeredSpacing(24, 18.6e-3)


ST=[] #Transverse Pitch
SL=[] #Longitudinal Pitch
x = 0.019
while x < 0.025:
    ST.append(x)
    x = x + 0.001
x = 0.0174
while x < .025:
    SL.append(x)
    x = x + 0.001

fig = plt.figure(figsize = (10,10))
ax = plt.axes(projection='3d')
X,Y = np.meshgrid(ST,SL)

h=[]
for j in SL:  # Outer loop for rows
    row = []  # Create an empty row
    for i in ST:  # Inner loop for columns
        # v, p = c.known_h_and_spacing(73, i, j)
        # row.append(v)  # Append an initial value (e.g., 0) to the row

        hh, p = c.known_v_and_spacing(2.189, i, j)
        row.append(hh/p)
    h.append(row)  # Append the filled row to the main array
#print(h)
Z=np.array(h)

surf = ax.plot_surface(X, Y, Z,  cmap = plt.cm.cividis)

# ax.scatter(.028,.022,22692, color = "red", s=50)

ax.set_xlabel('Transverse Pitch (m)', labelpad=20)
ax.set_ylabel('Longitudinal Pitch (m)', labelpad=20)
ax.set_zlabel('Convection Coefficient over Pressure Drop', labelpad=20)
ax.set_title('Cell Spacing vs Convection Coefficient')

plt.show()

# for i in range(15):
#     c.known_v_optimize_spacing(i+1,[.019, .025], [.0174, .025], 56)
#     print(f"Velocity: {i}")
#     print()

# for i in range(99):
#     print(i*5)
#     c.known_h_optimize_spacing(5*(i+1), [.022, .07], [.022, .07])
#     print()

t, l, _ = c.known_v_optimize_spacing(2.819, [.02, .025], [.0174, .025], 62.11, .0001, .0001)
print(c.known_v_and_spacing(2.819, t, l))
print(c._v_max(.0186, .02, .0174, 2.819))


# c.known_h_optimize_spacing(56,[.02, .03], [.0174, .03], .0001, .0001)
# v, p = c.known_h_and_spacing(56, .020, .0174)
# print(f"velocity: {v}")
# print(f"pressure drop: {p}")
# print(c.known_v_and_spacing(.2, .02, .0174))
# print(c._v_max(.0186, .02, .0174, v))
# print()
# print()

# c.known_v_optimize_spacing(.2, [.02, .0225], [.0174, .0225], 56)
