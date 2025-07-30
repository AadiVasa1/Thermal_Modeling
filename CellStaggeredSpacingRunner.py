import CellStaggeredSpacing as cs
import matplotlib.pyplot as plt
import numpy as np

c = cs.CellStaggeredSpacing(20, 21e-3)
# print(c.known_v_and_spacing(6,.0313,.0343))
# print(c.known_v_and_spacing(6,.0313,.0343, True))

ST=[] #Transverse Pitch
SL=[] #Longitudinal Pitch
x = 0.022
while x < 0.04:
    ST.append(x)
    SL.append(x)
    x=x+0.001

fig = plt.figure(figsize = (10,10))
ax = plt.axes(projection='3d')
X,Y = np.meshgrid(ST,SL)

h=[]
for j in SL:  # Outer loop for rows
    row = []  # Create an empty row
    for i in ST:  # Inner loop for columns
        row.append(c.known_v_and_spacing(2,i,j)/c.known_v_and_spacing(2,i,j, True))  # Append an initial value (e.g., 0) to the row
    h.append(row)  # Append the filled row to the main array
#print(h)
Z=np.array(h)

surf = ax.plot_surface(X, Y, Z,  cmap = plt.cm.cividis)

ax.scatter(.028,.022,22692, color = "red", s=50)

ax.set_xlabel('Transverse Pitch (m)', labelpad=20)
ax.set_ylabel('Longitudinal Pitch (m)', labelpad=20)
ax.set_zlabel('Convection Coefficient over Pressure Drop', labelpad=20)
ax.set_title('Cell Spacing vs Convection Coefficient')

plt.show()

c.known_v_optimize_spacing(2,[.022, .04], [.022, .04], minimum_h=100)
print(c.known_v_and_spacing(2,.028, .022))