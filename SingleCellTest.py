import math, numpy as np, matplotlib.pyplot as plt, pandas as pd
from scipy.optimize import minimize

"""
###############################################################################################
Cell Conditions
###############################################################################################
"""
m = .070 #kg, cell mass
rad = 21.55/2/1000 #m, cell radius
l = 70.15/1000 #m, cell height

r = .015 #ohm, cell internal resistance

area = math.pi*(rad**2) #m2, cell cross-sectional area
circumference = rad*2*math.pi
vol = area*l #m3, cell volume
density_cell = m/vol #kg/m3, cell avg density
    
"""
###############################################################################################
Testing Conditions
###############################################################################################
"""
def nu(re,pr):
    """
    Determines nusselt number given reynolds and prandtl numbers, unless velocity is 0, returns ambient value of 5
    :param re: reynolds number
    :param pr: prandtl number
    :return: nusselt number using hilpert cylindrical crossflow relationships
    """
    if v==0:
        return 5

    if pr<.7 or pr>50:
        raise RuntimeError("Invalid Prandlt Number")
    if re>400000 or re<.4:
        raise RuntimeError("Invalid Reynolds Number")

    if re>40000:
        return .027*(re**.805)*(pr**(1/3))
    if re>4000:
        return .193*(re**.618)*(pr**(1/3))
    if re>40:
        return .683*(re**.466)*(pr**(1/3))
    if re>4:
        return .911*(re**.385)*(pr**(1/3))
    if re>.4:
        return .989*(re**.330)*(pr**(1/3))

i = 70/3 #A, current through cell

runtime = 400 #s, length of test/sim

v = 0 #m/s, air velocity. If zero, h will be set to ambient value of 5
t_ambient = 33 #C, ambient temperature

density_air = 1.225 #kg/m3, air density at sea level
c_air = 1005 #J/kgC, specific heat capacity of air
k_air = .025 #W/mK, thermal conductivity of air
dyn_visc = 18.74e-6 #Pa*s, dynamic viscosity of air at 33C

l_char = rad*2 #m, characteristic length: cell diameter


reynolds = density_air*v*l_char/dyn_visc
prandtl = dyn_visc*c_air/k_air
nusselt = nu(reynolds, prandtl)

h = nusselt*k_air/l_char

"""
###############################################################################################
Testing Data
###############################################################################################
"""
time_measured = [] #Real test data, must be filled up by test_data_grabber
temp_measured = [] #Real test data at cell temp sensor location, must be filled up by test_data_grabber. If matching multiple curves, will default to top sensor values

temp_measured2 = [] #Real test data at middle of cell if matching with multiple sensors, must be filled up by test_data_grabber

temperature_sensor_location = -1 #-1 for top, 0 for middle, 1 for bottom, will be used if data given is only from one test sensor

def test_data_grabber_top_and_middle(filepath):
    """
    Grabs excel data and fills up three lists, used to develop best fitting curve using top and middle temp sensor data
    Assumes valid data in excel file
    :return: none, but creates time_measured and temp_measured_top/middle lists made with data points
    """
    global time_measured, temp_measured, t_ambient, temp_measured2
    df = pd.read_excel(filepath)

    time_measured = df.iloc[:, 0].to_numpy()
    temp_measured = df.iloc[:, 1].to_numpy()
    temp_measured2 = df.iloc[:,2].to_numpy()

    time_measured = time_measured - time_measured[0]
    t_ambient = temp_measured[0]

    if np.isnan(time_measured).any() or np.isnan(temp_measured).any() or np.isnan(temp_measured2).any():
        raise RuntimeError("NaN values detected while reading time and temp values from excel sheet, ensure column height of data is the maximum height of sheet")

def test_data_grabber_single(filepath):
    """
    Grabs excel data and fills up three lists, used to develop best fitting curve using only one sensor temp sensor data, may be less accurate
    Assumes valid data in excel file
    :return: none, but creates time_measured and temp_measured lists made with data points
    """
    global time_measured, temp_measured, t_ambient, temp_measured2
    df = pd.read_excel(filepath)

    time_measured = df.iloc[:, 0].to_numpy()
    temp_measured = df.iloc[:, 1].to_numpy()
    temp_measured2 = None

    time_measured = time_measured - time_measured[0]
    t_ambient = temp_measured[0]

    if np.isnan(time_measured).any() or np.isnan(temp_measured).any():
        raise RuntimeError("NaN values detected while reading time and temp values from excel sheet, ensure column height of data is the maximum height of sheet")

"""
###############################################################################################
Solver Conditions
###############################################################################################
"""
num_nodes = 2  #nodes across cell
dx = l / num_nodes

vol_sec = vol / num_nodes  #m3, volume of node-occupied space
lat_sa_sec = circumference * dx  #m2, lateral surface area of node-occupied space
m_sec = vol_sec * density_cell  #kg, mass of node-occupied space
p_sec = r * (i ** 2) / num_nodes  #heat generated in node-occupied space

def cell_simulation(k, c):
    """
    Algorithm
    Initial condition: cell starts at ambient temperature
    Boundary condition: bottom of cell held at ambient temperature
    :param k: Thermal conductivity of cell
    :param c: Specific heat capacity of cell
    :return: time and temp values at a particular location, can change based on what you want reported
    """
    if math.isnan(k):
        raise RuntimeError("k value reached NaN, check guess function, input or excel sheet")
    elif math.isnan(c):
        raise RuntimeError("c value reached NaN, check guess function, input, or excel sheet")

    dt = .5 * density_cell * c * (dx ** 2) / k  # s, length per timestep, uses fourier stability condition
    timesteps = math.ceil(runtime / dt)

    temp_cells = np.zeros((timesteps, num_nodes)) #array representing temperature across cell per timestep
    temp_cells[0] = t_ambient

    timevals = dt*np.arange(0, timesteps) #for plotting purposes, returns time at each timestep
    xvals = dx*np.arange(0, num_nodes) #for plotting purposes, 0=top of cell, last element = bottom of cell

    for t in range(temp_cells.shape[0]-1): #finding temperature at each timestep
        temp_cells[t+1,0] = (p_sec+h*(lat_sa_sec+area)*(t_ambient-temp_cells[t,0])+(k*area/dx)*(temp_cells[t,1]-temp_cells[t,0]))*dt/m_sec/c+temp_cells[t,0] #updating temperature at top of cell
        temp_cells[t+1,1:-1] = (p_sec+h*lat_sa_sec*(t_ambient-temp_cells[t,1:-1])+(k*area/dx)*(temp_cells[t,:-2]+temp_cells[t,2:]-2*temp_cells[t,1:-1]))*dt/m_sec/c+temp_cells[t,1:-1] #updating temperature for interior nodes of cell
        temp_cells[t+1,-1] = (p_sec+h*lat_sa_sec*(t_ambient-temp_cells[t,-1])+(k*area/dx)*(temp_cells[t,-2]+2*t_ambient-3*temp_cells[t,-1]))*dt/m_sec/c+temp_cells[t,-1] #updating temperature at bottom of cell

    return timevals, xvals, temp_cells

    # print(temp_cells)
    # _,axes = plt.subplots(1,2)
    #
    # axes[0].plot(timevals,temp_cells[:,num_nodes//2])
    # axes[1].plot(xvals, temp_cells[-1])
    # plt.show()

def thermal_prop_guess(k_guess, c_guess, bounds = [(1,20),(100,15000)]):
    """
    Finds the best fitting thermal conductivity and specific heat capacity of cells
    :param k_guess: guess for cell thermal conductivity
    :param c_guess: guess of cell specific heat
    :return: best fitting k and c values with the test data
    """

    guess = np.array([k_guess,c_guess]) #guesses for thermal conductivity and specific heat capacity of cell
    answer = minimize(loss, guess, bounds = bounds, method = "Powell") #returns best fitting combo of thermal conductivity and heat capacity of cell using data
    return answer.x[0], answer.x[1]

def loss(parameters):
    """
    helper for thermal_prop_guess
    :param parameters: initial guess for thermal conductivity and specific heat capacity of cell
    :return: least squares of curve (smaller is a better fitting pair of values)
    """
    k_guess, c_guess = parameters

    time_vals, x_vals, t_sim = cell_simulation(k_guess, c_guess) #run cell sim with guessed parameters

    if temp_measured2 is None: #will run if only matching one sensor
        location = temp_sensor_location() #uses actual temp sensor location to report temp at same point in sim and compare
        t_sim_interp = np.interp(time_measured, time_vals, t_sim[:,location]) #converts data points in sim to same domain as test
    else: #will run if we are matching both top and middle sensors
        t_sim_interp = np.interp(time_measured, time_vals, t_sim[:, 0])
        t2_sim_interp = np.interp(time_measured, time_vals, t_sim[:, num_nodes//2])

    if temp_measured2 is None: #returns least squares difference between sim and test, lower is better
        return np.mean((t_sim_interp-temp_measured)**2)
    else:
        return 0.5 * (np.mean((t_sim_interp - temp_measured)**2) + np.mean((t2_sim_interp - temp_measured2)**2)) #averages error if matching two curves
    
def temp_sensor_location():
    """
    returns location in temp_cells array to monitor for data matching based on what temp sensor is being looked at
    :param: temperature_sensor_location above, modify following instructions
    :return: location, index in temp_cells array
    """
    if temperature_sensor_location==-1:
        return 0
    elif temperature_sensor_location==0:
        return num_nodes//2
    elif temperature_sensor_location==1:
        return -1
    raise RuntimeError("Invalid value for temperature_sensor_location variable")

"""
###############################################################################################
Setter Functions for Above Properties, for situations where you're testing different situations at once
###############################################################################################
"""
def test_conditions_setter(current, time, air_velocity, ambient_temp, temp_sensor_location = -1, air_density = 1.225, air_specific_heat = 1005, air_thermal_conductivity = .025, air_dynamic_viscosity = 18.74e6):
    global i, runtime, v, t_ambient, density_air, c_air, k_air, dyn_visc, l_char, reynolds, prandtl, nusselt, h, temperature_sensor_location

    i = current
    runtime = time
    v = air_velocity
    t_ambient = ambient_temp

    density_air = air_density
    c_air = air_specific_heat
    k_air = air_thermal_conductivity
    dyn_visc = air_dynamic_viscosity

    l_char = rad*2

    reynolds = density_air*v*l_char/dyn_visc
    prandtl = dyn_visc*c_air/k_air
    nusselt = nu(reynolds, prandtl)

    h = nusselt*k_air/l_char

    temperature_sensor_location = temp_sensor_location

def cell_properties_setter(mass, diameter, height, internal_resistance):
    global m, rad, l, r, area, circumference, vol, density_cell
    
    m = mass
    rad = diameter/2
    l = height
    r = internal_resistance
    
    area = math.pi*(rad**2)
    circumference = rad*2*math.pi
    vol = area*l
    density_cell = m/vol

def solver_setter(number_nodes):
    global num_nodes, dx, vol_sec, lat_sa_sec, m_sec, p_sec
    num_nodes = number_nodes
    dx = l / number_nodes

    vol_sec = vol / number_nodes
    lat_sa_sec = circumference * dx
    m_sec = vol_sec * density_cell
    p_sec = r * (i ** 2) / number_nodes 

"""
###############################################################################################
Main Method, Request Outputs Here
###############################################################################################
"""

# test_data_grabber_top_and_middle("Cold_Plate_Test_Data_Filtered.xlsx")

# k_estimate, c_estimate = thermal_prop_guess(10, 1000)
# time_v, x_v, temp_c = cell_simulation(k_estimate,c_estimate)
# location = temp_sensor_location()
# print(f"Estimated Thermal Conductivity (W/mK): {k_estimate}")
# print(f"Estimated Specific Heat Capacity (J/kgC): {c_estimate}")
# print(f"Average Temperature of cell at end of run: {np.mean(temp_c[-1, :])}")

# plt.figure()
# plt.plot(x_v*1000, temp_c[-1,:])
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