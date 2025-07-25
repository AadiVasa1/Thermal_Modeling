import math, numpy as np, matplotlib.pyplot as plt, pandas as pd
from scipy.optimize import minimize
from scipy.integrate import solve_ivp

class SingleCellSimulation:
    """
    
    This class can create a transient thermal simulation of a single cell across its length. 
    Can also add in air cooling and see how transient temperature gradient looks.

    If given test data in an excel sheet it can predict what the likely specific heat and axial thermal conductivity of the cell is.
    Only the first 3 columns matter, look at Cold_Plate_Test_Data_Filtered.xlsx for an example. Also ensure that those columns span the entire height of the sheet
    Can try to find results if given only one temperature sensor's data over time and location on cell height (less accurate) or if given top and middle temperature sensor vs time data (more accurate).

    Assumptions:
    1) Cell is uniformly heating across its length due to internal resistance
    2) Radial Heat flows fast enough for radial thermal gradients to not be of concern (Biot Number <<1)
    3) The floor on which the cell lies on alway remains at ambient temperature
    4) If there is no air cooling (air velocity = 0 m/s), h = 5

    """

    def __init__(self):
        """
        Does nothing lol, simply for object creation
        """
        pass

    """
    ###############################################################################################
    Setters, Call in this Order when Creating
    ###############################################################################################
    """
    def cell_properties_setter(self, mass, diameter, height, internal_resistance):
        """
        Physical Properties of Cell
        """        
        self.m = mass #kg, mass of cell
        self.rad = diameter/2 #m, radius of the cell
        self.l = height #m, height of the cell
        self.r = internal_resistance #ohms, internal resistance of the cell
        
        self.area = math.pi*(self.rad**2) #m2, cross sectional area of the cell
        self.circumference = self.rad*2*math.pi #m, circumference of the cell
        self.vol = self.area*self.l #m3, volume of the cell
        self.density_cell = self.m/self.vol #kg/m3, density of the cell
    
    def test_conditions_setter(self, current, time, air_velocity = 0, ambient_temp = 33, temperature_sensor_location = -1, air_density = 1.225, air_specific_heat = 1005, air_thermal_conductivity = .025, air_dynamic_viscosity = 18.74e6):
        """
        Single Cell Heating Test Conditions
        """
        self.i = current #A, current through cell
        self.runtime = time #s, length of test/sim
        self.v = air_velocity #m/s, air velocity. If zero, h will be set to ambient value of 5
        self.t_ambient = ambient_temp #C, ambient temperature, but can be replaced given data set

        self.density_air = air_density #kg/m3, air density at sea level
        self.c_air = air_specific_heat #J/kgC, specific heat capacity of air
        self.k_air = air_thermal_conductivity #W/mK, thermal conductivity of air
        self.dyn_visc = air_dynamic_viscosity #Pa*s, dynamic viscosity of air

        l_char = self.rad*2 #m, characteristic length, typically cell diameter
        reynolds = self.density_air*self.v*l_char/self.dyn_visc
        prandtl = self.dyn_visc*self.c_air/self.k_air
        nusselt = self.nu(reynolds, prandtl)

        self.h = nusselt*self.k_air/l_char #used for heat convection calcs

        self.temp_sensor_location = temperature_sensor_location #-1 for top, 0 for middle, 1 for bottom, will be used if data given is only from one test sensor

    def solver_setter(self, number_nodes = 50):
        """
        Sets up finite difference properties for cell_simulations solver
        """
        self.num_nodes = number_nodes #number of nodes across the cell
        self.dx = self.l / number_nodes

        self.vol_sec = self.vol / number_nodes #m3, volume of node-occupied space
        self.lat_sa_sec = self.circumference * self.dx #m2, lateral surface area of node-occupied space
        self.m_sec = self.vol_sec * self.density_cell #kg, mass of node-occupied space
        self.p_sec = self.r * (self.i ** 2) / number_nodes #heat generated in node-occupied space

    """
    ###############################################################################################
    Testing Data Setter
    ###############################################################################################
    """
    def test_data_grabber_top_and_middle(self, filepath):
        """
        Grabs excel data and fills up three lists, used to develop best fitting curve using top and middle temp sensor data
        Assumes valid data in excel file
        :return: none, but creates time_measured and temp_measured_top/middle lists made with data points
        """
        df = pd.read_excel(filepath)

        self.time_measured = df.iloc[:, 0].to_numpy() #real test data, times recorded
        self.temp_measured = df.iloc[:, 1].to_numpy() #real test data, temperatures recorded at the cell top over time
        self.temp_measured2 = df.iloc[:,2].to_numpy() #real test data, temperatures recorded at the cell middle over time

        self.time_measured = self.time_measured - self.time_measured[0]
        self.t_ambient = self.temp_measured[0]

        if np.isnan(self.time_measured).any() or np.isnan(self.temp_measured).any() or np.isnan(self.temp_measured2).any():
            raise RuntimeError("NaN values detected while reading time and temp values from excel sheet, ensure column height of data is the maximum height of sheet")

    def test_data_grabber_single(self, filepath):
        """
        Grabs excel data and fills up three lists, used to develop best fitting curve using only one sensor temp sensor data, may be less accurate
        Assumes valid data in excel file
        :return: none, but creates time_measured and temp_measured lists made with data points
        """
        df = pd.read_excel(filepath)

        self.time_measured = df.iloc[:, 0].to_numpy() #real test data, times recorded
        self.temp_measured = df.iloc[:, 1].to_numpy() #real test data, temperatures recorded at the given location over time
        self.temp_measured2 = None

        self.time_measured = self.time_measured - self.time_measured[0]
        self.t_ambient = self.temp_measured[0]

        if np.isnan(self.time_measured).any() or np.isnan(self.temp_measured).any():
            raise RuntimeError("NaN values detected while reading time and temp values from excel sheet, ensure column height of data is the maximum height of sheet")

    """
    ###############################################################################################
    Helper Methods
    ###############################################################################################
    """    
    def nu(self, re, pr):
            """
            Determines nusselt number given reynolds and prandtl numbers, unless velocity is 0, returns ambient value of 5
            :param re: reynolds number
            :param pr: prandtl number
            :return: nusselt number using hilpert cylindrical crossflow relationships
            """
            if self.v==0:
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

    @property
    def temp_sensor_location(self):
        """
        returns location in temp_cells array to monitor for data matching based on what temp sensor is being looked at
        :param: temperature_sensor_location above, modify following instructions
        :return: location, index in temp_cells array
        """
        if self._temp_sensor_location==-1:
            return 0
        elif self._temp_sensor_location==0:
            return self.num_nodes//2
        elif self._temp_sensor_location==1:
            return -1
        raise RuntimeError("Invalid value for temperature_sensor_location variable")

    @temp_sensor_location.setter
    def temp_sensor_location(self, new_location):
        if new_location not in (-1, 0, 1):
            raise RuntimeError("Invalid temp sensor location, should be -1, 0 or 1")
        else:
            self._temp_sensor_location = new_location
        
    """
    ###############################################################################################
    Solver
    ###############################################################################################
    """
    def cell_simulation(self, k, c):
        """
        Algorithm
        Initial condition: cell starts at ambient temperature
        Boundary condition: bottom of cell held at ambient temperature
        :param k: Thermal conductivity of cell
        :param c: Specific heat capacity of cell
        :return: xvals array and timevals array, along with calculated cell temperature 2d array with rows = timeval and columns = node (0 is top)
        """
        if math.isnan(k):
            raise RuntimeError("k value reached NaN, check guess function, input or excel sheet")
        elif math.isnan(c):
            raise RuntimeError("c value reached NaN, check guess function, input, or excel sheet")

        dx = self.dx
        runtime = self.runtime
        num_nodes = self.num_nodes
        t_ambient = self.t_ambient
        p_sec = self.p_sec
        lat_sa_sec = self.lat_sa_sec
        area = self.area
        m_sec = self.m_sec
        h = self.h

        xvals = dx*np.arange(0, num_nodes) #for plotting purposes, 0=top of cell, last element = bottom of cell

        def ivp_helper(t,temp):
            """
            returns dT/dt at each node per timestep to aid solve_ivp
            """
            dTdt = np.zeros_like(temp)
            dTdt[0] = (p_sec+h*(lat_sa_sec+area)*(t_ambient-temp[0])+(k*area/dx)*(temp[1]-temp[0]))/m_sec/c #updating temperature at top of cell
            dTdt[1:-1] = (p_sec+h*lat_sa_sec*(t_ambient-temp[1:-1])+(k*area/dx)*(temp[:-2]+temp[2:]-2*temp[1:-1]))/m_sec/c #updating temperature for interior nodes of cell
            dTdt[-1] = (p_sec+h*lat_sa_sec*(t_ambient-temp[-1])+(k*area/dx)*(temp[-2]+2*t_ambient-3*temp[-1]))/m_sec/c #updating temperature at bottom of cell
            return dTdt
        
        t_initial = np.full(num_nodes, t_ambient)
        sol = solve_ivp(ivp_helper, (0,runtime), t_initial) #builds transient solution where sol.t = timesteps and sol.y.T = temps(t, num_nodes)
        
        return xvals, sol.t, sol.y.T

    """
    ###############################################################################################
    Optimizer Methods for finding k and c given test data
    ###############################################################################################
    """
    def thermal_prop_guess(self, k_guess, c_guess, bounds = [(1,20),(100,15000)]):
        """
        Finds the best fitting thermal conductivity and specific heat capacity of cells
        :param k_guess: guess for cell thermal conductivity
        :param c_guess: guess of cell specific heat
        :return: best fitting k and c values with the test data
        """

        guess = np.array([k_guess,c_guess]) #guesses for thermal conductivity and specific heat capacity of cell
        answer = minimize(self.loss, guess, bounds = bounds, method = "Powell") #returns best fitting combo of thermal conductivity and heat capacity of cell using data
        return answer.x[0], answer.x[1]

    def loss(self, parameters):
        """
        helper for thermal_prop_guess
        :param parameters: initial guess for thermal conductivity and specific heat capacity of cell
        :return: least squares of curve (smaller is a better fitting pair of values)
        """
        k_guess, c_guess = parameters

        x_vals, time_vals, t_sim = self.cell_simulation(k_guess, c_guess) #run cell sim with guessed parameters

        if self.temp_measured2 is None: #will run if only matching one sensor
            location = self.temp_sensor_location #uses actual temp sensor location to report temp at same point in sim and compare
            t_sim_interp = np.interp(self.time_measured, time_vals, t_sim[:,location]) #converts data points in sim to same domain as test
        else: #will run if we are matching both top and middle sensors
            t_sim_interp = np.interp(self.time_measured, time_vals, t_sim[:, 0])
            t2_sim_interp = np.interp(self.time_measured, time_vals, t_sim[:, self.num_nodes//2])

        if self.temp_measured2 is None: #returns least squares difference between sim and test, lower is better
            return np.mean((t_sim_interp-self.temp_measured)**2)
        else:
            return 0.5 * (np.mean((t_sim_interp - self.temp_measured)**2) + np.mean((t2_sim_interp - self.temp_measured2)**2)) #averages error if matching two curves     
    