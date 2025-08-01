#By Aadi Vasa
#aadivtx@gmail.com
import math
from scipy.integrate import solve_ivp
class LumpedCapacitanceModel:

    def __init__(self, total_current, num_series, num_parallel, cell_diameter, cell_height, cell_weight, cell_resistance, cell_specific_heat,
                  ambient_temp = 33, air_temp_out = 33, air_density = 1.225, air_specific_heat = 1005):
        """
        Initializes with battery pack or cell properties, as well as expected ambient temperature

        :param total_current: total current running through the battery, A
        :param num_series: total number of cells in series
        :param num_parallel: total number of cells in parallel
        :param cell_diameter: cell diameter, m
        :param cell_height: cell height, m
        :param cell_weight: cell weight, kg
        :param cell_resistance: cell typical impedance, ohms
        :param cell_specific_heat: cell specific heat capacity, J/kgC
        :param ambient_temp: ambient temperature (usually 33C), C
        :param temp_out: temperature of air on outlet, C
        :param air_density: density of air, kg/m3
        :param air_specific_heat: specific heat of air, J/kgC
        """
        self.current = total_current

        self.num_series = num_series
        self.num_parallel = num_parallel

        num_cells = num_series * num_parallel
        cell_surface_area = cell_diameter * math.pi * cell_height

        self.lumped_surface_area = num_cells * cell_surface_area
        self.lumped_weight = num_cells * cell_weight

        self.lumped_resistance = self.equivalent_resistance(num_series, num_parallel, cell_resistance)

        self.cell_specific_heat = cell_specific_heat

        self.ambient_temp = ambient_temp
        self.air_temp_out = air_temp_out
        self.air_density = air_density
        self.air_specific_heat = air_specific_heat

    def equivalent_resistance(self, num_series, num_parallel, cell_resistance):
        """
        Returns the equivalent resistance of the battery pack

        :param num_series: number of cells in series
        :param num_parallel: number of cells in parallel
        :param cell_resistance: typical impedance per cell
        :return: equivalent resistance of the pack
        """
        return num_series * cell_resistance / num_parallel

    def temperature_over_time(self, avg_convection_coefficient, runtime = 30*60):
        """
        Given a convection coefficient and desired runtime, graphs average battery temperature over time

        :param avg_convection_coefficient: average convection coefficient of the air around the cells
        :param runtime: total length of time to analyze temperature, s
        :return: arrays with time steps and temperature
        """
        self.__convection_coefficient = avg_convection_coefficient

        sol = solve_ivp(self.__ivp_helper, (0,runtime), [self.ambient_temp])
        return sol.t, sol.y[0]

    def __ivp_helper(self, t, temp):
        """
        Used to help solve ivp for temperature over time function

        :param t: time, s
        :param temp: previous temperature
        :return: dT/dt
        """
        dTdt = (self.current ** 2) * self.lumped_resistance + self.__convection_coefficient * self.lumped_surface_area * ((self.ambient_temp+self.air_temp_out)/2 - temp)
        dTdt /= (self.lumped_weight * self.cell_specific_heat)
        return dTdt

    def required_convection_coefficient(self, maximum_cell_temperature = 60):
        """
        Given the maximum temperature the cells can reach, gives the minimum avg convection coefficient to cool the cells. Advised to test value in cfd and to add a safety factor to result.

        :param maximum_cell_temperature: max temperature you want the cells to reach, 60C for rules, 80C for cell limit
        :return: minimum average h value
        """
        return (self.current ** 2) * self.lumped_resistance / self.lumped_surface_area / (maximum_cell_temperature - (self.ambient_temp+self.air_temp_out)/2)

    def required_volumetric_flow_rate(self):
        """
        Knowing the worst desired temperature rise of air from initialization, calculates what the minimum volumetric flow rate is in m3/s
        """
        if self.air_temp_out == self.ambient_temp:
            raise RuntimeError("Outlet air temperature cannot equal ambient temperature")
        return (self.current ** 2) * self.lumped_resistance / self.air_specific_heat / (self.air_temp_out - self.ambient_temp) / self.air_density

    @staticmethod
    def biot_number(characteristic_length, convection_coefficient, thermal_conductivity=None,
                    thermal_conductivity_axial=None, thermal_conductivity_radial=None):
        """
        Calculates the biot number.
        Kind of useless since we don't have a way currently to determine radial thermal conductivity, but can use axial.
        If below 1, lumped capacitance model is justified.

        :param characteristic_length: the characteristic length of the object, usually the diameter for cylinders
        :param convection_coefficient: convection coefficient of air
        :param thermal_conductivity: just a general thermal conductivity, can be used for any object. If you have both axial and radial data for cell, use below variables instead
        :param thermal_conductivity_axial: For cells, axial thermal conductivity
        :param thermal_conductivity_radial: For cells, radial thermal conductivity.
        :return: biot number
        """

        if thermal_conductivity is None and (thermal_conductivity_axial is None or thermal_conductivity_radial is None):
            raise RuntimeError("Invalid thermal conductivity combination")

        if thermal_conductivity is None:
            thermal_conductivity = min(thermal_conductivity_axial, thermal_conductivity_radial)

        return convection_coefficient * characteristic_length / thermal_conductivity
