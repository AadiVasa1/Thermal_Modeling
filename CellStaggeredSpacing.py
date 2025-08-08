import numpy as np
import math
from scipy.optimize import minimize_scalar
class CellStaggeredSpacing:

    def __init__(self, num_rows, cells_per_row, cell_diameter, cell_height, air_density=1.09, air_thermal_conductivity=.025, air_specific_heat_capacity = 1007, air_dynamic_viscosity = 19.35e-6, air_prandtl = .705, surface_prandtl = 0.703):
        self.N_l = num_rows #number of rows of cells/how many cells long from inlet to outlet
        self.N_t = cells_per_row #number of cells per row
        self.N = self.N_t * self.N_l #total number of cells in the segment

        self.D = cell_diameter #m, cell diameter
        self.l = cell_height #m, cell height

        self.rho = air_density #kg/m3, air density, optional
        self.k = air_thermal_conductivity #W/mK, air thermal conductivity, optional
        self.c_p = air_specific_heat_capacity #J/kgC, air specific heat capacity, optional
        self.mu = air_dynamic_viscosity #Pa s, air dynamic viscosity, optional

        self.Pr = air_prandtl #air prandtl number, defaulted at 33C value
        self.Pr_s = surface_prandtl #air prandtl number at cell surfaces, defaulted at 60C value

    """
    ###############################################################################################
    Helper Methods
    ###############################################################################################
    """
    @property
    def N_l(self):
        return self._N_l

    @N_l.setter
    def N_l(self, newVal):
        """
        Updates number of rows and the second correction factor, if necessary. Property decorated
        :param newVal: new number of rows
        :return:
        """
        if newVal >= 20:
            self.c2 = 1
        elif newVal >= 16:
            self.c2 = .99
        elif newVal >= 13:
            self.c2 = .98
        elif newVal >= 10:
            self.c2 = .97
        elif newVal >= 7:
            self.c2 = .95
        elif newVal >= 5:
            self.c2 = .92
        elif newVal >= 4:
            self.c2 = .89
        elif newVal >= 3:
            self.c2 = .84
        elif newVal >= 2:
            self.c2 = .76
        elif newVal >= 1:
            self.c2 = .64
        else:
            raise RuntimeError("Invalid N_l input")

        self._N_l = newVal

    def _s_d(self, s_t, s_l):
        """
        returns distance between cells in different rows
        :param s_t: transverse pitch (m)
        :param s_l: longitudinal pitch (m)
        :return: s_d (m)
        """
        return math.sqrt(s_t**2+(s_l/2)**2)

    def _v_max(self, D, s_t, s_l, inlet_v):
        """
        Calculates v_max given cell diameter, spacing, and inlet velocity

        :param D: cell diameter, m
        :param s_t: transverse pitch, m
        :param s_l: longitudinal pitch, m
        :param inlet_v: inlet velocity, m/s
        :return: max velocity, m/s
        """
        s_d = self._s_d(s_t, s_l)

        if s_d < D:
            raise RuntimeError("Invalid s_d or diameter, ensure geometry has space for air")
        elif 2 * (s_d - D) < (s_t - D):
            return s_t * inlet_v / (2 * (s_d - D)) #max velocity occurs on A2 (diagonal) plane
        else:
            return s_t * inlet_v / (s_t - D) #max velocity occurs on A1 (transverse) plane

    def _reynolds_max(self, v_max):
        """
        returns max reynolds number
        :param v_max: max velocity obtained via self.v_max()
        :return:
        """
        return self.rho * v_max * self.D / self.mu

    def _c1_m_table(self, s_t, s_l, re):
        """
        uses table 7.5 to find c1 and m values. only looks as the staggered portion
        :param s_t: transverse pitch, m
        :param s_l: longitudinal pitch, m
        :param re: reynolds number
        :return: [c1, m]
        """
        if re < 10 or re > 2e6:
            raise RuntimeError("Invalid reynolds number")

        if re < 1e2:
            return 0.90, 0.40
        elif re < 1e3:
            raise RuntimeError("Reynolds number between 1e2 and 1e3, approximate as a single isolated cylinder using SingleCellTest.py")
        elif s_t / s_l < 2 and re < 2e5:
            return .35 * ((s_t / s_l) ** (1/5)) , 0.60
        elif re < 2e5:
            return 0.40, 0.60
        elif re <= 2e6:
            return 0.022, 0.84

    def _p_correction_factor(self, P_t, P_l, re):
        """
        Uses figure 7.15 to get correction factor, might be a little off unless
        a better way to translate graph to code is found
        :param P_t: dimensionless transverse pitch, = s_t/D
        :param P_l: dimensionless longitudinal pitch, = s_l/D
        :param re: max reynolds number
        :return: correction factor
        """
        ptpl = P_t/P_l
        if ptpl <= 1.2:
            if re >= 1e5:
                return 1.685-0.571*ptpl
            elif re >= 1e4:
                return 1.343-.286*ptpl
            else:
                return 1.085-.071*ptpl
        else:
            if re >= 1e5:
                return 1.068-.056*ptpl
            elif re >= 1e4:
                return .867+.111*ptpl
            elif re >= 1e3:
                return .734+.222*ptpl
            else:
                return .445+.462*ptpl

    def _p_friction_factor(self, P_t, re):
        """
        Uses figure 7.15 to get friction factor, might be a little off unless
        a better way to translate graph to code is found
        :param P_t: dimensionless transverse pitch, = s_t/D
        :param re: max reynolds number
        :return: friction factor
        """
        if abs(P_t - 1.25) < abs(P_t - 1.5):
            if re <= 40:
                return -0.667 * re + 32.67
            elif re <= 60:
                return -0.1 * re + 10
            elif re <= 600:
                return -0.005556 * re + 4.33336
            elif re <= 2e3:
                return -0.000143 * re + 1.086
            elif re <= 6e3:
                return -0.00005 * re + 0.9
            elif re <= 4e4:
                return -0.000006 * re + 0.64
            elif re <= 2e5:
                return -0.000001 * re + 0.44
            else:
                return .2
        elif abs(P_t - 1.5) < abs(P_t - 2):
            if re <= 8e1:
                return -0.114286 * re + 11.14286
            elif re <= 2e2:
                return -0.008333 * re + 2.6666
            elif re <=1e3:
                return -0.0005 * re + 1.1
            elif re <= 2e4:
                return -0.000010526 * re + 0.61052
            elif re <= 1e5:
                return -0.0000025 * re + 0.45
            else:
                return .2
        elif abs(P_t - 2) < abs(P_t - 2.5):
            if re <= 2e2:
                return -0.002 * re + 1.2
            elif re <= 4e2:
                return -0.001 * re + 1
            elif re <= 6e3:
                return -0.0000357143 * re + 0.61428572
            elif re <= 8e4:
                return -0.0000027027 * re + 0.416216
            else:
                return .2
        else:
            if re <= 2e2:
                return -0.002 * re + 1
            elif re <= 2e3:
                return -0.0001111111 * re + 0.6222222
            elif re <= 6e4:
                return -0.0000034483 * re + 0.4068966
            else:
                return .2

    def _pressure_drop(self, s_t, s_l, v_max, reynolds):
        """
        Get pressure drop across bank of cells

        :param s_t: transverse pitch
        :param s_l: longitudinal pitch
        :param v_max: maximum velocity
        :param reynolds: max reynolds number
        :return: pressure drop in bars
        """
        P_t = s_t / self.D
        P_l = s_l / self.D
        chi = self._p_correction_factor(P_t, P_l, reynolds)
        f = self._p_friction_factor(P_t, reynolds)

        p_drop = self.N_l * chi * f * self.rho * (v_max**2) / 2 #N/m2
        return p_drop #Convert to bars

    def _re_guess_from_nu_req(self, nu_req, s_t, s_l):
        """
        given the required nusselt number, tries to figure out the closest integer reynolds number that'll achieve it
        uses c1_m_table method
        used with known_h_and_spacing method

        :param nu_req: required nusselt number
        :param s_t: transverse pitch
        :param s_l: longitudinal pitch
        :return: optimal reynolds number
        """

        def objective(re):
            # Reject unphysical Reynolds numbers
            if re < 100 or re > 2e6:
                return 1e9

            try:
                c1, m = self._c1_m_table(s_t, s_l, re)
            except RuntimeError:
                return 1e9

            nu = c1 * (re ** m) * (self.Pr ** .36) * ((self.Pr / self.Pr_s) ** .25) * self.c2
            return abs(nu - nu_req)

        result = minimize_scalar(objective, bounds=(100, 2e6), method='bounded')

        if not result.success:
            raise RuntimeError("Failed to converge on Reynolds number for given Nu")

        return result.x

    def _v_inlet_from_v_max(self, D, s_t, s_l, v_max):
        s_d = self._s_d(s_t, s_l)
        if 2 * (s_d - D) > (s_t - D):
            return v_max * (s_t - D) / s_t
        else:
            return v_max * 2 * (s_d - D) / s_t

    def heat_transfer_rate(self, inlet_velocity, h, s_t, cell_temp = 60, ambient_temp = 33):
        outlet_temp = (cell_temp-ambient_temp) * math.exp(-1 * math.pi * self.D * self.N * h / (self.rho * inlet_velocity * self.N_t * s_t * self.c_p))
        outlet_temp -= cell_temp
        outlet_temp *= -1
        # print(f"Outlet Temperature: {outlet_temp}")

        lmtd = (cell_temp - ambient_temp) - (cell_temp - outlet_temp)
        lmtd /= math.log((cell_temp - ambient_temp) / (cell_temp - outlet_temp))
        # print(lmtd)
        return self.N * h * math.pi * self.D * lmtd * self.l

    """
    ###############################################################################################
    Solver Methods
    ###############################################################################################
    """

    def known_v_and_spacing(self, inlet_velocity, s_t, s_l):
        """
        Returns the average convection coefficient or, if requested, the pressure drop
        through a bank of cells given a velocity and spacing
        :param inlet_velocity: velocity at the inlet, m/s
        :param s_t: transverse pitch, m
        :param s_l: longitudinal pitch, m
        :return: h and pressure drop
        """
        v_max = self._v_max(self.D, s_t, s_l, inlet_velocity)

        re_max = self._reynolds_max(v_max)

        c1, m = self._c1_m_table(s_t, s_l, re_max)

        nu_D = c1 * (re_max ** m) * (self.Pr ** .36) * ((self.Pr / self.Pr_s) ** .25) * self.c2

        h = nu_D * self.k / self.D

        pressure_drop = self._pressure_drop(s_t, s_l, v_max, re_max)

        return float(h), float(pressure_drop)

    def known_v_optimize_spacing(self, inlet_velocity, s_t_bounds, s_l_bounds, minimum_heat_transfer = 0, s_t_increment = .001, s_l_increment = .001, should_print = True):
        """
        given an inlet velocity, figures out an optimal spacing that maximizes h / pressure drop
        also can give a minimum heat transfer rate so that anything below gets filtered out
        give bounds to limit pitches

        :param inlet_velocity: velocity at the inlet, m/s
        :param s_t_bounds: an array characterizing transverse pitch bounds as [min, max], m
        :param s_l_bounds: an array characterizing longitudinal pitch bounds as [min, max], m
        :param minimum_heat_transfer: the minimum heat transfer rate that should be considered, default to 0
        :param s_t_increment: intervals at which to examine optimal s_t given bounds
        :param s_l_increment: intervals at which to examine optimal s_l given bounds
        :return:[optimal transverse spacing, optimal longitudinal spacing, maximum h_over_p value]
        """
        s_t_vals = [round(t,10) for t in np.arange(s_t_bounds[0], s_t_bounds[1], s_t_increment)]
        s_l_vals = [round(l,10) for l in np.arange(s_l_bounds[0], s_l_bounds[1], s_l_increment)]

        h_over_p = []
        for t in s_t_vals:
            row = [] #will be filled with h/p values with transverse pitch fixed and longitudinal pitch changing
            for l in s_l_vals:

                h, p = self.known_v_and_spacing(inlet_velocity, t, l)
                q = self.heat_transfer_rate(inlet_velocity, h, t)
                if q < minimum_heat_transfer:
                    h=0
                row.append(h/p)

            h_over_p.append(row)

        h_over_p = np.array(h_over_p)
        max_h_over_p = np.max(h_over_p)
        max_location_flattened = np.argmax(h_over_p)
        optimal_s_t_location, optimal_s_l_location = np.unravel_index(max_location_flattened, h_over_p.shape)

        optimal_s_t = s_t_vals[optimal_s_t_location]
        optimal_s_l = s_l_vals[optimal_s_l_location]

        if should_print:
            print(f"Transverse pitch: {optimal_s_t}")
            print(f"Longitudinal pitch: {optimal_s_l}")
            print(f"h/p val: {max_h_over_p}")

        return optimal_s_t, optimal_s_l, max_h_over_p

    def optimize_v_and_spacing(self, inlet_v_bounds, s_t_bounds, s_l_bounds, minimum_heat_transfer = 0, s_t_increment = .001, s_l_increment = .001, v_increment = .1, should_print = True):
        """
           given bounds for inlet velocity, figures out optimal spacing that maximizes h / (pressure drop * v * s_t)
           as this correlates to cooling effectiveness / fan power
           also can give a minimum h value so that anything below gets filtered out
           give bounds to limit pitches

           :param inlet_v_bounds: an array characterizing inlet velocity as [min, max], m/s
           :param s_t_bounds: an array characterizing transverse pitch bounds as [min, max], m
           :param s_l_bounds: an array characterizing longitudinal pitch bounds as [min, max], m
           :param minimum_heat_transfer: the minimum heat transfer rate that should be considered, default to 0
           :param s_t_increment: intervals at which to examine optimal s_t given bounds
           :param s_l_increment: intervals at which to examine optimal s_l given bounds
           :param v_increment: intervals at which to examine optimal inlet_v given bounds
           :return:[optimal transverse spacing, optimal longitudinal spacing, optimal inlet velocity, maximum h_over_p value]
           """
        v_vals = [round(t, 10) for t in np.arange(inlet_v_bounds[0], inlet_v_bounds[1], v_increment)]

        h_over_vp = []
        optimal_spacings = []
        for v in v_vals:
            try:
                s_t, s_l, h_over_p = self.known_v_optimize_spacing(v, s_t_bounds, s_l_bounds, minimum_heat_transfer, s_t_increment, s_l_increment, False)
            except RuntimeError:
                s_t = 1e9
                s_l = 1e9
                h_over_p = 0
            h_over_vp.append(h_over_p / v / s_t)
            optimal_spacings.append([s_t, s_l])

        h_over_vp = np.array(h_over_vp)
        max_h_over_vp = np.max(h_over_vp)
        max_location_index = np.argmax(h_over_vp)

        optimal_s_t = optimal_spacings[max_location_index][0]
        optimal_s_l = optimal_spacings[max_location_index][1]
        optimal_v = v_vals[max_location_index]

        if should_print:
            print(f"Transverse pitch: {optimal_s_t}")
            print(f"Longitudinal pitch: {optimal_s_l}")
            print(f"inlet_velocity: {optimal_v}")
            print(f"h/vp val: {max_h_over_vp}")

        return optimal_s_t, optimal_s_l, optimal_v, max_h_over_vp

    """
    ###############################################################################################
    Garbage Solver Methods
    I don't want to delete them </3
    ###############################################################################################
    """
    def known_h_and_spacing(self, req_h, s_t, s_l):
        """
        returns the required inlet velocity or pressure drop given the required h value and spacing
        :param req_h: required h value to reach
        :param s_t: transverse pitch, m
        :param s_l: longitudinal pitch, m
        :return: inlet velocity, m/s and pressure drop, bar
        """
        nu_req = req_h * self.D / self.k
        re_req = self._re_guess_from_nu_req(nu_req, s_t, s_l)

        v_max = re_req * self.mu / self.rho / self.D
        v_inlet = self._v_inlet_from_v_max(self.D, s_t, s_l, v_max)

        pressure_drop = self._pressure_drop(s_t, s_l, v_max, re_req)

        return float(v_inlet), float(pressure_drop)

    def known_h_optimize_spacing(self, req_h, s_t_bounds, s_l_bounds, s_t_increment = .001, s_l_increment = .001):
        """
        given a fixed h value and transverse and longitudinal bounds, determines the optimal spacing to hit that h value
        with a minimal flow rate and pressure drop (aka minimizing s_t * v * p)
        this is because this will minimize necessary fan power as

        Pumping Power is proportional to flow rate * pressure
        and flow rate is proportional to the inlet velocity * transverse pitch

        :param req_h: h value that must be hit
        :param s_t_bounds: an array characterizing transverse pitch bounds as [min, max], m
        :param s_l_bounds: an array characterizing longitudinal pitch bounds as [min, max], m
        :param s_t_increment: intervals at which to examine optimal s_t given bounds
        :param s_l_increment: intervals at which to examine optimal s_l given bounds
        :return: [optimal transverse spacing, optimal longitudinal spacing, minimum q_times_p value]
        """
        s_t_vals = [round(t, 10) for t in np.arange(s_t_bounds[0], s_t_bounds[1], s_t_increment)]
        s_l_vals = [round(l, 10) for l in np.arange(s_l_bounds[0], s_l_bounds[1], s_l_increment)]

        q_times_p = []
        for t in s_t_vals:
            row = []  # will be filled with h/p values with transverse pitch fixed and longitudinal pitch changing
            for l in s_l_vals:

                v, p = self.known_h_and_spacing(req_h, t, l)

                row.append(t*v*p)

            q_times_p.append(row)

        q_times_p = np.array(q_times_p)
        min_q_times_p = np.min(q_times_p)
        min_location_flattened = np.argmin(q_times_p)
        optimal_s_t_location, optimal_s_l_location = np.unravel_index(min_location_flattened, q_times_p.shape)

        optimal_s_t = s_t_vals[optimal_s_t_location]
        optimal_s_l = s_l_vals[optimal_s_l_location]

        print(f"Transverse pitch: {optimal_s_t}")
        print(f"Longitudinal pitch: {optimal_s_l}")
        print(f"q*p val: {min_q_times_p}")

        return optimal_s_t, optimal_s_l, min_q_times_p
