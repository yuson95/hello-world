
import CoolProp.CoolProp as CP
import numpy as np
import pandas as pd
import math
import json


def compressor_energy_demand(p_min, compressor_stage, prate_cascade, e_isentropic):
    """Compressor Energy Demand, Assumptions: T_in_2+=40°C"""

    w_spec_cascade = 0
    p_out_last_stage = p_min*100000
    for i in range(1, compressor_stage+1):
        if i == 1:
            T_in_stage = 298.15  # K
        else:
            T_in_stage = 313.15  # K

        p_in_stage = p_out_last_stage
        p_out_stage = p_in_stage*prate_cascade
        entropy_in = CP.PropsSI('S', 'T', T_in_stage,
                                'P', p_in_stage, "hydrogen")
        h_in = CP.PropsSI('H', 'S', entropy_in, 'P', p_in_stage, "hydrogen")
        h_out_is = CP.PropsSI('H', 'S', entropy_in, 'P',
                              p_out_stage, "hydrogen")
        h_out = (h_out_is-h_in)/e_isentropic+h_in
        w_comp = h_out - h_in
        p_out_last_stage = p_out_stage

        w_spec_cascade += w_comp/3600000

    return w_spec_cascade


class Liquid_H2:

    """Design of liquid h2 station"""

    def __init__(self, senario):
        self.station_type = senario['station type']  # Station Type
        # Vehicle Service Pressure[bar]
        self.vehicle_service_pressure = senario['vehicle service pressure']
        # Refueling Station Size [kg/day]
        self.station_size = senario['station size']
        # Hydrogen Source
        self.hydrogen_source = senario['hydrogen source']
        # Dispensing Options to vehicle Tank
        self.dispensing_option = senario['dispensing option']
        # Compressor Stages
        self.compressor_stage = senario['compressor stage']
        # Vehicle Lingering time (min)
        self.t_linger = senario['vehicle Lingering time']
        # Utilization Rate
        self.utilization_rate = senario['utilization rate']

        # total vehicle fill time (min)
        self.t_fill = senario['total vehicle fill time']
        # Max. Dispensed Amount per Vehicle (kg)
        self.m_dispensed = senario['max. dispensed amount per vehicle']

        # Maximum number of vehicle back-to-back fills per HOSE during peak hours
        self.n_max = 60/(self.t_linger+self.t_fill)
        self.n_hoses = 1 if self.station_size == 212 else (
            2 if self.station_size == 420 else 4)  # Number of Hoses
        self.capacity_truck = 3800  # Capacity of Liquid Truck Delivery [kg]

    # Storage Vessels
    def get_p_vessels(self):
        """get the pressure of Storage Vessels"""

        type = ['max', 'min', 'switch']
        # lp=low pressure mp= middle pressure hp=high pressure
        level = ['lp', 'mp', 'hp']
        p_vessels = pd.DataFrame(np.arange(9).reshape(
            (3, 3)), index=type, columns=level)  # create a dataframe
        # Maximum Pressure service pressure= Cascade Dispensing=vehicle service pressure*14,5*1,25+1000;350*14,5*1,25+1000)
        p_vessels.loc['max', :] = self.vehicle_service_pressure*1.25+70
        # for low pressure vessel: Minimum Pressure=0,35*maximum pressure
        p_vessels.loc['min', 'lp'] = 0.35*p_vessels.loc['max', 'lp']
        # for middle pressure vessel: Minimum Pressure=0,65*maximum pressure
        p_vessels.loc['min', 'mp'] = 0.65*p_vessels.loc['max', 'mp']
        # for hp pressure vessel: Minimum Pressure=0,85*maximum pressure
        p_vessels.loc['min', 'hp'] = 0.85*p_vessels.loc['max', 'hp']
        p_vessels.loc['switch', 'lp'] = 0.906 * \
            p_vessels.loc['min', 'lp']  # * 0,871;0,895; 0,906
        p_vessels.loc['switch', 'mp'] = 0.895*p_vessels.loc['min', 'mp']
        p_vessels.loc['switch', 'hp'] = 0.871*p_vessels.loc['min', 'hp']

        return p_vessels

    # Dispenser Input
    def get_input_dispenser(self):
        """Dispenser Input"""

        # Vehicle Tank maximum pressure at end of cascade dispensing [bar]#At 25oC
        p_max_dipenser = 1.25*self.vehicle_service_pressure
        t_fill = self.t_fill  # Vehicle fill time through cascade (min)
        p_vessels = self.get_p_vessels()
        p_hp_switch = p_vessels.loc['switch', 'hp']
        v_tank_vehicle = self.n_hoses*self.m_dispensed / \
            CP.PropsSI("D", "P", p_hp_switch*100000, "T", 25 +
                       273.15, "hydrogen")  # Vehicle tank size (m³)
        input_dispenser = {'p_max_dipenser': p_max_dipenser,
                           't_fill': t_fill,
                           'v_tank_vehicle': v_tank_vehicle}

        return input_dispenser

    def refrigeration_design(self):
        # Refrigeration Unit Energy demand [kWh/kg] estimated according to HRSAM model
        """ Refrigeration Unit Energy demand [kWh/kg]"""
        if self.vehicle_service_pressure == 350:
            T_dipensing_max = -20
        else:
            T_dipensing_max = -40
        if self.dispensing_option == "cascade compressor":
            p_max_dipenser = self.get_input_dispenser()['p_max_dipenser']
            T_ref_design = 35
            COP = 0.9
            h_1 = CP.PropsSI('H', 'T', T_ref_design+273.15,
                             'P', p_max_dipenser*100000, "hydrogen")
            h_2 = CP.PropsSI('H', 'T', T_dipensing_max+273.15,
                             'P', p_max_dipenser*100000, "hydrogen")
            # Refrigeration Unit Energy demand [kWh/kg]
            w_refrigeration = (h_1 - h_2)/COP/3600000
            P_refrigeration = w_refrigeration * \
                self.m_dispensed*60/(self.t_linger+self.t_fill)
            n_ref = self.n_hoses if T_dipensing_max < 25 else 0
        else:
            n_ref = 0
            w_refrigeration = 0
            P_refrigeration = 0

        refrigeration_design = {'n_ref': n_ref,
                                'T_dipensing_max': T_dipensing_max,
                                'w_refrigeration': w_refrigeration,
                                'P_refrigeration': P_refrigeration}
        return refrigeration_design

    def get_input_compressor(self):
        """Compressor Input"""
        # SAE J2601   #input("Compressor Input: Please enter the Compressor Discharge Temperature Before Intercooling (C): ")
        T_comp_before_cooling = 150
        # SAE J2601    #input("Compressor Input: Please enter Maximum vehicle tank temperature (C): ")
        T_vehicle_tank_max = 85
        # Operating Storage Temperature (degrees C)   #input("Compressor Input: Please enter Operating Storage Temperature (degrees C): ")
        T_oper = 25
        T_comp_after_cooling = 40
        # Main Compressor Discharge Pressure (bar),Corresponds to cascade vessels maximum pressure
        p_discharge = self.get_p_vessels().loc['max', 'hp']+25
        # self.compressor_stages
        # Isentropic Compressor Efficiency for Refueling Station Compressor (%)
        e_isentropic = 0.8
        p_dewar = 6

        f_sizing = 1.1  # Compressor Motor Sizing Factor, H2A base case value is 110% based on typical industrial practices
        # Hydrogen Lost During Compression (% of Feed H2)
        e_lost_compressor = 0.005
        input_compressor = {'p_discharge': p_discharge,
                            'e_isentropic': e_isentropic,
                            'f_sizing': f_sizing,
                            'e_lost_compressor': e_lost_compressor}

        return input_compressor

    def compressor_design(self):
        """ Compressor Design"""
        input_compressor = self.get_input_compressor()

        # Cm, Minimum Compressor Capacity (kg/hr)
        m_compressor_min = self.n_hoses * \
            self.m_dispensed*60 / \
            (self.t_linger+self.t_fill) if self.dispensing_option == "cascade compressor" else 0
        p_min = input_compressor['p_dewar']
        # Design Compressor Flow Rate (kg/hr)
        m_compressor_design = m_compressor_min*1.13
        # Number of Compressors in Operation at Any Time
        n_compressor = math.ceil(m_compressor_design/35)

        prate_compressor = 10**((math.log10(input_compressor['p_discharge'])-math.log10(
            p_min))/self.compressor_stage) if self.dispensing_option == "cascade compressor" else 0  # Pressure Ratio per stage cascade compressor
        # specific energy requirement  for the compression (kWh/kg)
        w_spec_compressor = compressor_energy_demand(
            p_min, self.compressor_stage, prate_compressor, input_compressor['e_isentropic']) if self.dispensing_option == "cascade compressor" else 0
        # Theoretical Power Requirement for Each Compressor (kW)
        P_compressor_theo = w_spec_compressor*m_compressor_design
        # Actual Shaft Power Requirement for Each Compressor (kW)
        P_compressor_actual = P_compressor_theo if self.dispensing_option == "cascade compressor" else 0
        e_motor = (0.00008*(math.log(P_compressor_actual))**4-0.0015*(math.log(P_compressor_actual))
                   ** 3+0.0061*(math.log(P_compressor_actual))**2+0.0311*(math.log(P_compressor_actual))+0.7617) if self.dispensing_option == "cascade compressor" else 0
        if P_compressor_actual == 0:  # Motor Rating per Compressor (kW)
            P_compressor_motor = P_compressor_actual
        else:
            P_compressor_motor = P_compressor_actual / \
                e_motor*input_compressor['f_sizing']

        # specific total energy requirement  (kWh/kg)
        w_sepc_compressor_motor = w_spec_compressor / \
            e_motor if self.dispensing_option == "cascade compressor" else 0

        compresor_design = {'m_compressor_design': m_compressor_design,
                            'prate_compressor': prate_compressor,
                            'n_compressor': n_compressor,
                            'P_compressor_actual': P_compressor_actual,
                            'P_compressor_motor': P_compressor_motor,
                            'w_sepc_compressor_motor': w_sepc_compressor_motor}
        return compresor_design

    def liquid_cryogenic_tank_design(self):
        capacity_cryo_tank_desired = 2*self.capacity_truck if self.capacity_truck < self.station_size else(
            self.capacity_truck/4 if 4*self.station_size < self.capacity_truck/4 else(
                self.capacity_truck/3 if 3*self.station_size < self.capacity_truck/3 else(
                    self.capacity_truck/2 if 2*self.station_size < self.capacity_truck/2 else
                    self.capacity_truck
                )))
        capacity_cryo_tank_available = 2*self.station_size
        e_boiloff_storage = 0.005
        liquid_cryogenic_tank_design = {'capacity_cryo_tank_desired': capacity_cryo_tank_desired,
                                        'capacity_cryo_tank_available': capacity_cryo_tank_available,
                                        'e_boiloff_storage': e_boiloff_storage}

        return liquid_cryogenic_tank_design
    # design pump

    def pump_design(self):
        m_pump = self.n_hoses*self.m_dispensed*60 / \
            (self.t_linger+self.t_fill) if self.dispensing_option == "pump" else 0
        n_pump = math.ceil(
            m_pump/100) if self.dispensing_option == "pump" else 0
        # Cryo-Pump Boil-off (kg/day)
        m_cryo_pump_boil_off = 4 if self.dispensing_option == "pump" else 0
        # specific energy requirement  for the pump (kWh/kg), estimate according to MYRDD
        e_boiloff_pump = n_pump*m_cryo_pump_boil_off / \
            self.station_size/self.utilization_rate
        w_sepc_pump = 0.6 if self.dispensing_option == "pump" else 0
        # energy requirement  for the pump (kWh)
        W_pump = w_sepc_pump*self.station_size

        pump_design = {'m_pump': m_pump,
                       'n_pump': n_pump,
                       'm_cryo_pump_boil_off': m_cryo_pump_boil_off,
                       'e_boiloff_pump': e_boiloff_pump,
                       'w_sepc_pump': w_sepc_pump,
                       'W_pump': W_pump}

        return pump_design

    def cascade_calculation(self):
        """Design of Cascade System"""
        # t_oper = 24  #Hours the Refueling Station Is Operating(h)
        # m_hour_ava = self.m_hour_ava  #Average Hourly Demand During a Day (kg/hr)
        T_oper = 25  # ℃
        v_tank_vehicle = self.get_input_dispenser()['v_tank_vehicle']
        n_cascade_system = 1  # Optimum Number of Cascade Systems
        v_cascade_per_storage = 0.05  # Volume of Cascade storage [m³]
        variable = ['density_at_pmax', 'density_at_pmin', 'delta_density',  # list of needed variates
                    'm_vehicle_tank', 'm_cascade_vessel', 'v_min_cascade',
                    'n_tank', 'v_cascade', 'capacity_at_pmax',
                    'capacity_at_pmin', 'delta_capacity', 'percent_cascade_useable']
        # lp=low pressure mp= middle pressure hp=high pressure
        level = ['lp', 'mp', 'hp']
        p_vessels = self.get_p_vessels()
        cascade_design = pd.DataFrame(np.arange(36).reshape(
            (12, 3)), index=variable, columns=level)
        for i in range(3):
            # ### Low Pressure Cascade Storage Vessel
            # Density at Maximum Specified Unit Pressure (kg/m³)
            cascade_design.iloc[0, i] = CP.PropsSI(
                "D", "T", T_oper+273.15, "P", p_vessels.iloc[0, i]*100000, "hydrogen")
            # Density at Minimum Specified Unit Pressure (kg/m³)
            cascade_design.iloc[1, i] = CP.PropsSI(
                "D", "T", T_oper+273.15, "P", p_vessels.iloc[1, i]*100000, "hydrogen")
            # Difference Between Density at Maximum and Minimum Specified Unit Pressure - Cascade Tank (kg/m³)
            cascade_design.iloc[2, i] = cascade_design.iloc[0,
                                                            i] - cascade_design.iloc[1, i]
            # Maximum mass (kg) inside vehicle tank
            cascade_design.iloc[3, i] = v_tank_vehicle*CP.PropsSI(
                "D", "T", T_oper+273.15, "P", p_vessels.iloc[2, i]*100000, "hydrogen")
            # Mass required from  cascade vessel (kg)
            cascade_design.iloc[4, i] = cascade_design.iloc[3, i] if i == 0 else (
                cascade_design.iloc[3, i]-cascade_design.iloc[3, i-1])
            # Minimum Cascade Volume (m³)
            cascade_design.iloc[5, i] = cascade_design.iloc[4, i] / \
                cascade_design.iloc[2, i]
            # Number of Storage Tanks
            cascade_design.iloc[6, i] = math.ceil(
                cascade_design.iloc[5, i]/v_cascade_per_storage)
            # Cascade Volume [m³]
            cascade_design.iloc[7, i] = cascade_design.iloc[6,
                                                            i]*v_cascade_per_storage
            # Tank Capacity at Maximum Specified Unit Pressure (kg)
            cascade_design.iloc[8, i] = cascade_design.iloc[7, i]*CP.PropsSI(
                "D", "T", T_oper+273.15, "P", p_vessels.iloc[0, i]*100000, "hydrogen")
            # Tank Capacity at Minimum Specified Unit Pressure (kg)
            cascade_design.iloc[9, i] = cascade_design.iloc[7, i]*CP.PropsSI(
                "D", "T", T_oper+273.15, "P", p_vessels.iloc[1, i]*100000, "hydrogen")
            # Difference Between Capacity at Maximum and Minimum Specified Unit Pressure -  Cascade Tank (kg)
            cascade_design.iloc[10, i] = cascade_design.iloc[8,
                                                             i]-cascade_design.iloc[9, i]
            # Useable Capacity at Peak Demane
            cascade_design.iloc[11, i] = cascade_design.iloc[10,
                                                             i]/cascade_design.iloc[8, i]

        # Maximum Dispensible Amount from Cascade (kg)
        m_dispensiable = cascade_design.iloc[10, :].sum()
        # Number of Storage Tanks
        n_tank_total = cascade_design.iloc[6, :].sum()
        # Cascade Vessel Capacity/Unit (kg)
        capacity_at_pmax_sum = cascade_design.iloc[8, :].sum()/n_cascade_system
        # Cascade Size as a Percent of Average Daily Demand
        percent_cascad_of_demand = n_cascade_system*capacity_at_pmax_sum / \
            self.station_size
        cascade_system = {'cascade_design': cascade_design,
                          'm_dispensiable': m_dispensiable,
                          'n_tank_total': n_tank_total,
                          'capacity_at_pmax_sum': capacity_at_pmax_sum,
                          'percent_cascad_of_demand': percent_cascad_of_demand}
        return cascade_system


def calcu_investment(s):
    """to calculat all the investment"""
    liquid_sation = Liquid_H2(s)
    # Investment of refrigeration equipment
    f_install = 1.3
    rate_exchange = 0.82  # 2018 conversion of exchange rate
    # Investment of Refrigeration Equipment
    refrigeration_design = liquid_sation.refrigeration_design()
    n_ref = refrigeration_design['n_ref']
    P_refrigeration = refrigeration_design['P_refrigeration']
    T_dipensing_max = refrigeration_design['T_dipensing_max']
    cost_refrigeration = 1.25*11092 * \
        (100*n_ref*P_refrigeration/(T_dipensing_max+273.15)
         )**0.8579*f_install*rate_exchange
    # Investment of Liquid Hydrogen Cryogenic Storage Tank
    capacity_cryo_tank = liquid_sation.liquid_cryogenic_tank_design()[
        'capacity_cryo_tank_available']
    cost_cryo_tank = 991.89*(capacity_cryo_tank**0.692)*f_install*rate_exchange
    # Investment of Dispenser
    cost_dispenser = liquid_sation.n_hoses * 65000*f_install*rate_exchange
    # Investment of Cascade system
    cascade_system = liquid_sation.cascade_calculation()
    capacity_at_pmax_sum = cascade_system['capacity_at_pmax_sum']
    if s['vehicle service pressure'] == 350:
        cost_cascade = 1200*capacity_at_pmax_sum*f_install*rate_exchange
    else:
        cost_cascade = 1800*capacity_at_pmax_sum*f_install*rate_exchange
    # Investment of Evaporator
    cost_evaporator = liquid_sation.station_size / \
        24*1000+15000*f_install*rate_exchange

    # Investment of High Pressure Cryogenic LH2 Pump
    n_pump = liquid_sation.pump_design()['n_pump']
    cost_pump = n_pump*425000*f_install * \
        rate_exchange if liquid_sation.vehicle_service_pressure == 350 else n_pump * \
        700000*f_install*rate_exchange

    # Investment of gas compressor
    compressor_design = liquid_sation.compressor_design()
    n_compressor = compressor_design['n_compressor']
    P_compressor_motor = compressor_design['P_compressor_motor']
    if s['vehicle service pressure'] == 700:
        cost_compressor = n_compressor*40035 * \
            (P_compressor_motor**0.6038)*f_install*rate_exchange
    else:
        cost_compressor = n_compressor*40528 * \
            (P_compressor_motor**0.4603)*f_install*rate_exchange

    # cost of Overall Control and Safety Equipment
    cost_control = 100000*rate_exchange
    # Total Capital Investment (€)
    investment = {'cost_refrigeration': cost_refrigeration,
                  'cost_cryo_tank': cost_cryo_tank,
                  'cost_dispenser': cost_dispenser,
                  'cost_cascade': cost_cascade,
                  'cost_compressor': cost_compressor,
                  'cost_evaporator': cost_evaporator,
                  'cost_pump': cost_pump,
                  'cost_control': cost_control}
    investment_total = 0
    for v in investment.values():
        investment_total += 1.3*v
    investment['investment_total'] = investment_total

    return investment


def calcu_w_spec(s):
    liquid_sation = Liquid_H2(s)

    compresor_design = liquid_sation.compressor_design()
    w_sepc_compressor_motor = compresor_design['w_sepc_compressor_motor']

    pump_design = liquid_sation.pump_design()
    w_sepc_pump = pump_design['w_sepc_pump']

    refrigeration_design = liquid_sation.refrigeration_design()
    w_refrigeration = refrigeration_design['w_refrigeration']

    w_spec_total = w_sepc_compressor_motor + w_sepc_pump + w_refrigeration

    w_spec = {'w_sepc_compressor_motor':w_sepc_compressor_motor,
                'w_sepc_booster_motor':w_sepc_booster_motor,
                'w_refrigeration':w_refrigeration,
                'w_spec_total':w_spec_total}

    return w_spec



def calcu_cost_spec(s):
    investment = calcu_investment(s)
    liquid_sation = Liquid_H2(s)

    q = 0.08  # Zinssazt
    n = 10  # a, Betrachtungszeitraum
    a = ((1+q)**n)*q/((1+q)**n-1)  # Annuity factor
    annuity = investment*a
    demand_year = liquid_sation.station_size*liquid_sation.utilization_rate * \
        365  # Sum of hydrogen refueled per year [kg]
    CAPEX = annuity/demand_year

    # fix OPEX

    fix_OPEX = investment*0.1/demand_year
    # Variable OPEX
    w_spec_total = calcu_w_spec(s)['w_spec_total']
    price_electric = 0.1509  # €/KWh
    cost_electric = w_spec_total*price_electric

    e_lost_compressor = liquid_sation.get_input_compressor()[
        'e_lost_compressor']
    pump_design = liquid_sation.pump_design()
    e_boiloff_pump = pump_design['e_boiloff_pump']
    e_boiloff_tank = 2 * \
        liquid_sation.liquid_cryogenic_tank_design()['e_boiloff_storage']
    price_h2 = 9.5  # €/kg
    cost_lost = price_h2 * \
        (e_lost_compressor + e_boiloff_pump+e_boiloff_tank)
    var_OPEX = cost_electric + cost_lost
    # total cost
    cost_spec_total = CAPEX + fix_OPEX + var_OPEX
    cost_spec = {'CAPEX': CAPEX,
                 'fix_OPEX': fix_OPEX,
                 'var_OPEX': var_OPEX,
                 'cost_spec_total': cost_spec_total}

    return cost_spec
