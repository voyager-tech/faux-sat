# Utilized Modules
from parts_list import Parts_List
import numpy as np
import Requirements as Req
from orbital_analyses import u


# Spacecraft State per step
def State_Determine(sc_state, batt_capacity, max_batt_cap,
                    min_batt, i, sun, GS, target, data_stored):
    """
    Determines the state of each spacecraft component based on previous state
    information and many conditional statements relating the components and
    outside data

    Parameters
    ----------
    sc_state : dict
        - All component states up to the current step
    batt_capacity : float
        - Previous step's stored power value
    max_batt_cap : float
        - Max battery capacity allowed by Requirements.py
    min_batt : float
        - Min battery percent allowed by Requirements.py
    i : int
        - Current step of propagation
    sun : int
        - 1 or 0 value for if the spacecraft is illuminated by the Sun
    GS : int
        - 1 or 0 value for if the spacecraft is in view of the Ground Station
    target : int
        - 1 or 0 value for if the spacecraft is within the target border
    data_stored : float
        - Previos step's stored data value

    Returns
    -------
    sc_state : dict
        - Component states for the current step, appended onto previous states
    """
    # Battery Condition
    if (batt_capacity.magnitude / max_batt_cap.magnitude) <= min_batt:
        batt_rule = 'off'
    elif (batt_capacity.magnitude / max_batt_cap.magnitude) > min_batt:
        batt_rule = 'on'

    # inSun Data
    if sun[0, i] == 1:  # Illuminated
        sc_state['solar_panel_mode'][0, i] = 'on'
    elif sun[0, i] == 0:  # Not illuminated
        sc_state['solar_panel_mode'][0, i] = 'off'

    # inContact Data
    # TODO: Add conditional for idle antenna with data flow
    if GS[0, i] == 1:  # Over ground station
        sc_state['magnetorquer_mode'][0, i] = 'on'
        sc_state['antenna_mode'][0, i] = 'idle'
        if data_stored >= Req.Data_min_store:
            sc_state['antenna_mode'][0, i] = 'on'
            sc_state['magnetorquer_mode'][0, i] = 'on'
    elif GS[0, i] == 0:  # Not over ground station
        sc_state['antenna_mode'][0, i] = 'off'

    # inCA Data
    # TODO: How to weight the times we can see target?
    # So we have bettery for pictures over pictures elsewhere?
    # 2 iterations? Make "picture" area larger than CA and weight that somehow?
    if GS[0, i] == 1:  # Over ground station
        # If over CA turn off radio and take pictures
        sc_state['payload_mode'][0, i] = 'off'
    elif target[0, i] == 1:  # Over target (CA)
        sc_state['payload_mode'][0, i] = 'on'
        sc_state['magnetorquer_mode'][0, i] = 'on'
    elif target[0, i] == 0:  # Not over target (CA)
        if batt_rule == 'on':  # More than 20% battery
            sc_state['payload_mode'][0, i] = 'idle'
            sc_state['magnetorquer_mode'][0, i] = 'on'
        elif batt_rule == 'off':  # Less than 20% battery
            sc_state['payload_mode'][0, i] = 'off'

    # Magnetorquer Conditional for 'off'
    if sc_state['payload_mode'][0, i] == 'off' and sc_state['antenna_mode'][0, i] == 'off':  # Not over GS or target (CA)
        sc_state['magnetorquer_mode'][0, i] = 'off'

    # Radio Dependency
    # TODO: When idle? if antenna is idle? never idle?
    # Currently idle if antenna is idle
    if sc_state['antenna_mode'][0, i] == 'on':  # Antenna is on
        sc_state['radio_mode'][0, i] = 'on'
    elif sc_state['antenna_mode'][0, i] == 'idle':  # Antenna is idle
        sc_state['radio_mode'][0, i] = 'idle'
    elif sc_state['antenna_mode'][0, i] == 'off':  # Antenna is off
        sc_state['radio_mode'][0, i] = 'off'

    # Processor Dependency
    if sc_state['payload_mode'][0, i] == 'on' or sc_state['payload_mode'][0, i] == 'idle' or sc_state['radio_mode'][0, i] == 'on' or sc_state['radio_mode'][0, i] == 'idle':  # If payload or radio are on / idle
        sc_state['processor_mode'][0, i] = 'on'
    elif sc_state['payload_mode'][0, i] == 'off' and sc_state['radio_mode'][0, i] == 'off':  # If payload and radio are off
        sc_state['processor_mode'][0, i] = 'idle'

    # EPS Dependency
    if sc_state['payload_mode'][0, i] == 'on' or sc_state['payload_mode'][0, i] == 'idle' or sc_state['solar_panel_mode'][0, i] == 'on' or sc_state['magnetorquer_mode'][0, i] == 'on' or sc_state['processor_mode'][0, i] == 'on' or sc_state['processor_mode'][0, i] == 'idle':  # If anything is on / idle
        sc_state['eps_mode'][0, i] = 'on'
    else:  # If everything is off
        sc_state['eps_mode'][0, i] = 'idle'

    # Battery Percentage: kill everything if below 20%, except solar panels
    if batt_rule == 'off':  # If battery charge is below 20%
        sc_state['eps_mode'][0, i] = 'off'
        sc_state['processor_mode'][0, i] = 'off'
        sc_state['magnetorquer_mode'][0, i] = 'off'
        sc_state['payload_mode'][0, i] = 'off'
        sc_state['antenna_mode'][0, i] = 'off'
        sc_state['radio_mode'][0, i] = 'off'
    return sc_state


# Power Consumption
def Power_Consume(state, pwr_step, max_cap, i):
    """
    Determines total power consumption at a single point in time, utilizing
    previous measurements to determine change in total stored power

    Parameters
    ----------
    state : dict
        - Component state info per step of propagation, from State_Determine
    pwr_step : float
        - Stored power from the previous step
    Max_cap : float
        - Maximum battery capacity, as defined in Requirements.py
    i : int
        - Current step

    Returns
    -------
    batt_capacity : float
        - Current stored power
    batt_percent : float
        - Current stored power as a percent of total capacity
    """
    pwr_cons = 0.0
    pwr_cons = pwr_cons + Parts_List['EPS']['power_consumption'][state['eps_mode'][0, i]]
    pwr_cons = pwr_cons + Parts_List['Payload']['power_consumption'][state['payload_mode'][0, i]]
    pwr_cons = pwr_cons + Parts_List['Radio']['power_consumption'][state['radio_mode'][0, i]]
    pwr_cons = pwr_cons + Parts_List['Antenna']['power_consumption'][state['antenna_mode'][0, i]]
    pwr_cons = pwr_cons + Parts_List['Magnetorquer']['power_consumption'][state['magnetorquer_mode'][0, i]]
    pwr_cons = pwr_cons + Parts_List['Processor']['power_consumption'][state['processor_mode'][0, i]]
    pwr_input = Parts_List['Solar_Panel']['power_output'][state['solar_panel_mode'][0, i]]

    # Battery Cacpacity Analysis
    # TODO: Check all units in these calculations - espically Wh
    batt_capacity = (pwr_step.magnitude +
                     ((pwr_input - pwr_cons).magnitude / 60)) * u.Wh
    if batt_capacity.magnitude > max_cap.magnitude:
        batt_capacity = max_cap
    if batt_capacity.magnitude <= 0.0:
        batt_capacity = 0.0 * u.Wh
    batt_percent = (batt_capacity / max_cap).magnitude
    return batt_capacity, batt_percent


# Volume Analysis
def Total_Volume():
    """
    Summs the volumes of all selected components and returns the total volume
    along with whether the max volume requirement has been met

    Returns
    -------
    sc_volume : list
        - First entry contains the total volume of the components
        - Second entry corresponds to if max volume requirement has been met

    See Also
    --------
    Total_Mass : Summation of component masses
    """
    # Component Volume Analysis
    sc_vol = []
    volume = 0
    for part in Parts_List:
        part_vol = (Parts_List[part]['volume'] * Parts_List[part]['quantity'])
        volume = volume + part_vol
    sc_vol.append(volume)

    # Determine if part volume fits into chassis volume
    if sc_vol[0].magnitude <= Req.Vol.magnitude:
        sc_vol.append('True')
    else:
        sc_vol.append('False')
    return sc_vol


# Mass Analysis
def Total_Mass():
    """
    Summs the masses of all selected components and returns the total mass
    along with whether the max mass requirement has been met

    Returns
    -------
    sc_mass : list
        - First entry contains the total mass of the components
        - Second entry corresponds to if max mass requirement has been met

    See Also
    --------
    Total_Volume : Summation of component volumes
    """
    # Component Mass Analysis
    sc_mass = []
    mass = 0
    for part in Parts_List:
        part_mass = (Parts_List[part]['mass'] * Parts_List[part]['quantity'])
        mass = mass + part_mass
    sc_mass.append(mass)

    # Determine if part mass exceeds the mass requirement
    if sc_mass[0].magnitude <= Req.Mass.magnitude:
        sc_mass.append('True')
    else:
        sc_mass.append('False')
    return sc_mass


def Link_Budget(sc_rad, gs_inert, i):
    # All taken from SMAD
    # Additional doccuments listed below
    # https://www.isispace.nl/product/full-ground-station-kit-for-vhfuhfs-band/  # Gound Antenna Description
    # https://www.isispace.nl/wp-content/uploads/2016/02/ISIS.GSK_.DS_.01.01_V2.2.1.pdf  # Ground Antenna Datasheet
    # https://www.isispace.nl/brochures/ISIS_TXS_Brochure_v.12.4.pdf  # ISIS S-Band Transmitter

    # A * indicates that the equations have been cheecked and are correct

    # %%############################# Uplink ##################################
    # Calculate values to be used by the rest of the function
    # Sender: Ground Station - ISIS S-Band Transmitter - in Madrid
    uG_Diam = 3.0                                           # m
    uG_Efficiency = 0.65                                    # TODO: Determine efficiency (decimal)
    uG_pwr = 4.0                                            # Watts - if mW need to change Eq. below
    uRb = 100.0e3                                           # bps - Date Rate

    # Reciever: Satellite
    uFreq_up = Parts_List['Antenna']['up_freq']             # GHz* - S-Band
    uS_Gain = Parts_List['Antenna']['gain']                 # dBi* - Max range
    uS_Noise = Parts_List['Antenna']['noise_figure']        # dBK <-- # TODO: Need to revamp with Eq. 16-24?

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Gateway Terminal - Type: Tracking
    uG_beamwidth = (21 / (uG_Diam * uFreq_up))              # deg*
    uG_Gain = (20.4 + (20 * np.log10(uFreq_up)) +
               (20 * np.log10(uG_Diam)) +
               (10 * np.log10(uG_Efficiency)))              # dBi*
    uG_Pwr_dB = (10 * np.log10(uG_pwr / 1))                 # dBW* - Input in W
#    uG_L_l = 6.0                                            # dB - Assumption*
    uG_EIRP = (uG_Pwr_dB + uG_Gain)  # - uG_L_l)                # dBW*
    # TODO: add EIRP per user

    # Propagation Range
    uDist = np.sqrt(((sc_rad[0, i] - gs_inert[0, i]) ** 2) +
                    ((sc_rad[1, i] - gs_inert[1, i]) ** 2) +
                    ((sc_rad[2, i] - gs_inert[2, i]) ** 2))  # km*
    uL_s = 92.45 + (20 * np.log10(uDist) +
                    (20 * np.log10(uFreq_up)))              # dB* - Space Loss
    uL_a = (10.0)                                           # dB - Atmospheric/Rain Losses - Assumption*
    # TODO: Find and pull from database for atmospheric losses
    uL_Net_Path = (uL_s + uL_a)                             # dB* - Net Path Loss

    # Satellite Antenna
    uS_L_l = 2.0                                            # dB - Line Loss on Satellite - Assumtion*
    uC_Power = (uG_EIRP + uS_Gain - uL_Net_Path - uS_L_l)   # dBW* - Recieved Carrier Power
    uG_T = (uS_Gain - uS_Noise)                             # dB/K*
    uL_comb = (uL_Net_Path + uS_L_l)                        # dB* - Combined Losses
    uC_No = (uG_EIRP + uG_T - uL_comb + 228.6)              # dB-Hz*
    uRb_l = (10 * np.log10(uRb))                            # dB-Hz*
    uEb_No = (uC_No - uRb_l)                                # dB* - Avaliable Uplink

    # %%############################ Downlink #################################
    # Calculate values to be used by the rest of the function
    # Sender: Satellite
    dS_Freq_down = Parts_List['Antenna']['down_freq']       # GHz* - S-Band
    dS_Gain = Parts_List['Antenna']['gain']                 # dBi* - Max of range
    dS_Pwr = Parts_List['Antenna']['power_consumption']['on']  # W

    # Reciever: Ground Station - ISIS S-Band Ground Station Kit - in Madrid
    dG_diam = 3.0                                           # m*
    dG_Efficiency = (0.2857)                                # decimal - Percent - Derived*
    dG_Noise = (22.75)                                      # dB/K - System Noise Temperature - Derived*
    dG_Gain = 31.35                                         # dBi* - S-Band
    dG_T = 8.6                                              # dB* - for 3m dish
    dRb = 115.2e3                                           # bps* - Max Data Rate, min rate = 14.4 kbps

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Satellite Antenna
    dL_bl = 5.5                                             # dB - Backoff & Line Loss - Assumption*
    dS_Pwr_dB = (10 * np.log10(dS_Pwr / 1))                 # dBW*
    dS_EIRP = (dS_Pwr_dB + dS_Gain - dL_bl)                 # dBW* - Per Channel
    # TODO: add EIRP per user

    # Propagation Range
    dDist = np.sqrt(((sc_rad[0, i] - gs_inert[0, i]) ** 2) +
                    ((sc_rad[1, i] - gs_inert[1, i]) ** 2) +
                    ((sc_rad[2, i] - gs_inert[2, i]) ** 2))  # km*
    dL_s = 92.45 + (20 * np.log10(dDist) +
                    (20 * np.log10(dS_Freq_down)))          # dB* - Space Losses
    dL_a = (7.0)                                            # dB - Atmospheric/Rain Losses - Assumption*
    # TODO: Find and pull from database for atmospheric losses
    dL_Net_Path = (dL_s + dL_a)                             # dB* - Net Path Losses

    # User Terminal
    dG_beamwidth = (21 / (dG_diam * dS_Freq_down))          # deg*
    dG_L_l = 2.0                                            # dB - Assumption*
    dC_Power = (dS_EIRP + dG_Gain - dL_Net_Path - dG_L_l)   # dBW* - Recieved Carrier Power
    dL_comb = (dL_Net_Path + dG_L_l)                        # dB*
    dC_No = (dS_EIRP + dG_T - dL_comb + 228.6)              # dB-Hz*
    dRb_l = (10 * np.log10(dRb))                            # dB-Hz*
    dEb_No = (dC_No - dRb_l)                                # dB*

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    udEb_No = -(10 * np.log10((10 ** (-(uEb_No) / 10)) +
                              (10 ** (-(dEb_No) / 10))))    # dB* - End-to-End Eb_No
    udL_mi = ()                                             # dB - Modem Implementation Loss - QPSK modulation
    ud_Eb_No_Req = ()
    Link_Margin = (udEb_No - udL_mi - ud_Eb_No_Req)         # dB

    # %%####################### Uplink / Downlink #############################

# =============================================================================
#
#   # Duplex Signal
#   # BPSK Downlink - Supported by GS
#   dData_Rate = (dRb * 10e-4)  # kbps
#   FEC = 0.75  # Forward Error Correction
#   ROF = 0.2  # Roll Off Factor
#   Mod_Ord = 1.0  # No of bits per symbol - BPSK
#   dSymbol_Rate = (dData_Rate / (FEC * Mod_Ord))
#   dBandwith = (((dSymbol_Rate + dSymbol_Rate * (ROF)) / 1000) * 2)  # MHz - *2 for Duplex
#
#   dBand = 90  # MHz
#   Chnl_Bandwith =()
#
# =============================================================================

    # Defined by Hardware
    Data_Down = 115.2e3 * u.byte / u.s                      # bps - max rate
    Data_Up = 100.0e3 * u.byte / u.s                        # bps - max rate
    return Data_Down, Data_Up
