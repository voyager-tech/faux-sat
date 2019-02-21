# Utilized Modules
import numpy as np
import sample.Requirements as Req
from sample.parts_list import Parts_List

from sample.Operational_Analyses import State_Determine
from sample.Operational_Analyses import Power_Consume, Total_Volume, Total_Mass

from sample.orbital_analyses import u

# Spacecraft State: lists operational conditions of subsystems/components


def bool_swap(state):
    # Change the state of a component into a boolean-ish value
    global boolish
    if state == 'on':
        boolish = 1.0
    if state == 'idle':
        boolish = 0.5
    if state == 'off':
        boolish = 0.0

    return boolish


def component_operations(rad, vel, jd, sun, GS, target, gs_inert):
    # TODO: Add check for array/matrix input size, return error if not valid
    sc_state = {}
    sc_state['solar_panel_mode'] = np.asmatrix(np.zeros((1, Req.Steps),
                                               np.dtype(object)))
    sc_state['radio_mode'] = np.asmatrix(np.zeros((1, Req.Steps),
                                         np.dtype(object)))
    sc_state['antenna_mode'] = np.asmatrix(np.zeros((1, Req.Steps),
                                           np.dtype(object)))
    sc_state['magnetorquer_mode'] = np.asmatrix(np.zeros((1, Req.Steps),
                                                np.dtype(object)))
    sc_state['payload_mode'] = np.asmatrix(np.zeros((1, Req.Steps),
                                           np.dtype(object)))
    sc_state['processor_mode'] = np.asmatrix(np.zeros((1, Req.Steps),
                                             np.dtype(object)))
    sc_state['eps_mode'] = np.asmatrix(np.zeros((1, Req.Steps),
                                       np.dtype(object)))
    sc_state['batt_capacity'] = np.asmatrix(np.zeros((1, Req.Steps),
                                            np.dtype(object)))
    sc_state['data_capacity'] = np.asmatrix(np.zeros((1, Req.Steps),
                                            np.dtype(object)))
    sc_state['time'] = jd

    # Initialize looping variables
    batt_capacity = Req.Batt_max_capacity
    batt_per = np.asmatrix(np.zeros((Req.Steps, 1), dtype=float))  # <-- get rid of eventually
    data_capacity = Parts_List['Payload']['storage'] + Parts_List['Processor']['storage']
    data_stored = 0.0 * u.byte
    data_per = np.asmatrix(np.zeros((Req.Steps, 1), dtype=float))  # <-- get rid of eventually

    for i in range(Req.Steps):
        # Determine the state of each component at each step
        # Based on orbital analysis data and component dependencies
        sc_state = State_Determine(sc_state, batt_capacity,
                                   Req.Batt_max_capacity, Req.Batt_min_percent,
                                   i, sun, GS, target, data_stored)

        # Power consumtion of the system:
        # Adds up all the power consumption and input and updates current value
        batt_capacity, batt_per[i, 0] = Power_Consume(sc_state, batt_capacity,
                                                      Req.Batt_max_capacity, i)
        sc_state['batt_capacity'][0, i] = batt_capacity

        # TODO: Link Budget
        data_uplink = 115.2 * u.kbyte / u.s  # bps - max rate
        data_downlink = 100.0 * u.kbyte / u.s  # bps - max rate

        # Set production of data from camera
        if sc_state['payload_mode'][0, i] != 'off':
            data_stored = data_stored + Parts_List['Payload']['data_picture']
        # Check if data_stored exceeds maximum storage
        if data_stored > data_capacity:
            data_stored = data_stored * 0.75  # delete 25% of data (CHANGE LATER)
        # Uplink / Downlink
        if sc_state['antenna_mode'][0, i] != 'off':
            # Downlink
            data_stored = (data_stored.to(u.byte) -
                           (data_downlink.to(u.byte / u.sec) * (60 * u.s)))
            # TODO: Uplink (currently incorrect)
#            data_stored = data_stored + (data_uplink * (60 * u.s))
        if data_stored.magnitude < 0.0:
            data_stored = (0.0 * u.byte)
        data_per[i, 0] = ((data_stored / data_capacity).to(u.byte)).magnitude
        sc_state['data_capacity'][0, i] = data_capacity

    # Component Volume Analysis
    sc_vol = Total_Volume()

    # Component Mass Analysis
    sc_mass = Total_Mass()

    # Set state data into "Boolean" for visalization
    boo = np.asmatrix(np.zeros((Req.Steps, 7), np.dtype(float)))
    for j in range(Req.Steps):
        boo[j, 0] = bool_swap(sc_state['eps_mode'][0, j])
        boo[j, 1] = bool_swap(sc_state['payload_mode'][0, j])
        boo[j, 2] = bool_swap(sc_state['radio_mode'][0, j])
        boo[j, 3] = bool_swap(sc_state['antenna_mode'][0, j])
        boo[j, 4] = bool_swap(sc_state['magnetorquer_mode'][0, j])
        boo[j, 5] = bool_swap(sc_state['processor_mode'][0, j])
        boo[j, 6] = bool_swap(sc_state['solar_panel_mode'][0, j])

    return sc_state, batt_per, data_per, boo, sc_vol, sc_mass










