# -*- coding: utf-8 -*-

# This file will contain the parts that will be added to the parts list
from sample.orbital_analyses import u

# Link to component data:
# https://docs.google.com/spreadsheets/d/1xwS2TIb4rqCcBVCQWCefzTvZD81DcXGUuBUhLh0Eu4M/edit#gid=64952105

# %% Payload (Camera)

Payload = {
        'name': 'NanoCam C1U 8mm',
        'manufacturer': 'GOMSpace',
        'quantity': 1,
        'modes': ['on', 'idle', 'off'],
        'mass': 167.0 * u.g,
        'volume': 434529.62 * (u.mm ** 3),
        'storage': 2.0 * u.Gbyte,
        'data_picture': 3.932160 * u.Mbyte,
        'power_consumption': {
                'idle': 1300.0 * u.mW,
                'on': 1300.0 * u.mW,
                'off': 0.0 * u.mW,
                },
        'voltage': 3.4 * u.V,
        'operational_alt': 650 * u.km
        }

# %% Solar Panel

Solar_Panel = {
        'name': '1U Cubesat Solar Panel 8V',
        'manufacturer': 'ISIS',
        'quantity': 4,
        'modes': ['on', 'off'],
        'power_output': {
                'on': 40.0 * u.W,  # <-- Go back and change back
                'off': 0.0 * u.W
                },
        'voltage_output': {
                'on': 8.0 * u.V,
                'off': 0.0 * u.V
                },
        'mass': 50.0 * u.g,
        'volume': 24010.0 * (u.mm ** 3),
        'cost': 2975.0,  # US Dollar
        'cell_efficiency': 0.3
        }

# %% Battery

Battery = {
        'name': 'Nanopower BP4 2P-2S',
        'manufacturer': 'GOMSpace',
        'quantity': 1,
        'mass': 258.0 * u.g,
        'volume': 181608.0 * (u.mm ** 3),
        'current_capacity': 5200.0 * u.mA * u.hour,
        'power_capacity': 38.5 * u.Wh,
        'voltage_nominal': 7.4 * u.V
        }

# %% EPS

EPS = {
       'name': 'Power Supply System EPSL',
       'manufacturer': 'NanoAvionics',
       'quantity': 1,
       'mass': 300 * u.g,
       'volume': 207399.25 * (u.mm ** 3),
       'cost': 3593.0,  # US Dollar
       'power_consumption': {
                'idle': 20.0 * u.W,  # <-- GO BACK AND CHECK THIS NUMBER
                'on': 25.0 * u.W,
                'off': 0.0 * u.W,
                },
       'input_efficiency': 0.9,
       'output_efficiency': 0.9
       }

# %% Magnetorquer

Magnetorquer = {
        'name': 'CubeADCS Magnetic',
        'manufacturer': 'CubeSpace',
        'quantity': 1,
        'mass': 225.0 * u.g,
        'volume': 293760.0 * (u.mm ** 3),
        'power_consumption': {
                'idle': 1.0 * u.mW,  # <-- CHECK IF THIS NUMBER IS NEEDED
                'on': 450.0 * u.mW,
                'off': 0.0 * u.mW,
                },
        }

# %% On-Board Computer (OBC)

Processor = {
        'name': 'NanoMind Z7000',
        'manufacturer': 'GOMSpace',
        'quantity': 1,
        'mass': 69.5 * u.g,
        'volume': 16900.0 * (u.mm ** 3),
        'storage': 3.2 * u.Gbyte,
        'clock_speed': 800.0,
        'power_consumption': {
                'idle': 2310.0 * u.mW,
                'on': 2520.0 * u.mW,
                'off': 0.0 * u.mW
                },
        }

# %% Chassis

Chassis = {
        'name': '1U NanoAvionics Structure',
        'manufacturer': 'NanoAvionics',
        'quantity': 1,
        'mass': 90.0 * u.g,
        'volume': 1000000 * (u.mm ** 3),
        'material': 'aluminum',
        }

# %% Radio

Radio = {
        'name': 'NanoCom AX 100U',
        'manufacturer': 'GOMSpace',
        'quantity': 1,
        'mass': 24.5 * u.g,
        'volume': 16900.0 * (u.mm ** 3),
        'power_consumption': {
                'idle': 3.96 * u.W,
                'on': 4.8 * u.W,
                'off': 0.0 * u.W
                },
        }

# %% Antenna

Antenna = {
        'name': 'NanoCom ANT 2000-DUP-2150',
        'manufacturer': 'GOMSpace',
        'quantity': 1,
        'mass': 110 * u.g,
        'volume': 193040 * (u.mm ** 3),  # <-- Go back and change
        'D_bandwidth': 80,  # MHz - Downlink
        'U_bandwidth': 80,  # MHz - Uplink
        'gain': 45,  # dBi - Decibel Isotropic
        'up_freq': 2.100,  # GHz
        'down_freq': 2.200,  # GHz
        'noise_figure': 2.2,  # dB
        'power_consumption': {
                'idle': 4.5 * u.W,  # <-- NEED TO GO BACK AND CHANGE
                'on': 10.0 * u.W,
                'off': 0.0 * u.W
                },
        }
