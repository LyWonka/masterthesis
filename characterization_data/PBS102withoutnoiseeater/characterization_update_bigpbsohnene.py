import matplotlib.pyplot as plt
# from IPython.display import display, clear_output
# import thorlabs_apt as apt
# import pyvisa as visa
import pandas as pd
import numpy as np
import time
import os
# from ThorlabsPM100 import ThorlabsPM100
from scipy.optimize import curve_fit
import logging

#TODO blocking = true setzen

# from ctypes import *
# from datetime import datetime
# lib = cdll.LoadLibrary("C:\Program Files\IVI Foundation\VISA\Win64\Bin\TLPMX_64.dll")

# def do_measurement(start, stop, step_size, motor, power_meter_ref, power_meter_path1, power_meter_path2, file_name): #Macht eine Messung 
#     df = pd.DataFrame(columns=["position", "pref", "P1", "P2"])
#     motor.move_to(start, blocking=True) #auspassen, dass er nicht immer zu 0 fährt
#     angle = start #pos
#     power = c_longdouble()

#     while angle <= stop:
#         measures_power_meter_ref = []
#         measures_power_meter_path1 = []
#         measures_power_meter_path2 = []
#         for _ in range(10):
#             lib.TLPMX_measPower(power_meter_ref, byref(power), 1)
#             measures_power_meter_ref.append(power.value)
#             lib.TLPMX_measPower(power_meter_path1, byref(power), 1)
#             measures_power_meter_path1.append(power.value)
#             lib.TLPMX_measPower(power_meter_path2, byref(power), 1)
#             measures_power_meter_path2.append(power.value)
#         power_meter_ref_value = np.mean(measures_power_meter_ref)
#         power_meter_path1_value = np.mean(measures_power_meter_path1)
#         power_meter_path2_value = np.mean(measures_power_meter_path2)
#         df.loc[len(df)] = {"position": motor.position, "Pref": power_meter_ref_value, "P1": power_meter_path1_value, "P2": power_meter_path2_value}
#         motor.move_by(step_size, blocking=True)
#         angle += step_size #pos + step_size 
    
#     df.to_csv(file_name)
#     return(df)

def get_min_values(measure_table, col_of_interest, plot_filename, initial_guess=None):
    """

    RETURN
    ---
    initial_guess, min_fit_pos, min_fit_value, min_pos, min_W_value
    """
    # Schätzungen für die Anfangswerte
    # initial_guess = [1e-4, 0.1, 0, 1e-5] #rad
    if initial_guess is None:
        initial_guess =  np.array([np.max(np.abs(measure_table[col_of_interest])) - np.mean(measure_table[col_of_interest]), 1, 0, np.mean(measure_table[col_of_interest])]) #grad

    # Parametergrenzen festlegen
    # bounds = ([0, 0, -np.inf, 0], [np.inf, np.inf, np.inf, np.inf]) # rad
    # bounds = ([0, 0, 0, 0], [np.inf, 360, 360, np.inf]) #grad

    sigma = None#pd.Series([np.linalg.norm(initial_guess[-1]-m) for m in measure_table['measure']])
    #sigma = sigma.max() - sigma
    #sigma = sigma / sigma.sum()
    #sigma = sigma / sigma.max()
    #plt.plot(measure_table['position'], sigma, label='Gewichte', color='red')
    #plt.savefig(f"output/gewichte_{plot_filename.split('_')[1]}")
    #plt.close()
    popt, pcov = curve_fit(sin2_func_in_degrees, measure_table['position'], measure_table[col_of_interest], p0=initial_guess, sigma=sigma, maxfev = 10000)#, bounds=bounds)

    # Curve Fit anwenden mit erhöhter maxfev und bounds
    

    # Daten für den Fit berechnen
    fit_x = np.linspace(measure_table['position'].min(), measure_table['position'].max(), 1000)
    fit_y = sin2_func_in_degrees(fit_x, *popt)

    # Minimum des Fits finden
    min_index = np.argmin(fit_y)
    min_fit_pos = fit_x[min_index]
    min_fit_value = fit_y[min_index]


    # Gegebener minimaler Wert und Position aus Messung
    min_W_value = measure_table[col_of_interest].min()
    min_pos_row = measure_table[measure_table[col_of_interest] == min_W_value]
    min_pos = min_pos_row['position'].iloc[0]
    


    # Plot der Originaldaten und des Fits
    plt.figure(figsize=(10, 6))
    plt.scatter(measure_table['position'], measure_table[col_of_interest], label='Originaldaten', color='blue')
    plt.plot(fit_x, fit_y, label='Gewichteter Sinus^2 Fit', color='red')
    plt.axvline(x=min_fit_pos, color='green', linestyle='--', label='Minimum Position')
    #plt.scatter(min_pos, min_W_value, color='orange', zorder=5, label='Measured Minimum Value')
    plt.axvline(x=min_pos, color='orange', linestyle='--', label='Measured Minimum Value')
    # plt.text(min_pos+45, min_W_value, f'Messung Min: {min_W_value:.2e} \nbei {min_pos:.2f}°', color='black')
    # plt.text(min_fit_pos-45, min_fit_value, f'Fit Min: {min_fit_value:.2e} \nbei {min_fit_pos:.2f}°', color='black')
    #plt.text(0, 0, f'Fit Min: {min_fit_value:.2e} \nbei {min_fit_pos:.2f}°; Messung Min: {min_W_value:.2e} \nbei {min_pos:.2f}°', color='black')
    plt.title(f'Fit Min: {min_fit_value:.2e} \nbei {min_fit_pos:.2f}°; Messung Min: {min_W_value:.2e} \nbei {min_pos:.2f}°')
    plt.xlabel('Position')
    plt.ylabel('Measure')
    plt.savefig(plot_filename)
    plt.show()
    plt.close()

    return(initial_guess, min_fit_pos, min_fit_value, min_pos, min_W_value)

def inextremo(measure_table, col_of_interest):
    """
    Funktion gibt die beiden Minima und Maxima aus der 360 grad Messung  
    """
    measure_table_cp = measure_table.copy() #copy um originaldaten nicht zu überschreiben
    min_values = {}
    max_values = {}
    min_W_value = measure_table_cp[col_of_interest].min()
    min_pos_row = measure_table_cp[measure_table_cp[col_of_interest] == min_W_value]
    min_pos = min_pos_row['position'].iloc[0]
    min_values.update({min_pos:min_W_value})
    measure_table_cp.loc[measure_table_cp[abs(measure_table['position'] - min_pos) <= 10].index,col_of_interest]=np.nan
    min_W_value = measure_table_cp[col_of_interest].min()
    min_pos_row = measure_table_cp[measure_table_cp[col_of_interest] == min_W_value]
    min_pos = min_pos_row['position'].iloc[0]
    min_values.update({min_pos:min_W_value})

    max_W_value = measure_table_cp[col_of_interest].max()
    max_pos_row = measure_table_cp[measure_table_cp[col_of_interest] == max_W_value]
    max_pos = max_pos_row['position'].iloc[0]
    max_values.update({max_pos:max_W_value})
    measure_table_cp.loc[measure_table_cp[abs(measure_table['position'] - max_pos) <= 10].index,col_of_interest]=np.nan
    max_W_value = measure_table_cp[col_of_interest].max()
    max_pos_row = measure_table_cp[measure_table_cp[col_of_interest] == max_W_value]
    max_pos = max_pos_row['position'].iloc[0]
    max_values.update({max_pos:max_W_value})

    return(min_values, max_values)

def sin2_func_in_degrees(x, A, B, C, D):
    return A * np.sin(np.deg2rad(B * x + C))**2 + D


if __name__ == "__main__":
    timestamp = time.strftime('%Y%m%d%H%M')
    timestamp_oi = "202411081636"
    ordner = "PBS102withoutnoiseeater"
    logger = logging.getLogger(__name__)
    logging.basicConfig(filename=f'characterization__{timestamp}.log', level=logging.INFO)
    operator_name = input("Please insert the name of the operator: ")
    logger.info(f"Operator: {operator_name}")
    logger.info(f"Timestamp: {timestamp}")
    logger.info(f"Data from Timestamp: {timestamp_oi}")

    isolator_id = input("Please insert the id of the isolator: ")
    logger.info(f"Isolator: {isolator_id}")

    wavelength = int(input("Wavelength in nm:"))
    logger.info(f'Chosen wavelength: {wavelength}')

    # devices = apt.list_available_devices()
    # print(devices)

    # load motor 1
    # fragt Benutzer nach Seriennummer für Motor 1, dieser blinkt anschließend
    # Benutzer bestätigt, ob der richtige Motor geblinkt hat
    # wenn nicht, wird der Vorgang wiederholt

    # while(True):
    #     motor1_name = int(input("Serial Number Motor 1:"))
    #     motor1 = apt.Motor(motor1_name)
    #     motor1.identify()
    #     if input("Serial Number correct [yN]") == "y":
    #         break
    # # Set Home Parameter (Werte von https://github.com/qpit/thorlabs_apt/issues/24)
    # logger.info(f"Motor 1: {motor1_name}")
    # motor1.set_move_home_parameters(2,1,10,4)
    # motor1.move_home(blocking=True)

    # while(True):
    #     motor2_name = int(input("Serial Number Motor 2:"))
    #     motor2 = apt.Motor(motor2_name)
    #     motor2.identify()
    #     if input("Serial Number correct [yN]") == "y":
    #         break
    # logger.info(f"Motor 2: {motor2_name}")
    # motor2.set_move_home_parameters(2,1,10,4)
    # motor2.move_home(blocking=True)

    # list available power meter
    # find out if there are devices connected
    # deviceCount = c_ulong()
    # lib.TLPMX_findRsrc(0, byref(deviceCount))
    # print(str(deviceCount.value)+" powemeters connected")

    # load all power meters
    # if there are more devices connected, determine their names
    # if deviceCount.value >= 2:
    #     meterName1 = create_string_buffer(256)
    #     meterName2 = create_string_buffer(256)
    #     meterName3 = create_string_buffer(256)
    #     lib.TLPMX_getRsrcName(0, 0, meterName1)
    #     print(f"Powermeter ref: {meterName1.value}")
    #     lib.TLPMX_getRsrcName(0, 1, meterName2)
    #     print(f"Powermeter path 1: {meterName2.value}")
    #     lib.TLPMX_getRsrcName(0, 2, meterName3)
    #     print(f"Powermeter path 2: {meterName3.value}")
    # input("Please check if order of powermeters is correct")
    # logger.info(f"Powermeter ref: {meterName1.value}, owermeter path 1: {meterName2.value}, Powermeter path 2: {meterName3.value}")

    # Initialize the devices
    # power_meter_ref = c_ulong(0)
    # lib.TLPMX_init(meterName1, 1, 0, byref(power_meter_ref))
    # power_meter_path1 = c_ulong(0)
    # lib.TLPMX_init(meterName2, 1, 0, byref(power_meter_path1))
    # power_meter_path2 = c_ulong(0)
    # lib.TLPMX_init(meterName3, 1, 0, byref(power_meter_path2))

    # Set Wavelength (given in nm)
    # lib.TLPMX_setWavelength(power_meter_ref, c_double(wavelength), 1);
    # lib.TLPMX_setWavelength(power_meter_path1, c_double(wavelength), 1);
    # lib.TLPMX_setWavelength(power_meter_path2, c_double(wavelength), 1);

    # Set Unit- below sets to Watts
    # lib.TLPMX_setPowerUnit(power_meter_ref, 0, 1)
    # lib.TLPMX_setPowerUnit(power_meter_path1, 0, 1)
    # lib.TLPMX_setPowerUnit(power_meter_path2, 0, 1)

    # Justage des Analysators auf maximalen Durchlass
    # power = c_longdouble()
    input("Justage des Analysators auf maximalen Durchlass. Please start with path 1")
    #measure_table = do_measurement(start=0, stop=360, step_size=1, motor=motor2, power_meter_ref=power_meter_ref, power_meter_path1=power_meter_path1, power_meter_path2=power_meter_path2, file_name=f"{timestamp}_justage_analysator.csv")
    measure_table = pd.read_csv(os.path.join(ordner, f"{timestamp_oi}_justage_analysator.csv"))
    measure_table["P1/Pref"]=measure_table["P1"]/measure_table["Pref"]
    _, min_fit_pos_just_analysator, _, min_pos_just_analysator, _ = get_min_values(measure_table=measure_table, col_of_interest="P1/Pref", plot_filename=f"{timestamp}_justage_analysator.png")
    pos_without_isolator = min(min_fit_pos_just_analysator, (min_fit_pos_just_analysator+180) % 360)
    measure_table["diff"] = abs(measure_table["position"] - pos_without_isolator)
    power_min_without_isolator = measure_table[measure_table['diff'] == measure_table['diff'].min()]["P1"].values[0]
    logger.info(f"Justage des Analysators auf maximalen Durchlass")
    logger.info(f"min fitted position: {pos_without_isolator}, min measured position: {min_pos_just_analysator}")
    logger.info(f"Minimum Power without isolator: {power_min_without_isolator}")
    # motor2.move_to(min_fit_pos_just_analysator+90, blocking=True)

    # Justage für eine Wellenlänge ohne Isolator
    # measures_power_meter_ref = []
    # measures_power_meter_path1 = []
    # input("Please start with path 1")
    # for _ in range(10):
    #     lib.TLPMX_measPower(power_meter_ref, byref(power), 1)
    #     measures_power_meter_ref.append(power.value)
    #     lib.TLPMX_measPower(power_meter_path1, byref(power), 1)
    #     measures_power_meter_path1.append(power.value)
    # power_meter_ref_value = np.mean(measures_power_meter_ref)
    # power_meter_path1_value = np.mean(measures_power_meter_path1)

    # # TODO was tun wenn input unterschiedlich berechnet
    # input_power_reflect_path1 = power_meter_ref_value / reflect_value
    # input_power_trans_path1 = power_meter_path1_value / trans_value
    # print(f"powermeter ref value: {input_power_reflect_path1}, powermeter path1 value: {input_power_trans_path1}")
    # logger.info(f"Justage")
    # logger.info(f"powermeter ref value: {input_power_reflect_path1}, powermeter path1 value: {input_power_trans_path1}")

    # Faktorbestimmung für eine Wellenlänge ohne Isolator
    # measures_power_meter_ref = []
    # measures_power_meter_path2 = []
    input("Please change to path 2")
    # for _ in range(10):
    #     lib.TLPMX_measPower(power_meter_ref, byref(power), 1)
    #     measures_power_meter_ref.append(power.value)
    #     lib.TLPMX_measPower(power_meter_path2, byref(power), 1)
    #     measures_power_meter_path2.append(power.value)
    power_meter_ref_value = 0.0008382260565999999 # TODO np.mean(measures_power_meter_ref)
    power_meter_path2_value = 0.0039518245490000005 # TODO np.mean(measures_power_meter_path2)
    #faktor = power_meter_ref_value / power_meter_path2_value
    faktor = 4.714509311520728 # TODO power_meter_path2_value / power_meter_ref_value
    print(f"powermeter ref value: {power_meter_ref_value}, powermeter path2 value: {power_meter_path2_value}, Faktor: {faktor}")
    logger.info(f"powermeter ref value: {power_meter_ref_value}, powermeter path2 value: {power_meter_path2_value}, Faktor: {faktor}")

    # Design Orientation bestimmen
    input("please insert PBS in measuring holder")
    # measure_table = do_measurement(start=0, stop=360, step_size=1, motor=motor1,power_meter_ref=power_meter_ref, power_meter_path1=power_meter_path1, power_meter_path2=power_meter_path2, file_name=f"{timestamp}_design_orientation.png")
    measure_table = pd.read_csv(os.path.join(ordner,f"{timestamp_oi}_design_orientation.csv"))
    measure_table["P2/(Pref*faktor)"]=measure_table["P2"]/(measure_table["Pref"]*faktor)
    _, min_fit_do_pos, min_fit_value, min_pos, min_W_value = get_min_values(measure_table=measure_table,col_of_interest="P2/(Pref*faktor)", plot_filename=f"{timestamp}_design_orientation.png")
    # motor1.move_to(min_fit_do_pos+90, blocking=True)
    # measures_power_meter_path2 = []
    # for _ in range(10):
    #     lib.TLPMX_measPower(power_meter_path2, byref(power), 1)
    #     measures_power_meter_path2.append(power.value)
    power_do = 0.003336294834 # TODO np.mean(measures_power_meter_path2)
    logger.info("Design Orientation")
    logger.info(f"position: {min_fit_do_pos}, power: {power_do}")

    # Messung (Isolation)
    input("Please insert the isolator in path 2 now")
    # measure_table = do_measurement(start=0, stop=360, step_size=1, motor=motor1, power_meter_ref=power_meter_ref, power_meter_path1=power_meter_path1, power_meter_path2=power_meter_path2, file_name=f"{timestamp}_measurement.csv")
    measure_table = pd.read_csv(os.path.join(ordner,f"{timestamp_oi}_measurement.csv"))
    measure_table["P2/(Pref*faktor)"]=measure_table["P2"]/(measure_table["Pref"]*faktor)
    #plotet iso für jede position
    _, min_fit_iso_pos, min_fit_value, min_pos, min_W_value = get_min_values(measure_table=measure_table, col_of_interest="P2/(Pref*faktor)", plot_filename=f"{timestamp}_measurement.png")
   
    #motor1.move_to(min_fit_iso_pos+90)
    #measures_power_meter_ref = []
    #measures_power_meter_path2 = []
    #for _ in range(10):
    #    lib.TLPMX_measPower(power_meter_ref, byref(power), 1)
    #    measures_power_meter_ref.append(power.value)
    #    lib.TLPMX_measPower(power_meter_path2, byref(power), 1)
    #
    #    measures_power_meter_path2.append(power.value)
    power_meter_ref_value = 0.0008382260565999999 # TODO np.mean(measures_power_meter_ref)
    #power_meter_path2_value = np.mean(measures_power_meter_path2)
    iso = max(measure_table["P2/(Pref*faktor)"])
    #isolation = 10 * np.log(power_meter_path2_value / (power_meter_ref_value * faktor))
    isolation = 10 * np.log10(iso)
    logger.info("Isolation")
    logger.info(f"Isolation min: {isolation}, Power ref: {power_meter_ref_value}, Verhältnis Iso: {iso}")
    input(f"Does the value {isolation} for iso makes sense?")

    # Isolation at Design Orientation
    # motor1.move_to(min_fit_do_pos+90, blocking=True)
    # measures_power_meter_ref = []
    # measures_power_meter_path2 = []
    # for _ in range(10):
    #     lib.TLPMX_measPower(power_meter_ref, byref(power), 1)
    #     measures_power_meter_ref.append(power.value)
    #     lib.TLPMX_measPower(power_meter_path2, byref(power), 1)
    #     measures_power_meter_path2.append(power.value)
    power_meter_ref_value1 = 0.0007554080222 # TODO np.mean(measures_power_meter_ref)
    # power_meter_path2_value1 = np.mean(measures_power_meter_path2)
    iso_do1 = 0.9297273982736343 # TODO(power_meter_path2_value1 / (power_meter_ref_value1 * faktor))

    # motor1.move_to(min_fit_do_pos+180, blocking=True)
    # measures_power_meter_ref = []
    # measures_power_meter_path2 = []
    # for _ in range(10):
    #     lib.TLPMX_measPower(power_meter_ref, byref(power), 1)
    #     measures_power_meter_ref.append(power.value)
    #     lib.TLPMX_measPower(power_meter_path2, byref(power), 1)
    #     measures_power_meter_path2.append(power.value)
    power_meter_ref_value2 = 0.0007584648497999999 # TODO np.mean(measures_power_meter_ref)
    # power_meter_path2_value2 = np.mean(measures_power_meter_path2)
    iso_do2 = 1.1051613754032665 # TODO (power_meter_path2_value2 / (power_meter_ref_value2 * faktor))
    iso_do = max(iso_do1, iso_do2)
    #isolation_do = 10 * np.log(power_meter_path2_value / (power_meter_ref_value * faktor))
    isolation_do = 10 * np.log10(iso_do)
    logger.info(f"Iso design orientation: {iso_do}, iso1: {iso_do1}, iso2:{iso_do2} ,Isolation design orientation: {isolation_do} , Power ref do1:{power_meter_ref_value1},  Power ref do2:{power_meter_ref_value2}")


    # Messung (Insertion loss + Transmission)
    input("Transmission : Put Isolator in Transmission direction and press Enter")
    # motor1.move_to(min_fit_do_pos+90, blocking=True)
    # measures_power_meter_ref = []
    # measures_power_meter_path2 = []
    # for _ in range(10):
    #     lib.TLPMX_measPower(power_meter_ref, byref(power), 1)
    #     measures_power_meter_ref.append(power.value)
    #     lib.TLPMX_measPower(power_meter_path2, byref(power), 1)
    #     measures_power_meter_path2.append(power.value)
    power_meter_ref_value1 = 0.0007668456820999999 # TODO np.mean(measures_power_meter_ref)
    power_meter_path2_value1 = 0.003362683464 # TODO np.mean(measures_power_meter_path2)

    # motor1.move_to(min_fit_do_pos+180, blocking=True)
    # measures_power_meter_ref = []
    # measures_power_meter_path2 = []
    # for _ in range(10):
    #     lib.TLPMX_measPower(power_meter_ref, byref(power), 1)
    #     measures_power_meter_ref.append(power.value)
    #     lib.TLPMX_measPower(power_meter_path2, byref(power), 1)
    #     measures_power_meter_path2.append(power.value)
    power_meter_ref_value2 = 0.0007583540222999999 # TODO np.mean(measures_power_meter_ref)
    power_meter_path2_value2 = 7.673531361e-05 # TODO np.mean(measures_power_meter_path2)

    insertion_loss1 = 10 * np.log10((power_meter_ref_value1 * faktor) / power_meter_path2_value1)
    insertion_loss2 = 10 * np.log10((power_meter_ref_value2 * faktor) / power_meter_path2_value2)
    insertion_loss = max(insertion_loss1, insertion_loss2)
    print(f"insertion loss: {insertion_loss1}, insertion loss: {insertion_loss2}, pref1 = {power_meter_ref_value1},pref2 = {power_meter_ref_value2}, power2_1 = {power_meter_path2_value1}, power2_2 = {power_meter_path2_value2}")
    logging.info(f"insertion loss:{insertion_loss}")
    logging.info(f"insertion loss1: {insertion_loss1}, insertion loss2: {insertion_loss2}, pref1 = {power_meter_ref_value1},pref2 = {power_meter_ref_value2}, power2_1 = {power_meter_path2_value1}, power2_2 = {power_meter_path2_value2}")

    # trasmission an do
    transmission_design1 = power_meter_path2_value1 / (power_meter_ref_value1 * faktor)
    transmission_design2 = power_meter_path2_value2 / (power_meter_ref_value2 * faktor)
    if transmission_design1 > transmission_design2:
        power_meter_ref_value = power_meter_ref_value1
        transmission_design = transmission_design1
    else:
        power_meter_ref_value = power_meter_ref_value2
        transmission_design = transmission_design2
    print(f"Transmission at design orientation: {transmission_design}")
    logging.info(f"Transmission at design orientation: {transmission_design}, Power ref: {power_meter_ref_value}") #TODO print

    #maximale transmission
    # measure_table = do_measurement(start=0, stop=360, step_size=1, motor=motor1, power_meter_ref=power_meter_ref, power_meter_path1=power_meter_path1, power_meter_path2=power_meter_path2, file_name=f"{timestamp}_transmission.csv")
    measure_table = pd.read_csv(os.path.join(ordner,f"{timestamp_oi}_transmission.csv"))
    measure_table["P2/(Pref*faktor)"]=measure_table["P2"]/(measure_table["Pref"]*faktor)
    _, min_fit_trans_pos, _, min_trans_pos, _ = get_min_values(measure_table=measure_table,col_of_interest="P2/(Pref*faktor)", plot_filename=f"{timestamp}_transmission.png")
    min_trans_values, max_trans_values = inextremo(measure_table, "P2/(Pref*faktor)")
    print(f"min values: {min_trans_values}, max values: {max_trans_values}, error: {[minv+90 - maxv for minv,maxv in zip(min_trans_values.keys(), max_trans_values.keys())]}")
    logger.info("Transmission")
    logger.info(f"measured min values: {min_trans_values}, measured max values: {max_trans_values}, error: {[minv+90 - maxv for minv,maxv in zip(min_trans_values.keys(), max_trans_values.keys())]}")
    logger.info(f"fitted min values: {min_fit_trans_pos}")

    transmission_max = max(max_trans_values.values())
    delta_t = transmission_design / transmission_max
    print(f"Transmission Max: {transmission_max}, Delta T: {delta_t}")
    logger.info(f" Transmission Max: {transmission_max}, Delta T: {delta_t}, Power ref: {power_meter_ref_value}")

    # c) PER
    input(f"start PER, please change to path1 and press Enter")

    # motor1.move_to(min_fit_do_pos+90, blocking=True)
    # measures_power_meter_path1 = []
    # for _ in range(10):
    #     lib.TLPMX_measPower(power_meter_path1, byref(power), 1)
    #     measures_power_meter_path2.append(power.value)
    # power_meter_path1_value1 = np.mean(measures_power_meter_path1)

    # motor1.move_to(min_fit_do_pos+180, blocking=True)
    # measures_power_meter_path2 = []
    # for _ in range(10):
    #     lib.TLPMX_measPower(power_meter_path1, byref(power), 1)
    #     measures_power_meter_path1.append(power.value)
    # power_meter_path1_value2 = np.mean(measures_power_meter_path1)

    # if power_meter_path1_value1 > power_meter_path1_value2:
    #     power_meter_path1_value = power_meter_path1_value1
    #     design_oriantation = min_fit_do_pos+90
    # else:
    #     power_meter_path1_value = power_meter_path1_value2
    #     design_oriantation = min_fit_do_pos+180

    # motor1.move_to(design_oriantation, blocking = True)

    # 360 grad power max und power min * 2
    # measure_table = do_measurement(start=0, stop=360, step_size=1, motor=motor2, power_meter_ref=power_meter_ref, power_meter_path1=power_meter_path1, power_meter_path2=power_meter_path2, file_name=f"{timestamp}_per.csv")
    measure_table = pd.read_csv(os.path.join(ordner,f"{timestamp_oi}_per.csv"))
    measure_table["P1/Pref"]=measure_table["P1"]/measure_table["Pref"]
    min_values_per, max_values_per = inextremo(measure_table, "P1/Pref")
    _, min_fit_per_pos, _, _, _ = get_min_values(measure_table=measure_table, col_of_interest="P1/Pref", plot_filename=f"{timestamp}_per.png")
    print(f"measured min values: {min_values_per}, max values: {max_values_per}, error: {[minv+90 - maxv for minv,maxv in zip(min_values_per.keys(), max_values_per.keys())]}")
    print(f"fitted min values: {min_fit_per_pos}")
    logger.info("PER")
    logger.info(f"min values: {min_values_per}, max values: {max_values_per}, error: {[minv+90 - maxv for minv,maxv in zip(min_values_per.keys(), max_values_per.keys())]}")
    logger.info(f"fitted min values: {min_fit_per_pos}")

    # per berechnen TODO auf fit positionen aendern
    per = 10 * np.log10(max(max_values_per.values())/min(min_values_per.values()))
    
    print(f"PER: {per}")
    logger.info(f"PER: {per}")    # d) Pol rotation
    print(f"start Pol rotation measurement")

    delta_pos = abs(min(min_values_per) - pos_without_isolator)
    delta_pow = power_do / power_min_without_isolator
    print(f"delta_pos: {delta_pos}, delta_pow: {delta_pow}")
    logger.info(f"delta_pos: {delta_pos}, delta_pow: {delta_pow}")

    # close
    # lib.TLPMX_close(power_meter_ref)
    # lib.TLPMX_close(power_meter_path1)
    # lib.TLPMX_close(power_meter_path2)


    