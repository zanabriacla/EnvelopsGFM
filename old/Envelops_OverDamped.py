
import pandas as pd
import scipy.io
import matplotlib.pyplot as plt
import math
import numpy as np


def EventOnZgrid_up(Z1,Z2):
    #INITIALLY Z1 and Z2 are connected in parallel and then Z2 is disconnected
    # Code block
    Zgrid_initial = Z1 * Z2 / (Z1 + Z2)  # Initially Z1 and Z2 are connected in paralell
    Zgrid_final = Z1  # After the event ONLY Z1 is connected


    Delta_ZGrid = Zgrid_initial - Zgrid_final  # Variation in Zgrid
    print("EventOnZgrid_up")
    print("SCR initial", 1 / Zgrid_initial)
    print("SCR final", 1 / Zgrid_final)
    print("Delta X GRID ", Delta_ZGrid)
    return Delta_ZGrid,Zgrid_initial   # (optional)

def EventOnZgrid_down(Z1,Z2):
    # INITIALLY only Z1 is connected and then Z1 and Z2 are connected in parallel
    # Code block
    Zgrid_initial = Z1  # Z1 is connected
    Zgrid_final = Z1 * Z2 / (Z1 + Z2)  # Z1 and Z2 are connected
    Delta_ZGrid = (Zgrid_initial - Zgrid_final)


    print("EventOnZgrid_down")
    print("SCR initial", 1 / Zgrid_initial)
    print("SCR final", 1 / Zgrid_final)
    print("Delta X GRID ", Delta_ZGrid)
    return Delta_ZGrid,Zgrid_initial  # (optional)


def delay_signal(delay_ms,fs,signal):
    delay_samples = int((delay_ms / 1000) * 1 / fs)
    signal_delayed = np.concatenate((np.full(delay_samples, signal[0]), signal))[:len(signal)]
    return signal_delayed

def GetDeltaP(D,H,Xtotal_initial,P0):
    print("Second Order Function OverDAMPED")
    alpha = D / (2 * H)
    betha = wb / (2 * H * Xtotal_initial)
    p1 = (alpha - np.sqrt(alpha ** 2 - 4 * betha)) / 2
    p2 = (alpha + np.sqrt(alpha ** 2 - 4 * betha)) / 2

    epsilon = (p1 + p2) / (2 * np.sqrt(p1 * p2))  # Damping ratio of the second order response
    wn = np.sqrt(p1 * p2)
    wd = wn * math.sqrt(epsilon ** 2 - 1)

    print("epsilon", epsilon)
    print("Natural freq. (rad/s)", wn)
    print("Damped freq. (rad/s)", wd)

    A = (2 * H * (-p1) + D) / (p2 - p1) / (2 * H)
    B = (2 * H * (-p2) + D) / (p1 - p2) / (2 * H)

    Ppeak = P0 * DeltaXtotal / (Xtotal_initial * (np.sqrt((D ** 2) - 8 * H * wb / Xtotal_initial)))
    DeltaP = Ppeak * ((D - 2 * H * p1) * np.exp(-p1 * t_DeltaP) - (D - 2 * H * p2) * np.exp(-p2 * t_DeltaP))
    DeltaP1 = DeltaXtotal * P0 / Xtotal_initial * (
                A * np.exp(-p1 * t_DeltaP) + B * np.exp(-p2 * t_DeltaP))  # Same as Delta P only another way to get it
    P1 = P0 - DeltaP1
    P = P0 - DeltaP
    return DeltaP,Ppeak,epsilon

def GetTunnel(Ppeak):
    Tunnel = max(0.02, 0.05 * Ppeak)
    return Tunnel
#Variables needed to be fulfilled in order to implement the envelops

EventOnZgrid= "up" # Can be "up" or "down"
Z1=0.5 # Impedance Z1
Z2=0.5 # Impedance Z2

D=284.73 #Damping constant of the VSM control
H=10 #Inertia constant (s)
wb=314 # Base angular frequency(rad/s)
xtr=0.06 #Transformer reactance (pu)
Ugrid=1 # RMS voltage Ugrid (pu)
Uconv=1 # RMS voltage Uconverter (pu)
Xeff=xtr # effective reactance (pu)
EMT= True  # Can be "True" or "False" EMT is activated (20ms for the measures)
P0=0.5 # Initial power (pu)
Delta_ZGrid = 0 #Initial Value for Delta_ZGrid

if EventOnZgrid == "up":
    Delta_ZGrid = EventOnZgrid_up(Z1,Z2)[0]
elif EventOnZgrid == "down":
    Delta_ZGrid = EventOnZgrid_down(Z1, Z2)[0]
else:
    print("Invalid mode")


#second Order system

# Define the time vector for simulation
Start_Time = 0
End_Time = 2
NbPoints = 10000
t_DeltaP = np.linspace(Start_Time, End_Time, NbPoints)  # From Start_Time to End_Time
fs = (End_Time - Start_Time) / NbPoints  # Sampling frequency (Hz)

# Calculating VARIABLES that need to be defined to calculate DeltaP
Xtotal_initial = Xeff + EventOnZgrid_up(Z1, Z2)[1]  # X total initial that is equal to Xeff+Xgrid inital
DeltaXtotal = Delta_ZGrid  # Variation in Xtotal

#Defining margins for H and D

Ratio_H_D_UP = 0.1 # Use to have two more values for D and H: D*(1+Ratio_H_D_UP), H*(1+Ratio_H_D_UP)
Ratio_H_D_Down = 0.1 # Use to have two more values for D and H: D*(1-Ratio_H_D_Down), H*(1-Ratio_H_D_Down)

# Defining arrays to consider DeltaP for different H and D
DeltaP_array = []
Ppeak_array = []
Tunnel = []

D_array=[D,D*(1+Ratio_H_D_UP),D*(1-Ratio_H_D_Down)]
H_array=[H,H*(1-Ratio_H_D_UP),H*(1+Ratio_H_D_Down)]
print("Set of D values to be considered:",D_array)
print("Set of H values to be considered", H_array)

#Retrieving the second order response and the Tunnel that will be used in the Margins
DeltaP, Ppeak, epsilon = GetDeltaP(D_array[0], H_array[0], Xtotal_initial, P0)
Tunnel = GetTunnel(Ppeak)


#Creating Envelops
MargeUp=0.1 # This is the Margin up used in DeltaP*(1+-MarginUp)+Tunnel
MargeDown=0.1 # This is the Margin down used in DeltaP*(1+-MargeDown)+Tunnel


if DeltaP[0]>0:
    P_up_anal = DeltaP*(1+MargeUp)+Tunnel + P0
    P_down_anal = DeltaP*(1-MargeDown)-Tunnel+P0

else:
    P_up_anal = DeltaP*(1-MargeUp)+Tunnel+P0
    P_down_anal = DeltaP*(1+MargeDown)-Tunnel+P0

#Envelop of 50% of Delta before t=10ms and after that it takes DeltaP
P_50Prc = P0+ np.where(t_DeltaP < 0.01, DeltaP*0.5 , DeltaP)

#From different possibilities , we select the max value for Envelop UP and the MIN value for envelop DOWN
P_up_finale = np.maximum(P_up_anal,P_50Prc)
P_down_finale = np.minimum(P_down_anal,P_50Prc)
P_PCC= P0+DeltaP

#Adding delays when the simulation is done in EMT
if EMT:
    # Delay settings
    delay_ms = 20 # 20 ms of delay for EMT simulations
    shift_Time = 0
    # Pad with the initial value instead of zero
    initial_value = P_up_finale[0]
    P_up_finale = delay_signal(delay_ms, fs, P_up_finale)
    P_down_finale = delay_signal(delay_ms, fs, P_down_finale)
    P_PCC = delay_signal(delay_ms, fs, P_PCC)
    P_up_anal =delay_signal(delay_ms, P_up_anal)
    P_down_anal = delay_signal(delay_ms, P_down_anal)
else:
    shift_Time = 0
    delay_samples=0
    initial_value = P_up_finale[0]

# Create the plot
plt.figure(figsize=(8, 5))  # Set figure size
plt.plot(t_DeltaP+shift_Time,P_PCC, label="Ppcc", color='b', linestyle='-')  # First plot

plt.plot(t_DeltaP+shift_Time,P_up_anal, label="Pup_analytical", color='r', linestyle='--')  # First plot
plt.plot(t_DeltaP+shift_Time,P_down_anal, label="Pdown_analytical", color='m', linestyle='--')  # First plot
plt.plot(t_DeltaP+shift_Time,P_50Prc, label="P_50%", linewidth='3')  # First plot
plt.plot(t_DeltaP+shift_Time,P_down_finale, label="Pdown_final")  # First plot
plt.plot(t_DeltaP+shift_Time,P_up_finale, label="Pup_final")  # First plot


# Add vertical line at t = 10 seconds
plt.axvline(x=0.010, color='black', linestyle='--', label='t = 10ms')

# Add labels, title, and legend
plt.xlabel("sec")
plt.ylabel("P at PCC (pu)")
plt.title("SCR= "+str(1/Z1) + " Epsilon= " + str(round(epsilon,3)) + " ωd= " + " D= " + str(D) + " H= " +str(H) + " Xeff= " + str(Xeff))
plt.legend()  # Show legend

# Show the plot
plt.grid(True)  # Add grid for better visualization
plt.show()

