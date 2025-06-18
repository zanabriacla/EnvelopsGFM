
import pandas as pd
import scipy.io
import matplotlib.pyplot as plt
import math
import numpy as np
from sympy.abc import s, t
from sympy import symbols, Function, laplace_transform, inverse_laplace_transform, exp, simplify
from scipy.integrate import quad
from sympy import symbols, Function, Heaviside




def delay_signal(delay_ms,fs,signal):
    delay_samples = int((delay_ms / 1000) * 1 / fs)
    signal_delayed = np.concatenate((np.full(delay_samples, signal[0]), signal))[:len(signal)]
    return signal_delayed

def GetDeltaP_NotDELAYED(D_Damping,H,Xtotal_initial,P0,t_DeltaP):

    omega_n = np.sqrt(wb * Uconv * Ugrid / (2 * H * Xtotal))  # Natural frequency (rad/s)
    xi = omega_n * D_damping * Xtotal / (2 * wb * Uconv * Ugrid)

    if xi > 1:
        xi_over = xi  # Overdamped damping ratio (ξ > 1)
        #print("xi_over", xi_over)
    else:
        xi_under = xi  # Underdamped damping ratio (ξ < 1)
        #print("xi_under", xi_under)

    #print("xi", xi)

    #print("natural frq", omega_n)




    t= t_DeltaP

    # Common terms
    alpha = 2 * H * T_pll * RoCoF
    betha = (2 * H + D_damping * T_pll) / (2 * H * T_pll)

    D = (T_pll ** 2 * alpha * omega_n ** 2 * (T_pll * betha - 1)) / (
                1 - 2 * xi * omega_n * T_pll + omega_n ** 2 * T_pll ** 2) * -1
    A = alpha * betha
    B = alpha * (2 * T_pll * betha * omega_n * xi - T_pll * omega_n ** 2 - betha) / (
                1 - 2 * xi * omega_n * T_pll + omega_n ** 2 * T_pll ** 2)
    C = alpha * (
                4 * T_pll * betha * omega_n ** 2 * xi ** 2 - T_pll * betha * omega_n ** 2 - 2 * T_pll * omega_n ** 3 * xi - 2 * betha * omega_n * xi + omega_n ** 2) / (
                    1 - 2 * xi * omega_n * T_pll + omega_n ** 2 * T_pll ** 2)

    # A
    A_val = -RoCoF * (2 * H + D_damping * T_pll)

    # D
    D_val = RoCoF * D_damping * omega_n ** 2 / (omega_n ** 2 + 1 / T_pll ** 2 - 2 * xi * omega_n / T_pll)

    # B
    B_val = -A - (D / T_pll)

    # C
    C_val = -RoCoF * 2 * H * T_pll * omega_n ** 2 \
            - A * (2 * xi * omega_n + omega_n ** 2 * T_pll) \
            - D * omega_n ** 2

    alpha1 = omega_n * (xi_over + np.sqrt(xi_over ** 2 - 1))
    alpha2 = omega_n * (xi_over - np.sqrt(xi_over ** 2 - 1))
    term1 = (B * alpha1 - C) * np.exp(-alpha1 * t) / (alpha1 - alpha2)
    term2 = (B * alpha2 - C) * np.exp(-alpha2 * t) / (alpha1 - alpha2)
    DeltaP = A + term1 - term2 + D*np.exp(-t/T_pll)/T_pll
    #Ppeak=DeltaP[-1] # In this case DeltaP steady state will be calculated in another function "GetDeltaP"
    return DeltaP,xi

#Here we are going to form DeltaP considering the RoCof that stops increasing or decreasing in order to see that P comes back to a steady state value
def GetDeltaP(D_Damping,H,Xtotal_initial,P0):



    DeltaP = np.zeros_like(t_DeltaP)
    DeltaP, xi = GetDeltaP_NotDELAYED(D_Damping, H, Xtotal_initial, P0, t_DeltaP)  # <-- your actual function here



    DeltaP=np.where(Event_Time<t_DeltaP ,DeltaP,0)
    DeltaPSteadyState = DeltaP[-1] #Ppeak is the


    DeltaP_Recovered = GetDeltaP_NotDELAYED(D_Damping, H, Xtotal_initial, P0, t_DeltaP-RoCofStop_Time)[0]*-1
    DeltaP_Recovered = np.where(RoCofStop_Time < t_DeltaP, DeltaP_Recovered, 0)

    DeltaP= DeltaP + DeltaP_Recovered
    DeltaP=DeltaP*-1

    #Deactivate only to see the analitical response of Delta¨P
    # Create the plot
    #plt.figure(figsize=(8, 5))  # Set figure size
    #plt.plot(t_DeltaP, DeltaP, label="Ppcc", color='black', linestyle='--')  # First plot


    print("PsteadyState",DeltaPSteadyState)
    return DeltaP,DeltaPSteadyState,xi

#gets a tolerance band equals to the maximum value among 0.02Pn and 5% of DeltaPSteadyState (DeltaP in steady state)
def GetTunnel(DeltaP):
    Tunnel = max(0.02, 0.05 * DeltaP)
    print("Tunnel=",Tunnel)
    return Tunnel

#Cuts a signal between a min and a max value
def Cutsignal(Valuemin,Signal,Valuemax):
    Signal = np.where(Signal < Valuemin, Valuemin, Signal)
    Signal = np.where(Signal > Valuemax, Valuemax, Signal)
    #print("Value Min:", Valuemin)
    #print("Value Max:", Valuemax,"-")
    return Signal




def GetEnvelops(MargeUp,MargeDown,Signal,Tunnel):


    # Creating Envelops
    if RoCoF <= 0:
        print(' RoCoF <= 0 ')

        ########################upper_envelope#########################
        margin=AddMargin(MargeUp, Tunnel, RoCofStop_Time)
        # Envelope curves
        upper_envelope = Signal + margin+P0 #adding margin to the envelop
        # Find value of upper_envelope at RoCofStop_Time and use this value to be kept in the mask range time
        mask = (t_DeltaP >= RoCofStop_Time) & (t_DeltaP <= End_Time)
        upper_envelope = KeepTheValueatSpecificTimeUpper(upper_envelope, RoCofStop_Time, mask)


        ########################lower_envelope#########################
        margin = AddMargin(MargeDown, Tunnel, RoCofStop_Time)
        # Envelope curves
        lower_envelope = Signal - margin+P0 #adding margin to the envelop
        # Find value "P0 - Tunnel" and use this value to be kept in the mask range time
        mask = (t_DeltaP >= Event_Time) & (t_DeltaP <= RoCofStop_Time)
        lower_envelope = KeepTheValueatSpecificTimeLower(lower_envelope, Event_Time, mask)



        #condition = mask & (lower_envelope < (P0 - Tunnel))
        #lower_envelope = np.where(condition, P0 - Tunnel, lower_envelope)



    else:
        print(' RoCoF > 0 ')
        ########################upper_envelope#########################
        margin=AddMargin(MargeUp, Tunnel, RoCofStop_Time)
        # Envelope curves
        upper_envelope = Signal + margin+P0 #adding margin to the envelop
        # Find value "P0 + Tunnel" and use this value to be kept in the mask range time
        mask = (t_DeltaP >= Start_Time) & (t_DeltaP <= RoCofStop_Time)
        upper_envelope = KeepTheValueatSpecificTimeUpper(upper_envelope, Event_Time, mask)




        #mask = (t_DeltaP >= Start_Time) & (t_DeltaP <= RoCofStop_Time)
        #condition = mask & (upper_envelope > (P0 + Tunnel))
        #upper_envelope = np.where(mask & (upper_envelope > (P0 + Tunnel)), P0 + Tunnel, upper_envelope)


        #Use to debug conditions
        #print(upper_envelope)
        #print(np.any(condition))  # Should be True if any value exceeds P0 + Tunnel during the window
        #print(np.max(upper_envelope[mask]))  # Check actual values
        #print(P0 + Tunnel)  # Compare to see if condition could ever be True

        ########################lower_envelope#########################

        margin=AddMargin(MargeDown, Tunnel, RoCofStop_Time)
        #Envelope curves
        lower_envelope = Signal - margin + P0 #adding margin to the envelop
        # Find value of lower_envelope at RoCofStop_Time and use this value to be kept in the mask range time
        mask = (t_DeltaP >= RoCofStop_Time) & (t_DeltaP <= End_Time)
        lower_envelope = KeepTheValueatSpecificTimeLower(lower_envelope, RoCofStop_Time, mask)


    # Putting a limit to the active power "Signal DOWN" in case of OverCurrent
    mask = (t_DeltaP >= Event_Time) & (t_DeltaP <= End_Time)
    condition = mask & (lower_envelope > (Pmax_MoisTunnel))
    lower_envelope = np.where(condition, Pmax_MoisTunnel, lower_envelope)

    # Putting a limit to the active power "SIGNAL UP" in case of OverCurrent

    mask = (t_DeltaP >= Event_Time) & (t_DeltaP <= End_Time)
    condition = mask & (upper_envelope < (Pmin_MoisTunnel))
    upper_envelope = np.where(condition, Pmin_MoisTunnel, upper_envelope)

    upper_envelope = Cutsignal(Pmin_,upper_envelope,Pmax_)
    lower_envelope = Cutsignal(Pmin_, lower_envelope, Pmax_)
    #print(type(Signal_up_anal),"the type of up_anal")
    return upper_envelope,lower_envelope

def AddMargin(MargeUp,Tunnel,RoCofStop_Time):
     initial_margin = MargeUp
     final_margin = Tunnel
     decay_rate = 3  # tune this to control how fast the margin narrows
     RoCofStop_Time = RoCofStop_Time
     #margin = final_margin + (initial_margin - final_margin) * np.exp(-decay_rate * t_DeltaP) + (0.1) * np.exp(-decay_rate * (t_DeltaP - 2))
     margin = initial_margin * ((RoCofStop_Time>=t_DeltaP) &(t_DeltaP>= Event_Time)) * np.exp(-decay_rate * t_DeltaP) + initial_margin * (t_DeltaP >= RoCofStop_Time) * np.exp(-decay_rate * (t_DeltaP - RoCofStop_Time)) + Tunnel
     # Deactivate only to see the analitical response of Delta¨P
     # Create the plot
     #plt.figure(figsize=(8, 5))  # Set figure size
     #plt.plot(t_DeltaP, margin, label="Ppcc", color='black', linestyle='--')  # First plot
     return margin


def GetValueatSpecificTime(SelectedTime,Signal):

    index = np.argmin(np.abs(t_DeltaP - (SelectedTime - 0.01)))  # taking the value of P 10ms before RoCofStop_Time
    # Get value from the signal
    value_at_RoCofStop_Time = Signal[index]
    return  value_at_RoCofStop_Time

def KeepTheValueatSpecificTimeUpper(Signal,SpecificTime,maskTime):
    # Find value of upper_envelope at RoCofStop_Time
    Signal_value_at_SpecificTime_Time = GetValueatSpecificTime(SpecificTime, Signal)

    print("Signal_value_at_SpecificTime_Time", str(Signal_value_at_SpecificTime_Time))
    #mask = (t_DeltaP >= RoCofStop_Time) & (t_DeltaP <= End_Time)
    condition = maskTime & (Signal > Signal_value_at_SpecificTime_Time)
    Signal = np.where(condition, Signal_value_at_SpecificTime_Time, Signal)
    return Signal

def KeepTheValueatSpecificTimeLower(Signal,SpecificTime,maskTime):
    # Find value of upper_envelope at RoCofStop_Time
    Signal_value_at_SpecificTime_Time = GetValueatSpecificTime(SpecificTime, Signal)

    print("Signal_value_at_SpecificTime_Time", str(Signal_value_at_SpecificTime_Time))
    #mask = (t_DeltaP >= RoCofStop_Time) & (t_DeltaP <= End_Time)
    condition = maskTime & (Signal < Signal_value_at_SpecificTime_Time)
    Signal = np.where(condition, Signal_value_at_SpecificTime_Time, Signal)
    return Signal

def DelayEnvelops(P_up_finale,P_down_finale,P_PCC,shift_Time):
    TimeTODelay_All_Signals = shift_Time  # ms

    P_up_finale = delay_signal(TimeTODelay_All_Signals, fs, P_up_finale)
    P_down_finale = delay_signal(TimeTODelay_All_Signals, fs, P_down_finale)
    P_PCC = delay_signal(TimeTODelay_All_Signals, fs, P_PCC)


    return P_up_finale,P_down_finale,P_PCC

RoCoF = -2.5/50  # Rate of Change of Frequency (Hz/s) ou pu ?
H = 2.2       # Inertia constant (s)
T_pll = 0.01    # PLL time constant (s)
SCR=2
D_damping=165#Damping constant of the VSM control
wb=314 # Base angular frequency(rad/s)
xtr=0.15 #Transformer reactance (pu)
Ugrid=1 # RMS voltage Ugrid (pu)
Uconv=1 # RMS voltage Uconverter (pu)
Xeff=0.25 # effective reactance (pu)
EMT= True # Can be "True" or "False" EMT is activated (20ms for the measures)
P0= 0.5 # Initial power (pu)
Pmax_=1.2 #Pmax
Pmin_=-1.2 #Pmin
Pmax_MoisTunnel= Pmax_*0.95 #Considered for current limitation
Pmin_MoisTunnel=Pmin_*0.95 #Considered for current limitation

Z_grid=1/SCR

print("Final DeltaP",RoCoF*(2*H+D_damping*T_pll))


#second Order system

# Define the time vector for simulation
Start_Time = -1 # sec
Event_Time=0 #keep this value to "0"
RoCofStop_Time= 0.5 # sec
End_Time = 8 # sec



NbPoints = 10000
t_DeltaP = np.linspace(Start_Time, End_Time, NbPoints)  # From Start_Time to End_Time
fs = (End_Time - Start_Time) / NbPoints  # Sampling frequency (Hz)

# Calculating VARIABLES that need to be defined to calculate DeltaP
Xtotal = Xeff + Z_grid  # X total initial that is equal to Xeff+Xgrid inital

#Defining margins for H and D

Ratio_H_D_UP = 0.1 # Use to have two more values for D and H: D*(1+Ratio_H_D_UP), H*(1+Ratio_H_D_UP)
Ratio_H_D_Down = 0.1 # Use to have two more values for D and H: D*(1-Ratio_H_D_Down), H*(1-Ratio_H_D_Down)

# Defining arrays to consider DeltaP for different H and D
DeltaP_array = []
DeltaPSteadyState_array = []
Tunnel_array = []
Epsilon_array = []
P_up_anal_array = []
P_down_anal_array = []

#3 different values for D and H are considered
D_array=[D_damping,D_damping*(1+Ratio_H_D_UP),D_damping*(1-Ratio_H_D_Down)]
H_array=[H,H*(1-Ratio_H_D_UP),H*(1+Ratio_H_D_Down)]
print("Set of D values to be considered:",D_array)
print("Set of H values to be considered", H_array)

#Retrieving the second order response and the Tunnel that will be used in the Margins

results = [GetDeltaP(D_array[i], H_array[i], Xtotal, P0) for i in range(len(D_array))]
DeltaP_array, DeltaPSteadyState_array , Epsilon_array= map(np.array, zip(*results))
Tunnel_array = [GetTunnel(DeltaPSteadyState_array[i]) for i in range(len(D_array))]


#Creating Envelops
MargeUp=0.2 # This is the Margin up used in an exponential function around DeltaP
MargeDown=0.2 # This is the Margin down used in an exponential function around DeltaP
DeltaP = DeltaP_array[0]
Tunnel = Tunnel_array[0]
epsilon = Epsilon_array[0]

print("Set of Tunnel values to be considered", Tunnel_array)
print("Set of epsilon values to be considered", Epsilon_array)


results = [GetEnvelops(MargeUp,MargeDown,DeltaP_array[i],Tunnel_array[i]) for i in range(len(D_array))]
P_up_anal_array, P_down_anal_array = map(np.array, zip(*results))

#Theoretical Value is delimited
P_PCC=Cutsignal(Pmin_,P0+DeltaP,Pmax_)


# Stack the arrays into a 2D NumPy array
stacked = np.vstack((P_up_anal_array ))
# Compute the element-wise max, ignoring NaNs
#P up is created from the maximum values of P_up arrays
P_up_finale = np.nanmax(stacked, axis=0)

# Stack the arrays into a 2D NumPy array
stacked = np.vstack((P_down_anal_array ))
# Compute the element-wise max, ignoring NaNs
#P down is created from the minimum values of P_down arrays
P_down_finale = np.nanmin(stacked, axis=0)



#Adding delays when the simulation is done in EMT
if EMT:
    # Delay settings
    delayEMT_ms = 20 # 20 ms of delay for EMT simulations
    shift_Time = 0
    # Pad with the initial value instead of zero
    initial_value = P_up_finale[0]
    P_up_finale = delay_signal(delayEMT_ms, fs, P_up_finale)
    P_down_finale = delay_signal(delayEMT_ms, fs, P_down_finale)
    P_PCC = delay_signal(delayEMT_ms, fs, P_PCC)

    #Delay P_up_anal_array
    results = [delay_signal(delayEMT_ms,fs, P_up_anal_array[i]) for i in range(len(D_array))]
    P_up_anal_array = results
    print("size P_up_anal_array:", len(P_up_anal_array))
    # Delay P_down_anal_array
    results = [delay_signal(delayEMT_ms,fs, P_down_anal_array[i]) for i in range(len(D_array))]
    P_down_anal_array = results

else:
    shift_Time = 0
    delay_samples=0
    initial_value = P_up_finale[0]


##############plot PCC from Open Modelica

# Path to the CSV file
BaseLocation= "RMSsimulations/"

#OverDAMPED
#csv_file_path_Gabarits=BaseLocation + "gabarit_overdamped.csv"
#csv_file_path_OM=BaseLocation + "OM_DeltaP_OverDampedSCR2H3D109Angle3.6.csv"
#csv_file_path_OM=BaseLocation +"H=2.2,D=133,Xeff=0.25,Imax=1.2,P0=0.5,SCR=20,Imax=1.2.csv"
csv_file_path_OM=BaseLocation +"P0=0.5,H=2.2,D=133,Xeff=0.25,Imax=1.2,P0=0.5,SCR=20,Imax=1.2.csv"
#Name of the Columns
NameColumnsDataFrame = ["Time","Pup","Pdown"]
# Read the CSV file into a DataFrame
#dataUseCase_Gabarits = pd.read_csv(csv_file_path_Gabarits,sep=";")

dataUseCase_OM = pd.read_csv(csv_file_path_OM)



P_Pcc="gFM_VSM_cc.measurementPcc.PGenPu"

#filtering times
TimeInit=9
TimeFinal=11

TimeInit=14
TimeFinal=19

#getting time
t=dataUseCase_OM['time']
y=dataUseCase_OM[P_Pcc]

# Extract data from t = 5 to t = 10
mask = (t >= TimeInit) & (t <= TimeFinal)
t_selected = t[mask]  # Time values in range 10-15
y_selected = y[mask]  # Corresponding function values



# Shift time so that it starts at t = 0
t_shifted = t_selected -t_selected.iloc[0]  # Subtract the first value to start from 0
#print(t_shifted)

#filtering over a time range
filtered_dataUseCase_OM = dataUseCase_OM[(dataUseCase_OM['time'] >= TimeInit) & (dataUseCase_OM['time'] <= TimeFinal)]

#Taking the axis X
axisX = filtered_dataUseCase_OM['time']



# Create the plot
Title= "P0="+ str(P0) +", RoCoF=" + str(RoCoF) +  ", Epsilon= " + str(round(epsilon,3)) + ", ωd= " + ", D= " + str(D_damping) + ", H= " +str(H) + ", Xeff= " + str(Xeff)+ ", EMT=" +str(EMT)

plt.figure(figsize=(8, 5))  # Set figure size
plt.plot(t_shifted-1, y_selected, label="P_pcc from Open Modelica", color='b', linestyle='-')  # First plot

plt.plot(t_DeltaP+shift_Time,P_PCC, label="Ppcc", color='black', linestyle='--')  # First plot
plt.plot(t_DeltaP+shift_Time,P_down_anal_array[0], label="Pdown_analytical", color='m', linestyle='--')  # First plot

plt.plot(t_DeltaP+shift_Time,P_up_anal_array[0], label="Pup_analytical", color='m', linestyle='-')  # First plot
plt.plot(t_DeltaP+shift_Time,P_down_anal_array[1], label="Pdown_analytical", color='r', linestyle='--')  # First plot
plt.plot(t_DeltaP+shift_Time,P_up_anal_array[1], label="Pup_analytical", color='r', linestyle='-')  # First plot
plt.plot(t_DeltaP+shift_Time,P_down_anal_array[2], label="Pdown_analytical", color='g', linestyle='--')  # First plot
plt.plot(t_DeltaP+shift_Time,P_up_anal_array[2], label="Pup_analytical", color='g', linestyle='-')  # First plot

#plt.plot(t_DeltaP+shift_Time,P_50Prc, label="P_50%", linewidth='3')  # First plot
#plt.plot(t_DeltaP+shift_Time,P_down_finale, label="Pdown_final", linewidth=2)  # First plot
#plt.plot(t_DeltaP+shift_Time,P_up_finale, label="Pup_final", linewidth=2)  # First plot


# Add vertical line at t = 20 ms
plt.axvline(x=0.020, color='black', linestyle='--', label='t = 20ms')

# Add labels, title, and legend
plt.xlabel("sec")
plt.ylabel("P at PCC (pu)")
plt.title(Title)
plt.legend()  # Show legend

# Show the plot
plt.grid(True)  # Add grid for better visualization


#Delay the signals that will be save in a csv file
shift_Time = 0 #ms
P_up_finale,P_down_finale,P_PCC = DelayEnvelops(P_up_finale,P_down_finale,P_PCC,shift_Time)

# Save to CSV
LocationFile=  "AnalyticalEnvelops/"+Title +".csv"
# Combine into a DataFrame with custom column names
df = pd.DataFrame({
    "Time (s)": t_DeltaP,
    "P_PCC (pu)": P_PCC,
    "P_down (pu)": P_down_finale,
    "P_up (pu)": P_up_finale
})

# Export to CSV

df.to_csv(LocationFile, index=False)

# Create the plot
plt.figure(figsize=(8, 5))  # Set figure size
plt.plot(t_DeltaP,P_PCC, label="P_PCC analytical", linewidth='3')  # First plot
plt.plot(t_DeltaP,P_down_finale, label="Pdown_final", linewidth=2)  # First plot
plt.plot(t_DeltaP,P_up_finale, label="Pup_final", linewidth=2)  # First plot
plt.plot(t_shifted-1, y_selected, label="P_pcc from Open Modelica", color='b', linestyle='-')  # First plot

# Add labels, title, and legend
plt.xlabel("sec")
plt.ylabel("P at PCC (pu)")
plt.title(Title)
# Add vertical line at t = x seconds
plt.axvline(x=Event_Time+shift_Time/1000, color='black', linestyle='--', label='t at Event Time')
plt.legend(loc='lower right')  # Show legend
# Show the plot
plt.grid(True)  # Add grid for better visualization

# Save the figure with specific size and resolution
Path = "RMSsimulations/PNGResults/"+Title + ".png"
plt.savefig(Path, bbox_inches='tight', dpi=300)

plt.show()



