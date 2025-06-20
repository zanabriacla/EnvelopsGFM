
import pandas as pd
import scipy.io
import matplotlib.pyplot as plt
import math
import numpy as np





def delay_signal(delay_ms,fs,signal):
    delay_samples = int((delay_ms / 1000) * 1 / fs)
    signal_delayed = np.concatenate((np.full(delay_samples, signal[0]), signal))[:len(signal)]
    return signal_delayed

def GetDeltaP(D,H,Xtotal_initial,P0):
    print("Second Order Function UNDER DAMPED")
    # second Order system
    # Calculating the epsilon and Wn of the second order ystsem

    wn = math.sqrt(wb * Uconv * Ugrid / (2 * H * Xtotal_initial))
    epsilon = D / (4 * H * wn)
    print("epsilon", epsilon)
    wd = wn * math.sqrt(1 - epsilon ** 2)

    print("epsilon: ", epsilon, "natural frequency: ", wn, "Damped natural frequency: ", wd)

    # Assuming DeltaAngleGrid, Uconv, Ugrid, and Xtotal are already defined
    Ppeak = DeltaXtotal * P0 / (Xtotal_initial)

    print("Delta X total", DeltaXtotal)
    print("peak power change", Ppeak)
    # Define the time vector for simulation
    #t_DeltaP = np.linspace(0, End_Time, 10000)  # From 0 to 2 seconds

    # Assuming epsilon, wn, wd, t, and Ppeak are already defined as NumPy arrays or scalars
    DeltaP = Ppeak * (
                np.exp(-epsilon * wn * t_DeltaP) * np.cos(wd * t_DeltaP) + (D / (2 * H) - epsilon * wn) / wd * np.exp(
            -epsilon * wn * t_DeltaP) * np.sin(wd * t_DeltaP)) * -1
    P = P0 + DeltaP

    numerator = D - 2 * H * epsilon * wn
    denominator = 2 * H * wd
    A = np.sqrt(1 + (numerator / denominator) ** 2)

    AmplitudeEnvelops = np.sqrt(1 + ((D - 2 * H * epsilon * wn) / (2 * H * wd)) ** 2)

    DeltaP_up_anal_array =  np.abs(AmplitudeEnvelops * Ppeak * np.exp(-epsilon * wn * t_DeltaP))
    DeltaP_down_anal_array = np.abs(AmplitudeEnvelops * Ppeak * np.exp(-epsilon * wn * t_DeltaP))*-1

    #Before Event the values of the signals are set to "0"
    DeltaP = np.where(t_DeltaP < Event_Time, 0, DeltaP)
    DeltaP_up_anal_array = np.where(t_DeltaP < Event_Time, 0, DeltaP_up_anal_array)
    DeltaP_down_anal_array = np.where(t_DeltaP < Event_Time, 0, DeltaP_down_anal_array)
    return DeltaP,Ppeak,epsilon, DeltaP_up_anal_array, DeltaP_down_anal_array

def GetTunnel(Ppeak):
    Tunnel = max(0.02, 0.05 * Ppeak) # For the tunnel we need to take the max between 0.02In and 5%DeltaP (at
    return Tunnel

def Cutsignal(Valuemin,Signal,Valuemax):
    Signal = np.where(Signal < Valuemin, Valuemin, Signal)
    Signal = np.where(Signal > Valuemax, Valuemax, Signal)
    print("Value Min:", Valuemin)
    print("Value Max:======", Valuemax,"-")
    return Signal

def GetEnvelops(MargeUp,MargeDown,Signal,Tunnel,DeltaP):

    # Creating Envelops
    if DeltaPAtEventTime > 0:
        print("DeltaP>0")
        Signal_up_anal = Signal * (1 + MargeUp) + Tunnel + P0
        Signal_down_anal = Signal * (1 - MargeDown) - Tunnel + P0

        # Putting a limit to the active power "Signal DOWN"
        mask = (t_DeltaP >= Event_Time) & (t_DeltaP <= End_Time)
        condition = mask & (Signal_down_anal > (Pmax_MoisTunnel))
        Signal_down_anal = np.where(condition, Pmax_MoisTunnel, Signal_down_anal)


    else:
        print("DeltaP<=0")
        #Signal_up_anal = Signal * (1 - MargeUp) + Tunnel + P0

        Signal_up_anal = Signal * (1 - MargeUp) + Tunnel + P0

        # Putting a limit to the active power "SIGNAL UP"

        mask = (t_DeltaP >= Event_Time) & (t_DeltaP <= End_Time)
        condition = mask & (Signal_up_anal  < (Pmin_MoisTunnel))
        Signal_up_anal = np.where(condition, Pmin_MoisTunnel, Signal_up_anal)


        Signal_down_anal = Signal * (1 + MargeDown) - Tunnel + P0

    Signal_up_anal = Cutsignal(Pmin_,Signal_up_anal,Pmax_)
    Signal_down_anal = Cutsignal(Pmin_, Signal_down_anal, Pmax_)
    print(type(Signal_up_anal),"the type of up_anal")
    return Signal_up_anal,Signal_down_anal

def LimitingReversePower(P_up_finale,P_down_finale, P0,Tunnel,DeltaP):

        # Creating Envelops
    if (DeltaPAtEventTime) > 0:
        print("DeltaP>0")

        # Putting a limit to the active power only for 100ms
        mask = (t_DeltaP >= Event_Time) & (t_DeltaP <= Event_Time+0.1)
        condition = mask & (P_down_finale < (P0-Tunnel))
        P_down_finale = np.where(condition, P0-Tunnel, P_down_finale)


    else:
        print("DeltaP<=0")
        # Signal_up_anal = Signal * (1 - MargeUp) + Tunnel + P0
        print("Event_Time",Event_Time)

        # Putting a limit to the active power only for 100ms
        mask = (t_DeltaP >= Event_Time) & (t_DeltaP <= Event_Time+0.1)
        condition = mask & (P_up_finale > (P0 + Tunnel))
        P_up_finale = np.where(condition, P0 + Tunnel, P_up_finale)
        # Create the plot
        plt.figure(figsize=(8, 5))  # Set figure size
        plt.plot(t_DeltaP - 1, P_up_finale, label="P_pcc from Open Modelica", color='b',
                 linestyle='-')  # Adding simulation

        return P_up_finale,P_down_finale

def DelayEnvelops(P_up_finale,P_down_finale,P_PCC,shift_Time):
    TimeTODelay_All_Signals = shift_Time  # ms
    TimeTODelay_LowerBoundATBeggining = 10  # ms

    P_up_finale = delay_signal(TimeTODelay_All_Signals, fs, P_up_finale)
    P_down_finale = delay_signal(TimeTODelay_All_Signals, fs, P_down_finale)

    P_up_finale = np.where(t_DeltaP < Event_Time, P0 + Tunnel, P_up_finale)
    P_down_finale = np.where(t_DeltaP < Event_Time, P0 - Tunnel, P_down_finale)

    if (P0 > 0 and DeltaPAtEventTime>0) or (P0 < 0 and DeltaPAtEventTime>0):



        #P down needs to be delayed even more
        P_down_finale = delay_signal(TimeTODelay_LowerBoundATBeggining, fs, P_down_finale)


    else:

        # P up needs to be delayed even more
        P_up_finale = delay_signal(TimeTODelay_LowerBoundATBeggining, fs, P_up_finale)

    P_PCC = delay_signal(TimeTODelay_All_Signals, fs, P_PCC)
    P_PCC = np.where(t_DeltaP < Event_Time, P0, P_PCC)

    return P_up_finale,P_down_finale,P_PCC
def GetValueatSpecificTime(SelectedTime,Signal):

    index = np.argmin(np.abs(t_DeltaP - (SelectedTime - 0.01)))  # taking the value of P 10ms before RoCofStop_Time
    # Get value from the signal
    value_at_RoCofStop_Time = Signal[index]
    return  value_at_RoCofStop_Time

#Variables needed to be fulfilled in order to implement the envelops

SCR_init=2 #SCR ini
SCR_final=10 #SCR final

Z_init=1/SCR_init
Z_final=1/SCR_final


Delta_ZGrid = Z_final-Z_init #DeltaZgrid

print("DeltaZgrid",Delta_ZGrid)


D=220#Damping constant of the VSM control
H=7 #Inertia constant (s)
wb=314 # Base angular frequency(rad/s)
xtr=0.06 #Transformer reactance (pu)
Ugrid=1 # RMS voltage Ugrid (pu)
Uconv=1 # RMS voltage Uconverter (pu)
Xeff=0.25 # effective reactance (pu)
EMT= True # Can be "True" or "False" EMT is activated (20ms for the measures)
P0= 0.5 # Initial power (pu)
Pmax_=1.2 #Pmax
Pmax_MoisTunnel= Pmax_*0.95 #Pmax
Pmin_=-1.2 #Pmin
Pmin_MoisTunnel=Pmin_*0.95




'''
if EventOnZgrid == "up":
    Delta_ZGrid = EventOnZgrid_up(Z1,Z2)[0]
elif EventOnZgrid == "down":
    Delta_ZGrid = EventOnZgrid_down(Z1, Z2)[0]
else:
    print("Invalid mode")

'''


#second Order system

# Define the time vector for simulation


Start_Time = -1
Event_Time = 0
End_Time = 4
NbPoints = 10000
t_DeltaP = np.linspace(Start_Time, End_Time, NbPoints)  # From Start_Time to End_Time
fs = (End_Time - Start_Time) / NbPoints  # Sampling frequency (Hz)


# Calculating VARIABLES that need to be defined to calculate DeltaP
Xtotal_initial = Xeff + Z_final  # X total initial that is equal to Xeff+Xgrid final
DeltaXtotal = Delta_ZGrid  # Variation in Xtotal

#Defining margins for H and D

Ratio_H_D_UP = 0.2 # Use to have two more values for D and H: D*(1+Ratio_H_D_UP), H*(1+Ratio_H_D_UP)
Ratio_H_D_Down = 0.2 # Use to have two more values for D and H: D*(1-Ratio_H_D_Down), H*(1-Ratio_H_D_Down)

# Defining arrays to consider DeltaP for different H and D
DeltaP_array = []
Ppeak_array = []
Tunnel_array = []
Epsilon_array = []
P_up_anal_array = []
P_down_anal_array = []

D_array=[D,D*(1+Ratio_H_D_UP),D*(1-Ratio_H_D_Down)]
H_array=[H,H*(1-Ratio_H_D_UP),H*(1+Ratio_H_D_Down)]
print("Set of D values to be considered:",D_array)
print("Set of H values to be considered", H_array)

#Retrieving the second order response and the Tunnel that will be used in the Margins

results = [GetDeltaP(D_array[i], H_array[i], Xtotal_initial, P0) for i in range(len(D_array))]
DeltaP_array, Ppeak_array , Epsilon_array, DeltaPSecond_up_anal_array, DeltaPSecond_down_anal_array= map(np.array, zip(*results))
Tunnel_array = [GetTunnel(Ppeak_array[i]) for i in range(len(D_array))]

#Creating Envelops
MargeUp=0.2 # This is the Margin up used in DeltaP*(1+-MarginUp)+Tunnel
MargeDown=0.2 # This is the Margin down used in DeltaP*(1+-MargeDown)+Tunnel
DeltaP = DeltaP_array[0]
Tunnel = Tunnel_array[0]
epsilon = Epsilon_array[0]

DeltaPAtEventTime = GetValueatSpecificTime(Event_Time+0.01,DeltaP)

results = [GetEnvelops(MargeUp,MargeDown,DeltaP_array[i],Tunnel_array[i],DeltaP) for i in range(len(D_array))]
results_PSecond_up_anal_array= [GetEnvelops(MargeUp,MargeDown,DeltaPSecond_up_anal_array[i],Tunnel_array[i],DeltaP) for i in range(len(D_array))]
results_PSecond_down_anal_array= [GetEnvelops(MargeUp,MargeDown,DeltaPSecond_down_anal_array[i],Tunnel_array[i],DeltaP) for i in range(len(D_array))]

#Second order limit
P_up_anal_array, P_down_anal_array = map(np.array, zip(*results))
#First order limit for P up
PSecond_up_up_anal_array, PSecond_up_down_anal_array = map(np.array, zip(*results_PSecond_up_anal_array))
#First order limit for P down
PSecond_down_up_anal_array, PSecond_down_down_anal_array  = map(np.array, zip(*results_PSecond_down_anal_array))

#Theoretical Value
#P_PCC= P0+DeltaP
P_PCC=Cutsignal(Pmin_,P0+DeltaP,Pmax_)
#P_PCC = np.where(P_PCC < -1, -1, P_PCC)



#Envelop of 50% of Delta before t=10ms and after that it takes DeltaP
P_50Prc = P0+ np.where(t_DeltaP >= Start_Time, DeltaP*0.5 , DeltaP)
P_50Prc=Cutsignal(Pmin_,P_50Prc,Pmax_)

#From different possibilities , we select the max value for Envelop UP and the MIN value for envelop DOWN
#P_up_finale = np.maximum(P_up_anal_array[0], P_up_anal_array[1],P_50Prc)

# Element-wise max across all


#P_up_finale = np.maximum.reduce(P_up_anal_array + [P_50Prc])
#P_down_finale = np.minimum.reduce(P_down_anal_array + [P_50Prc])


# Stack the arrays into a 2D NumPy array
#stacked = np.vstack((P_up_anal_array , [P_50Prc], P_up_anal_array, Cutsignal(Pmin_,PSecond_up_anal_array,Pmax_), Cutsignal(Pmin_,PSecond_down_anal_array,Pmax_) ))
stacked = np.vstack(( P_up_anal_array, [P_50Prc] ))
# Compute the element-wise max, ignoring NaNs
P_up_finale = np.nanmax(stacked, axis=0)


# Stack the arrays into a 2D NumPy array
#stacked = np.vstack((P_down_anal_array , [P_50Prc], P_up_anal_array, PSecond_up_anal_array, PSecond_down_anal_array ))
stacked = np.vstack((P_down_anal_array, [P_50Prc]))
# Compute the element-wise max, ignoring NaNs
P_down_finale = np.nanmin(stacked, axis=0)

print("P_50Prc",P_50Prc)



#cutting Pdown finale and Pup finale in order to avoir reverse power

#P_up_finale ,  P_down_finale = LimitingReversePower(P_up_finale,P_down_finale, P0,Tunnel,DeltaP)

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
    P_50Prc = delay_signal(delay_ms, fs, P_50Prc)


    results = [delay_signal(delay_ms,fs, P_up_anal_array[i]) for i in range(len(D_array))]
    P_up_anal_array = results
    print("size P_up_anal_array:", len(P_up_anal_array))
    results = [delay_signal(delay_ms,fs, P_down_anal_array[i]) for i in range(len(D_array))]
    #P_down_anal_array = map(np.array, zip(*results))
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
#csv_file_path_OM=BaseLocation +"H=10,D=133,Xeff=0.25,Imax=1.2,P0=0.5,SCRini=4,SCRfinal=2,Imax=1.2.csv"
csv_file_path_OM=BaseLocation+"H=5,D=140,Xeff=0.06,Imax=1.2,P0=0.5,SCRini=2,SCRmax=4,Imax=1.2.csv"
#Name of the Columns
NameColumnsDataFrame = ["Time","Pup","Pdown"]
# Read the CSV file into a DataFrame
#dataUseCase_Gabarits = pd.read_csv(csv_file_path_Gabarits,sep=";")

dataUseCase_OM = pd.read_csv(csv_file_path_OM)



P_Pcc="gFM_VSM_cc.measurementPcc.PGenPu"

#filtering times
TimeInit=9
TimeFinal=11

TimeInit=4
TimeFinal=8

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




Title =  "P0="+ str(P0) +"pu, SCRinit=" + str(SCR_init) + ", SCRfinal= "+str(SCR_final) + ", Epsilon= " + str(round(epsilon,3))  +", D= " + str(D) + ", H= " +str(H) + ", Xeff= " + str(Xeff)+ ", Pmax="+ str(Pmax_) +"pu"+ ", EMT=" +str(EMT)









shift_Time= 0

# Create the plot
plt.figure(figsize=(8, 5))  # Set figure size
plt.plot(t_shifted-1, y_selected, label="P_pcc from Open Modelica", color='b', linestyle='-')  # Adding simulation
plt.plot(t_DeltaP+shift_Time,P_PCC, label="Ppcc", color='black', linestyle='--')  # First plot

plt.plot(t_DeltaP+shift_Time,P_down_anal_array[0], label="Pdown_analytical", color='m', linestyle='-')  # First plot
plt.plot(t_DeltaP+shift_Time,P_up_anal_array[0], label="Pup_analytical", color='m', linestyle='--')  # First plot
plt.plot(t_DeltaP+shift_Time,P_down_anal_array[1], label="Pdown_analytical", color='b', linestyle='-')  # First plot
plt.plot(t_DeltaP+shift_Time,P_up_anal_array[1], label="Pup_analytical", color='b', linestyle='--')  # First plot
plt.plot(t_DeltaP+shift_Time,P_down_anal_array[2], label="Pdown_analytical", color='g', linestyle='-')  # First plot
plt.plot(t_DeltaP+shift_Time,P_up_anal_array[2], label="Pup_analytical", color='g', linestyle='--')  # First plot

#Plottijg the first order responses envelops of the second order response
plt.plot(t_DeltaP+shift_Time,PSecond_up_up_anal_array[0], label="PSecond_up_anal_array", color='r', linestyle=':')  # First plot
plt.plot(t_DeltaP+shift_Time,PSecond_down_down_anal_array[0], label="PSecond_down_anal_array", color='r', linestyle='-')  # First plot
plt.plot(t_DeltaP+shift_Time,PSecond_up_up_anal_array[1], label="PSecond_up_anal_array", color='cyan', linestyle=':')  # First plot
plt.plot(t_DeltaP+shift_Time,PSecond_down_down_anal_array[1], label="PSecond_down_anal_array", color='cyan', linestyle='-')  # First plot
plt.plot(t_DeltaP+shift_Time,PSecond_up_up_anal_array[2], label="PSecond_up_anal_array", color='brown', linestyle=':')  # First plot
plt.plot(t_DeltaP+shift_Time,PSecond_down_down_anal_array[2], label="PSecond_down_anal_array", color='brown', linestyle='-')  # First plot


plt.plot(t_DeltaP+shift_Time,P_50Prc, label="P_50%", linewidth='3')  # First plot
#plt.plot(t_DeltaP+shift_Time,P_down_finale, label="Pdown_final", linewidth=2)  # First plot
#plt.plot(t_DeltaP+shift_Time,P_up_finale, label="Pup_final", linewidth=2)  # First plot


# Add vertical line at t = 10 seconds
plt.axvline(x=0.010, color='black', linestyle='--', label='t = 10ms')

LocationFile = "AnalyticalEnvelops/"+Title +".csv"

# Add labels, title, and legend
plt.xlabel("sec")
plt.ylabel("P at PCC (pu)")
plt.title(Title)
plt.legend(loc='lower right')  # Show legend

# Show the plot
plt.grid(True)  # Add grid for better visualization

#Delay the signals that will be save in a csv file
shift_Time=0 #ms
P_up_finale,P_down_finale,P_PCC = DelayEnvelops(P_up_finale,P_down_finale,P_PCC,shift_Time)



# Save to CSV
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
plt.plot(t_shifted-1, y_selected, label="P_pcc from Open Modelica", color='b', linestyle='-')  # Adding simulation
plt.plot(t_DeltaP,P_PCC, label="P_PCC analytical", linewidth='3')  # First plot
plt.plot(t_DeltaP,P_down_finale, label="Pdown_final", linewidth=2)  # First plot
plt.plot(t_DeltaP,P_up_finale, label="Pup_final", linewidth=2)  # First plot
# Add labels, title, and legend
plt.xlabel("sec")
plt.ylabel("P at PCC (pu)")
plt.title(Title)
plt.legend(loc='lower right')  # Show legend
# Add vertical line at event time
plt.axvline(x=Event_Time, color='black', linestyle='--', label='t at Event Time')

# Show the plot
plt.grid(True)  # Add grid for better visualization

# Save the figure with specific size and resolution
Path = "RMSsimulations/PNGResults/"+Title + ".png"
#plt.savefig(Path, bbox_inches='tight', dpi=300)

plt.show()

