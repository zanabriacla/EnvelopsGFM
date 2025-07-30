
import pandas as pd
import scipy.io
import matplotlib.pyplot as plt
import math
import numpy as np





def delay_signal(delay_ms,fs,signal):
    delay_samples = int((delay_ms / 1000) * 1 / fs)
    signal_delayed = np.concatenate((np.full(delay_samples, signal[0]), signal))[:len(signal)]
    return signal_delayed

#DeltaP fonction in case of SCR jump, we return also Ppeak, epsilon and two exponential up (DeltaP_up_anal_array) and low (DeltaP_down_anal_array)
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

#The tunnel to be considered at the end is either 0.02In or 5% of Ppeak
def GetTunnel(Ppeak):
    Tunnel = max(0.02, 0.05 * Ppeak) # For the tunnel we need to take the max between 0.02In and 5%DeltaP (at
    return Tunnel

#A Signal is cutted bettwen a max value and a min value
def Cutsignal(Valuemin,Signal,Valuemax):
    Signal = np.where(Signal < Valuemin, Valuemin, Signal)
    Signal = np.where(Signal > Valuemax, Valuemax, Signal)
    print("Value Min:", Valuemin)
    print("Value Max:======", Valuemax,"-")
    return Signal

# The upper and lower envelopes are defined here.
# For DeltaP > 0, the upper bound is calculated as: (DeltaP * (1 + MarginUp) + Tunnel + P0).
# The lower bound is calculated as: (DeltaP * (1 - MargeDown) - Tunnel + P0)
# The lower bound follows an exponential behavior, but this is not applied during the initial instants.
# The DeltaP input corresponds to DeltaP in the case of D and H base variants.
# It is important to use 50% of the DeltaP value at the beginning of the response to reflect the expected behavior.
# Note: 'Signal' refers to the DeltaP input for both D and H variants.

def GetEnvelops(MargeUp,MargeDown,Signal,Tunnel,DeltaP):

    # Creating Envelops
    if DeltaPAtEventTime > 0:
        print("DeltaP>0")
        # upper and lower bound built from the MargeUp and MargeDown values
        Signal_up_anal = (Signal * (1 + MargeUp) + Tunnel + P0)
        Signal_down_anal = (Signal * (1 - MargeDown) - Tunnel + P0)

        # Apply a limit to the active power "Signal DOWN" when it exceeds Pmax_Mois_Tunnel.
        # This ensures that the upper and lower bounds are not confused in cases of saturation,
        # and maintains a margin between the bounds for better distinction and stability.
        mask = (t_DeltaP >= Event_Time) & (t_DeltaP <= End_Time)
        condition = mask & (Signal_down_anal > (Pmax_MoisTunnel))
        Signal_down_anal = np.where(condition, Pmax_MoisTunnel, Signal_down_anal)

        #Modification of the signal down to consider DeltaP*50% at the beginning
        # Adding 0.005, If DeltaP = 0 and the power was previously at P = 0.02, the power is reduced back to 0.
        # This prevents an unintended power offset from being retained when no DeltaP is applied.

        P_50Prc = P0 + np.where(t_DeltaP >= Event_Time, DeltaP * 0.5+0.005, Signal)
        Signal_down_anal = EnvelopDowModify(Signal_down_anal, P_50Prc)


    else:
        print("DeltaP<=0")
        # upper and lower bound built from the MargeUp and MargeDown values
        Signal_up_anal = Signal * (1 - MargeUp) + Tunnel + P0
        Signal_down_anal = Signal * (1 + MargeDown) - Tunnel + P0

        # Apply a limit to the active power "Signal UP" when it is lower than Pmin_MoisTunnel.
        # This ensures that the upper and lower bounds are not confused in cases of saturation,
        # and maintains a margin between the bounds for better distinction and stability.

        mask = (t_DeltaP >= Event_Time) & (t_DeltaP <= End_Time)
        condition = mask & (Signal_up_anal  < (Pmin_MoisTunnel))
        Signal_up_anal = np.where(condition, Pmin_MoisTunnel, Signal_up_anal)



        # Modification of the signal up to consider DeltaP*50% at the beginning
        P_50Prc = P0 + np.where(t_DeltaP >= Start_Time, DeltaP * 0.5+0.005, Signal)
        Signal_up_anal = EnvelopDowModify(Signal_up_anal, P_50Prc)

    #Limiting upper bound to Pmin and Pmax
    #Limiting lower bound to Pmin and Pmax

    Signal_up_anal = Cutsignal(Pmin_,Signal_up_anal,Pmax_)
    Signal_down_anal = Cutsignal(Pmin_, Signal_down_anal, Pmax_)
    print(type(Signal_up_anal),"the type of up_anal")
    return Signal_up_anal,Signal_down_anal

# In the case of a positive DeltaP, we need to account for the fact that DeltaP increases during the initial moment.
# This is why we do not immediately apply the exponential decay, as it would not reflect the expected behavior.
# Instead, we hold DeltaP at 50% of its value for a duration defined by "DeltaTInitToKeepPeakP_50Prc",
# before transitioning to the lower exponential decay.

def EnvelopDowModify(Signal_Down,P_50Prc):
    DeltaTInitToKeepPeakP_50Prc = 0.03  # 50ms we kept the value of P_50Prc after that we consider the value of the PSecond_down_up_anal_array
    mask = (t_DeltaP >= Event_Time) & (t_DeltaP <= Event_Time+DeltaTInitToKeepPeakP_50Prc)
    Signal_Down = np.where(mask, P_50Prc, Signal_Down)
    #in case reaching the max value we consider that the envelop that is not saturated will be at the level saturated-0.2pu
    #this is applied only during the mask time, it means when we have cosidereded 50% of DeltaP
    Signal_Down = np.where(Signal_Down*mask < Pmin_, Pmin_ +0.2, Signal_Down)
    Signal_Down = np.where(Signal_Down*mask > Pmax_, Pmax_ - 0.2, Signal_Down)
    return Signal_Down

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

# At the beginning of the simulation, we need to account for measurement delays—even when simulations are performed in RMS mode.
# In this case, we apply a delay equal to "TimeTODelay_ATBeginningms".
# For EMT simulations, the total delay is "TimeTODelay_ATBeginningms" plus "delay_EMT_ms".
# In EMT mode, measurements follow the IEC 61400-21-1 standard, where a filter is applied to extract the positive sequence values,
# with a typical window duration of 20 ms.

def DelayEnvelops(P_up_finale,P_down_finale,P_PCC,shift_Time):
    TimeTODelay_All_Signals = shift_Time  # ms
    TimeTODelay_ATBegginingms = 10  # ms

    P_up_finale = delay_signal(TimeTODelay_All_Signals, fs, P_up_finale)
    P_down_finale = delay_signal(TimeTODelay_All_Signals, fs, P_down_finale)




    if (P0 > 0 and DeltaPAtEventTime>0) or (P0 < 0 and DeltaPAtEventTime>0):



        #P down needs to be delayed even more
        P_down_finale = delay_signal(TimeTODelay_ATBegginingms, fs, P_down_finale)

        if EMT:
            #P_up_finale = np.where(t_DeltaP < Event_Time+TimeTODelay_ATBegginingms/1000+delay_EMT_ms/1000, P0 + Tunnel, P_up_finale)
            P_down_finale = np.where(t_DeltaP < Event_Time+TimeTODelay_ATBegginingms/1000+delay_EMT_ms/1000, P0 - Tunnel, P_down_finale)

    else:

        # P up needs to be delayed even more
        P_up_finale = delay_signal(TimeTODelay_ATBegginingms, fs, P_up_finale)

        if EMT:
            P_up_finale = np.where(t_DeltaP < Event_Time+TimeTODelay_ATBegginingms/1000+delay_EMT_ms/1000, P0 + Tunnel, P_up_finale)
            #P_down_finale = np.where(t_DeltaP < Event_Time+TimeTODelay_ATBegginingms/1000+delay_EMT_ms/1000, P0 - Tunnel, P_down_finale)

    P_PCC = delay_signal(TimeTODelay_All_Signals, fs, P_PCC)
    P_PCC = np.where(t_DeltaP < Event_Time, P0, P_PCC)

    return P_up_finale,P_down_finale,P_PCC
def GetValueatSpecificTime(SelectedTime,Signal):

    index = np.argmin(np.abs(t_DeltaP - (SelectedTime - 0.01)))  # taking the value of P 10ms before RoCofStop_Time
    # Get value from the signal
    value_at_RoCofStop_Time = Signal[index]
    return  value_at_RoCofStop_Time

#Variables needed to be fulfilled in order to implement the envelops

SCR_init=10 #SCR ini
SCR_final=2 #SCR final

Z_init=1/SCR_init
Z_final=1/SCR_final

print(Z_init)
print(Z_final)

Delta_ZGrid = Z_final-Z_init #DeltaZgrid

print("DeltaZgrid",Delta_ZGrid)


D=80#Damping constant of the VSM control
H=3 #Inertia constant (s)
wb=314 # Base angular frequency(rad/s)
xtr=0.06 #Transformer reactance (pu)
Ugrid=1 # RMS voltage Ugrid (pu)
Uconv=1 # RMS voltage Uconverter (pu)
Xeff=0.25 # effective reactance (pu)
EMT= True # Can be "True" or "False" EMT is activated (20ms for the measures)
P0= 0.95 # Initial power (pu)
Pmax_=1.4 #Pmax
Pmax_MoisTunnel= Pmax_*0.95 #Pmax
Pmin_=-1.4 #Pmin
Pmin_MoisTunnel=Pmin_*0.95
delay_EMT_ms = 20 # 20 ms of delay for EMT simulations



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

Ratio_H_D_UP = 0.15 # Use to have two more values for D and H: D*(1+Ratio_H_D_UP), H*(1+Ratio_H_D_UP)
Ratio_H_D_Down = 0.15 # Use to have two more values for D and H: D*(1-Ratio_H_D_Down), H*(1-Ratio_H_D_Down)

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

#We get the value of P at Event_Time + DeltaSmallT
DeltaPAtEventTime = GetValueatSpecificTime(Event_Time+0.01,DeltaP)


results = [GetEnvelops(MargeUp,MargeDown,DeltaP_array[i],Tunnel_array[i],DeltaP) for i in range(len(D_array))]
results_PSecond_up_anal_array= [GetEnvelops(MargeUp,MargeDown,DeltaPSecond_up_anal_array[i],Tunnel_array[i],DeltaP) for i in range(len(D_array))]
results_PSecond_down_anal_array= [GetEnvelops(MargeUp,MargeDown,DeltaPSecond_down_anal_array[i],Tunnel_array[i],DeltaP) for i in range(len(D_array))]



#the reponse is based on DeltaP
P_up_anal_array, P_down_anal_array = map(np.array, zip(*results))
# The response is based on DeltaPSecond_up_anal_array, representing the upper exponential behavior.
# This applies in the case where DeltaP > 0.
# Both the upper and lower bounds of the response are calculated as follows:
# PSecond_up_up_anal_array   = (DeltaPSecond_up_anal_array * (1 + MarginUp) + Tunnel + P0)
# PSecond_up_down_anal_array = (DeltaPSecond_up_anal_array * (1 + MarginUp) + Tunnel + P0)
# A modification to PSecond_up_down_anal_array is done as clarified in "GetEnvelops"
PSecond_up_up_anal_array, PSecond_up_down_anal_array = map(np.array, zip(*results_PSecond_up_anal_array))
# The response is based on DeltaPSecond_down_anal_array, representing the lower exponential behavior.
# This applies in the case where DeltaP > 0.
# Both the upper and lower bounds of the response are calculated as follows:
# PSecond_down_up_anal_array   = (DeltaPSecond_down_anal_array * (1 + MarginUp) + Tunnel + P0)
# PSecond_down_down_anal_array = (DeltaPSecond_down_anal_array * (1 + MarginUp) + Tunnel + P0)
# A modification to PSecond_up_down_anal_array is done as clarified in "GetEnvelops"
PSecond_down_up_anal_array, PSecond_down_down_anal_array  = map(np.array, zip(*results_PSecond_down_anal_array))



#P_PCC is limited to Pmin and Pmax , is calculated here to be plot as the expected analytical behavior
P_PCC=Cutsignal(Pmin_,P0+DeltaP,Pmax_)
#Calculating 50% of P_PCC
P_50Prc = P0+ np.where(t_DeltaP >= Start_Time, DeltaP*0.5 , DeltaP)
P_50Prc=Cutsignal(Pmin_,P_50Prc,Pmax_)



print("Start_Time+0.2",Start_Time+0.2)




# Stack the arrays into a 2D NumPy array
# Selecting the maximum value among P_up_anal_array, [P_50Prc], PSecond_up_up_anal_array
# I don't believe we still need P_up_anal_array or P_50Prc,
# but I’ve left them in place to avoid re-testing the entire workflow at this stage.

stacked = np.vstack(( P_up_anal_array, [P_50Prc], PSecond_up_up_anal_array ))
# Compute the element-wise max, ignoring NaNs
P_up_finale = np.nanmax(stacked, axis=0)


# Stack the arrays into a 2D NumPy array
# Selecting the minimum value among P_down_anal_array, [P_50Prc], PSecond_down_down_anal_array
# I don't believe we still need P_down_anal_array or P_50Prc,
# but I’ve left them in place to avoid re-testing the entire workflow at this stage.

stacked = np.vstack((P_down_anal_array, [P_50Prc],PSecond_down_down_anal_array))
# Compute the element-wise max, ignoring NaNs
P_down_finale = np.nanmin(stacked, axis=0)


print("P_50Prc",P_50Prc)

#cutting Pdown finale and Pup finale in order to avoir reverse power
#P_up_finale ,  P_down_finale = LimitingReversePower(P_up_finale,P_down_finale, P0,Tunnel,DeltaP)

#Adding delays "delay_EMT_ms" when the simulation is done in EMT
if EMT:
    # Delay settings

    shift_Time = 0
    # Pad with the initial value instead of zero
    initial_value = P_up_finale[0]
    P_up_finale = delay_signal(delay_EMT_ms, fs, P_up_finale)
    P_down_finale = delay_signal(delay_EMT_ms, fs, P_down_finale)
    P_PCC = delay_signal(delay_EMT_ms, fs, P_PCC)
    P_50Prc = delay_signal(delay_EMT_ms, fs, P_50Prc)


    results = [delay_signal(delay_EMT_ms,fs, P_up_anal_array[i]) for i in range(len(D_array))]
    P_up_anal_array = results
    print("size P_up_anal_array:", len(P_up_anal_array))
    results = [delay_signal(delay_EMT_ms,fs, P_down_anal_array[i]) for i in range(len(D_array))]
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

#Plot EMT
BaseLocation= "EMTSimulations/Ppos_GFM_EMT_model_P0=0.5.csv"
csv_file_path_EMT = BaseLocation
#Name of the Columns
NameColumnsDataFrame = ["Time","Signal"]
# Read the CSV file into a DataFrame
dataUseCase_EMT = pd.read_csv(csv_file_path_EMT)
# Access the columns
time1 = dataUseCase_EMT['Time']
Ppos_GFM_EMT_model = dataUseCase_EMT['Signal']


#Plot EMT
BaseLocation= "EMTSimulations/Ppos_GFM_ideal_P0=0.5.csv"
csv_file_path_EMT = BaseLocation
#Name of the Columns
NameColumnsDataFrame = ["Time","Signal"]
# Read the CSV file into a DataFrame
dataUseCase_EMT = pd.read_csv(csv_file_path_EMT)
# Access the columns
time2 = dataUseCase_EMT['Time']
Ppos_GFM_ideal = dataUseCase_EMT['Signal']






shift_Time= 0

# Create the plot
plt.figure(figsize=(8, 5))  # Set figure size
plt.plot(t_shifted-1, y_selected, label="P_pcc from Open Modelica", color='red', linestyle='-')  # Adding simulation

plt.plot(t_DeltaP+shift_Time,P_PCC, label="Ppcc", color='m', linestyle='-')  # First plot
plt.plot(t_DeltaP+shift_Time,P_down_anal_array[0], label="Pdown_analytical", color='m', linestyle='--')  # First plot
plt.plot(t_DeltaP+shift_Time,P_up_anal_array[0], label="Pup_analytical", color='m', linestyle=':')  # First plot

plt.plot(t_DeltaP+shift_Time,P_down_anal_array[1], label="Pdown_analytical", color='black', linestyle='--')  # First plot
plt.plot(t_DeltaP+shift_Time,P_up_anal_array[1], label="Pup_analytical", color='black', linestyle=':')  # First plot

plt.plot(t_DeltaP+shift_Time,P_down_anal_array[2], label="Pdown_analytical", color='r', linestyle='--')  # First plot
plt.plot(t_DeltaP+shift_Time,P_up_anal_array[2], label="Pup_analytical", color='r', linestyle=':')  # First plot

#Plottijg the first order responses envelops of the second order response

plt.plot(t_DeltaP+shift_Time,PSecond_up_up_anal_array[0], label="PSecond_up_anal_array", color='m', linestyle=':')  # First plot
plt.plot(t_DeltaP+shift_Time,PSecond_down_down_anal_array[0], label="PSecond_down_anal_array", color='m', linestyle='-')  # First plot

plt.plot(t_DeltaP+shift_Time,PSecond_up_up_anal_array[1], label="PSecond_up_anal_array", color='black', linestyle=':')  # First plot
plt.plot(t_DeltaP+shift_Time,PSecond_down_down_anal_array[1], label="PSecond_down_anal_array", color='black', linestyle='-')  # First plot

plt.plot(t_DeltaP+shift_Time,PSecond_up_up_anal_array[2], label="PSecond_up_anal_array", color='r', linestyle=':')  # First plot
plt.plot(t_DeltaP+shift_Time,PSecond_down_down_anal_array[2], label="PSecond_down_anal_array", color='r', linestyle='-')  # First plot

plt.plot(t_DeltaP+shift_Time,P_50Prc, label="P_50%", linewidth='3')  # First plot

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
#plt.plot(t_shifted-1, y_selected, label="P_pcc from Open Modelica", color='b', linestyle='-')  # Adding simulation
plt.plot(t_DeltaP,P_PCC, label="P_PCC analytical", linewidth='3')  # First plot
plt.plot(t_DeltaP,P_down_finale, label="Pdown_final", linewidth=2)  # First plot
plt.plot(t_DeltaP,P_up_finale, label="Pup_final", linewidth=2)  # First plot
#plt.plot(time1, Ppos_GFM_EMT_model, label="Ppos_GFM_EMT_model from EMTP", color='gold', linestyle='-')  # EMT plot
#plt.plot(time2, Ppos_GFM_ideal, label="Ppos_GFM_ideal from EMTP", color='red', linestyle='-')  # EMT plot
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
plt.savefig(Path, bbox_inches='tight', dpi=300)

plt.show()

