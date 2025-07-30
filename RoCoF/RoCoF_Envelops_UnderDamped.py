
import pandas as pd
import scipy.io
import matplotlib.pyplot as plt
import math
import numpy as np





def delay_signal(delay_ms,fs,signal):
    delay_samples = int((delay_ms / 1000) * 1 / fs)
    signal_delayed = np.concatenate((np.full(delay_samples, signal[0]), signal))[:len(signal)]
    return signal_delayed

def GetDeltaP(D_Damping,H,Xtotal_initial,P0):



    DeltaP = np.zeros_like(t_DeltaP)
    DeltaP, xi ,TresponseTime = GetDeltaP_NotDELAYED(D_Damping, H, Xtotal_initial, P0, t_DeltaP)  # <-- your actual function here



    DeltaP=np.where(Event_Time<t_DeltaP ,DeltaP,0)
    DeltaPSteadyState = DeltaP[-1] #Ppeak is the


    DeltaP_Recovered = GetDeltaP_NotDELAYED(D_Damping, H, Xtotal_initial, P0, t_DeltaP-RoCofStop_Time)[0]*-1
    DeltaP_Recovered = np.where(RoCofStop_Time < t_DeltaP, DeltaP_Recovered, 0)

    DeltaP= DeltaP + DeltaP_Recovered
    DeltaP=DeltaP*-1





    print("PsteadyState",DeltaPSteadyState)
    return DeltaP,DeltaPSteadyState,xi,TresponseTime

def GetDeltaP_NotDELAYED(D_Damping,H,x_total_initial,P0,t_DeltaP):

    epsilon = (D / 2 * np.sqrt(x_total_initial / (2 * H * wb * u_prod)) )
    wn = np.sqrt(wb * u_prod / (2 * H * x_total_initial))

    wd = wn * np.sqrt(1 - epsilon ** 2)


    alpha = 2 * H * t_pll * RoCoF
    beta = (2 * H + D * t_pll) / (2 * H * t_pll)

    A_coeff = alpha * beta
    common_denom = 1 - 2 * epsilon * wn * t_pll + wn ** 2 * t_pll ** 2
    B_coeff = -(t_pll ** 2 * alpha * wn ** 2 * (t_pll * beta - 1)) / common_denom
    C_coeff = (alpha * (2 * t_pll * beta * epsilon * wn - t_pll * wn ** 2 - beta)) / common_denom
    D_coeff = (
                      alpha
                      * (
                              4 * t_pll * beta * wn ** 2 * epsilon ** 2
                              - t_pll * beta * wn ** 2
                              - 2 * t_pll * wn ** 3 * epsilon
                              - 2 * beta * wn * epsilon
                              + wn ** 2
                      )
              ) / common_denom

    term2 = np.exp(-t_DeltaP / t_pll)
    term3 = np.exp(-epsilon * wn * t_DeltaP) * np.cos(wd * t_DeltaP)
    term4 = np.exp(-epsilon * wn * t_DeltaP) * np.sin(wd * t_DeltaP)

    delta_p1 = (
            A_coeff
            + B_coeff / t_pll * term2
            + C_coeff * term3
            + C_coeff * (D_coeff / C_coeff - epsilon * wn) / wd * term4
    )
    print("A_coeff",A_coeff)
    print("B_coeff", B_coeff)
    print("C_coeff", C_coeff)
    print("D_coeff", D_coeff)
    print("epsilon", epsilon)
    print("wd", wd)
    print("wn", wn)

    C2_coeff = (D_coeff-C_coeff*epsilon*wn)/wd
    R_coeff = np.sqrt(C_coeff ** 2 + ((D_coeff - C_coeff * epsilon * wn) / wd) ** 2)  # Amplitude factor

    DeltaPSecond_up = A_coeff + B_coeff/t_pll*np.exp(-t_DeltaP/t_pll)+R_coeff*np.exp(-epsilon*wn*t_DeltaP)
    DeltaPSecond_down = A_coeff +B_coeff/t_pll*np.exp(-t_DeltaP/t_pll)-R_coeff*np.exp(-epsilon*wn*t_DeltaP)

    Ppeak = A_coeff+B_coeff/t_pll+np.sqrt(C_coeff**2+(C_coeff*(D_coeff/C_coeff-epsilon*wn)/wd))
    TresponseTime = 4/(epsilon*wn)


    exp_decay_pll = (B_coeff / t_pll) * np.exp(-t_DeltaP / t_pll)
    R_coeff = np.sqrt(C_coeff ** 2 + ((D_coeff - C_coeff * epsilon * wn) / wd) ** 2)  # Amplitude factor
    upper_env = A_coeff + exp_decay_pll + R_coeff * np.exp(-epsilon * wn * t_DeltaP)
    lower_env = A_coeff + exp_decay_pll - R_coeff * np.exp(-epsilon * wn * t_DeltaP)

    return -delta_p1,Ppeak,epsilon,-upper_env,-lower_env,TresponseTime


#DeltaP fonction in case of SCR jump, we return also Ppeak, epsilon and two exponential up (DeltaP_up_anal_array) and low (DeltaP_down_anal_array)
#etDeltaP_NotDELAYED(D_Damping,H,x_total_initial,P0,time_array):
def GetDeltaP(D_Damping,H,x_total_initial,P0):
    delta_p1 = np.zeros_like(t_DeltaP)
    DeltaPSecond_up= np.zeros_like(t_DeltaP)
    DeltaPSecond_down = np.zeros_like(t_DeltaP)

    delta_p1,Ppeak,epsilon,delta_p1_up,delta_p1_down,TresponseTime = GetDeltaP_NotDELAYED(D_Damping, H, x_total_initial, P0, t_DeltaP)  # <-- your actual function here

    delta_p1 = np.where(Event_Time < t_DeltaP, delta_p1, 0)
    DeltaPSteadyState = delta_p1[-1]  # Ppeak is the
    delta_p1_up = np.where(Event_Time <t_DeltaP, delta_p1_up ,0 )
    delta_p1_down =np.where(Event_Time<t_DeltaP,delta_p1_down,0)


   # -delta_p2, Ppeak, epsilon, -delta_p2_up, -delta_p2_down, TresponseTime = GetDeltaP_NotDELAYED(D_Damping, H, x_total_initial,P0,t_DeltaP-RoCofStop_Time)
    delta_p2 =  GetDeltaP_NotDELAYED(D_Damping, H, x_total_initial,P0,t_DeltaP-RoCofStop_Time)[0]*-1
    delta_p2_up = GetDeltaP_NotDELAYED(D_Damping, H, x_total_initial, P0, t_DeltaP - RoCofStop_Time)[3] * -1
    delta_p2_down = GetDeltaP_NotDELAYED(D_Damping, H, x_total_initial, P0, t_DeltaP - RoCofStop_Time)[4] * -1

    delta_p2 = np.where(RoCofStop_Time < t_DeltaP, delta_p2, 0)
    delta_p2_up = np.where(RoCofStop_Time < t_DeltaP, delta_p2_up, 0)
    delta_p2_down = np.where(RoCofStop_Time < t_DeltaP, delta_p2_down, 0)

    DeltaP = delta_p1 + delta_p2


    DeltaP_up = delta_p1_up + delta_p2_up


    DeltaP_down = delta_p1_down + delta_p2_down


    TresponseTime= GetDeltaP_NotDELAYED(D_Damping, H, x_total_initial, P0, t_DeltaP)[5]

    print("TresponseTime", TresponseTime)
    return DeltaP, Ppeak, epsilon, DeltaP_up, DeltaP_down

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
    Signal_Down = np.where(Signal_Down < Pmin_, Pmin_ +0.2, Signal_Down)
    Signal_Down = np.where(Signal_Down > Pmax_, Pmax_ - 0.2, Signal_Down)
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

Z_init=1/SCR_init



RoCoF = -0.5/50  # Rate of Change of Frequency (Hz/s) ou pu ?
D=50#Damping constant of the VSM control
H=3 #Inertia constant (s)
t_pll = 0.01    # PLL time constant (s)
wb=314 # Base angular frequency(rad/s)
xtr=0.06 #Transformer reactance (pu)
Ugrid=1 # RMS voltage Ugrid (pu)
Uconv=1 # RMS voltage Uconverter (pu)
Xeff=0.25 # effective reactance (pu)
EMT= True # Can be "True" or "False" EMT is activated (20ms for the measures)
P0= 0 # Initial power (pu)
Pmax_=1.4 #Pmax
Pmax_MoisTunnel= Pmax_*0.95 #Pmax
Pmin_=-1.4 #Pmin
Pmin_MoisTunnel=Pmin_*0.95
delay_EMT_ms = 20 # 20 ms of delay for EMT simulations
u_prod = Uconv*Ugrid

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

# Define the time vector for simulation
Start_Time = -1 # sec
Event_Time=0 #keep this value to "0"
RoCofDuration=3 # duration of RoCof after that RoCof=0
RoCofStop_Time= Event_Time+RoCofDuration # sec
End_Time = 10 # sec
NbPoints = 10000
t_DeltaP = np.linspace(Start_Time, End_Time, NbPoints)  # From Start_Time to End_Time
fs = (End_Time - Start_Time) / NbPoints  # Sampling frequency (Hz)


# Calculating VARIABLES that need to be defined to calculate DeltaP
Xtotal_initial = Xeff + Z_init  # X total initial that is equal to Xeff+Xgrid final


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
DeltPUperExpBound = DeltaPSecond_up_anal_array[0]
DeltPLowerExpBound = DeltaPSecond_down_anal_array[0]

#We get the value of P at Event_Time + DeltaSmallT
DeltaPAtEventTime = GetValueatSpecificTime(Event_Time+0.01,DeltaP)


results = [GetEnvelops(MargeUp,MargeDown,DeltaP_array[i],Tunnel_array[i],DeltaP) for i in range(len(D_array))]
#results_PSecond_up_anal_array= [GetEnvelops(MargeUp,MargeDown,DeltaPSecond_up_anal_array[i],Tunnel_array[i],DeltaP) for i in range(len(D_array))]
#results_PSecond_down_anal_array= [GetEnvelops(MargeUp,MargeDown,DeltaPSecond_down_anal_array[i],Tunnel_array[i],DeltaP) for i in range(len(D_array))]



#the reponse is based on DeltaP
P_up_anal_array, P_down_anal_array = map(np.array, zip(*results))
# The response is based on DeltaPSecond_up_anal_array, representing the upper exponential behavior.
# This applies in the case where DeltaP > 0.
# Both the upper and lower bounds of the response are calculated as follows:
# PSecond_up_up_anal_array   = (DeltaPSecond_up_anal_array * (1 + MarginUp) + Tunnel + P0)
# PSecond_up_down_anal_array = (DeltaPSecond_up_anal_array * (1 + MarginUp) + Tunnel + P0)
# A modification to PSecond_up_down_anal_array is done as clarified in "GetEnvelops"
#PSecond_up_up_anal_array, PSecond_up_down_anal_array = map(np.array, zip(*results_PSecond_up_anal_array))
# The response is based on DeltaPSecond_down_anal_array, representing the lower exponential behavior.
# This applies in the case where DeltaP > 0.
# Both the upper and lower bounds of the response are calculated as follows:
# PSecond_down_up_anal_array   = (DeltaPSecond_down_anal_array * (1 + MarginUp) + Tunnel + P0)
# PSecond_down_down_anal_array = (DeltaPSecond_down_anal_array * (1 + MarginUp) + Tunnel + P0)
# A modification to PSecond_up_down_anal_array is done as clarified in "GetEnvelops"
#PSecond_down_up_anal_array, PSecond_down_down_anal_array  = map(np.array, zip(*results_PSecond_down_anal_array))



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

#stacked = np.vstack(( P_up_anal_array, [P_50Prc], PSecond_up_up_anal_array ))
# Compute the element-wise max, ignoring NaNs
#P_up_finale = np.nanmax(stacked, axis=0)


# Stack the arrays into a 2D NumPy array
# Selecting the minimum value among P_down_anal_array, [P_50Prc], PSecond_down_down_anal_array
# I don't believe we still need P_down_anal_array or P_50Prc,
# but I’ve left them in place to avoid re-testing the entire workflow at this stage.

#stacked = np.vstack((P_down_anal_array, [P_50Prc],PSecond_down_down_anal_array))
# Compute the element-wise max, ignoring NaNs
#P_down_finale = np.nanmin(stacked, axis=0)


print("P_50Prc",P_50Prc)

#cutting Pdown finale and Pup finale in order to avoir reverse power
#P_up_finale ,  P_down_finale = LimitingReversePower(P_up_finale,P_down_finale, P0,Tunnel,DeltaP)

#Adding delays "delay_EMT_ms" when the simulation is done in EMT
if EMT:
    # Delay settings

    shift_Time = 0
    # Pad with the initial value instead of zero
#    initial_value = P_up_finale[0]
 #   P_up_finale = delay_signal(delay_EMT_ms, fs, P_up_finale)
  #  P_down_finale = delay_signal(delay_EMT_ms, fs, P_down_finale)
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





# Create the plot



# Create the plot
Title= "P0="+ str(P0) +", RoCoF=" + str(RoCoF) +"pu"+ ", Duration="+str(RoCofDuration) +"s" +  ", Epsilon= " + str(round(epsilon,3)) + ", ωd= " + ", D= " + str(D) + ", H= " +str(H) + ", Xeff= " + str(Xeff)+ ", EMT=" +str(EMT) +", SCR="+str(SCR_init)



shift_Time= 0

# Create the plot
plt.figure(figsize=(8, 5))  # Set figure size

plt.plot(t_DeltaP+shift_Time,DeltPUperExpBound, label="UperBound", color='r', linestyle=':')  # First plot
plt.plot(t_DeltaP+shift_Time,DeltPLowerExpBound, label="LowerBound", color='r', linestyle='-')  # First plot
plt.plot(t_DeltaP+shift_Time,DeltaP, label="DeltaP", color='r', linestyle='-')  # First plot


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
#P_up_finale,P_down_finale,P_PCC = DelayEnvelops(P_up_finale,P_down_finale,P_PCC,shift_Time)




# Create the plot
plt.figure(figsize=(8, 5))  # Set figure size
#plt.plot(t_shifted-1, y_selected, label="P_pcc from Open Modelica", color='b', linestyle='-')  # Adding simulation
plt.plot(t_DeltaP,P_PCC, label="P_PCC analytical", linewidth='3')  # First plot

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

