
import pandas as pd
import scipy.io
import matplotlib.pyplot as plt
import math
import numpy as np
from sympy.abc import s, t
from sympy import symbols, Function, laplace_transform, inverse_laplace_transform, exp, simplify
from scipy.integrate import quad
from sympy import symbols, Function, Heaviside

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

def GetDeltaP_NotDELAYED(D_Damping,H,Xtotal_initial,P0,t_DeltaP):
    print("time",t_DeltaP)
    omega_n = np.sqrt(wb * Uconv * Ugrid / (2 * H * Xtotal))  # Natural frequency (rad/s)
    xi = omega_n * D_damping * Xtotal / (2 * wb * Uconv * Ugrid)

    if xi > 1:
        xi_over = xi  # Overdamped damping ratio (ξ > 1)
        print("xi_over", xi_over)
    else:
        xi_under = xi  # Underdamped damping ratio (ξ < 1)
        print("xi_under", xi_under)

    print("xi", xi)

    print("natural frq", omega_n)




    t= t_DeltaP
    print("time_delay", t)
    print("T_Delay",t)
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
    Ppeak=DeltaP[0]
    return DeltaP,Ppeak,xi


def GetDeltaP(D_Damping,H,Xtotal_initial,P0):




    # Define time variable and shift amount
    #t, T = symbols('t T', real=True, positive=True)


    DeltaP = np.zeros_like(t_DeltaP)
    f_t,Ppeak, xi = GetDeltaP_NotDELAYED(D_Damping, H, Xtotal_initial, P0, t_DeltaP)  # <-- your actual function here

    print(f_t)
    f_t=np.where(t1<t_DeltaP ,f_t,0)
    Ppeak = f_t[-1]
    f_t_delayed = GetDeltaP_NotDELAYED(D_Damping, H, Xtotal_initial, P0, t_DeltaP-t2)[0]*-1
    f_t_delayed = np.where(t2 < t_DeltaP, f_t_delayed, 0)
    DeltaP= f_t + f_t_delayed
    DeltaP=DeltaP*-1

    print("Ppeak",Ppeak)
    return DeltaP,Ppeak,xi

def GetTunnel(Ppeak):
    Tunnel = max(0.02, 0.05 * Ppeak)
    return Tunnel

def Cutsignal(Valuemin,Signal,Valuemax):
    Signal = np.where(Signal < Valuemin, Valuemin, Signal)
    Signal = np.where(Signal > Valuemax, Valuemax, Signal)
    print("Value Min:", Valuemin)
    print("Value Max:======", Valuemax,"-")
    return Signal




def GetEnvelops(MargeUp,MargeDown,Signal,Tunnel):
    tx=0.5

    print('DeltaP ',Signal[2000])
    # Creating Envelops
    if RoCoF <= 0:
        print('I am here ')

        margin=AddMargin(MargeUp, Tunnel, RoCofStop_Time)

        # Envelope curves
        upper_envelope = Signal + margin+P0
        # Find index of closest time
        index = np.argmin(np.abs(t_DeltaP - (RoCofStop_Time-0.01))) # taking the value of P 10ms before RoCofStop_Time

        # Get value from the signal
        value_at_RoCofStop_Time = upper_envelope[index]
        print("value_at_RoCofStop_Time", str(value_at_RoCofStop_Time))
        mask = (t_DeltaP >= RoCofStop_Time) & (t_DeltaP <= End_Time)
        condition = mask & (upper_envelope > (value_at_RoCofStop_Time))
        upper_envelope = np.where(condition, value_at_RoCofStop_Time, upper_envelope)

        #################################################################################


        margin = AddMargin(MargeDown, Tunnel, RoCofStop_Time)

        lower_envelope = Signal - margin+P0
        # Plot
        mask = (t_DeltaP >= Start_Time) & (t_DeltaP <= RoCofStop_Time)
        condition = mask & (lower_envelope < (P0 - Tunnel))
        lower_envelope = np.where(condition, P0 - Tunnel, lower_envelope)


        #Signal_up_anal = Signal * (1 + MargeUp) + Tunnel + P0
        Signal_up_anal = upper_envelope

        #Signal_down_anal = Signal * (1 - MargeDown) - Tunnel + P0
        Signal_down_anal = lower_envelope


    else:
###############################################################
       # Signal_up_anal = Signal * (1 - MargeUp) + Tunnel + P0
################################################################""
        margin=AddMargin(MargeUp, Tunnel, RoCofStop_Time)
        # Envelope curves
        upper_envelope = Signal + margin+P0
        mask = (t_DeltaP >= Start_Time) & (t_DeltaP <= RoCofStop_Time)

        condition = mask & (upper_envelope > (P0 + Tunnel))
        print(upper_envelope)
        print(np.any(condition))  # Should be True if any value exceeds P0 + Tunnel during the window
        print(np.max(upper_envelope[mask]))  # Check actual values
        print(P0 + Tunnel)  # Compare to see if condition could ever be True

        upper_envelope = np.where(mask & (upper_envelope > (P0 + Tunnel)), P0 + Tunnel, upper_envelope)


#####################################################""
       # Signal_down_anal = Signal * (1 + MargeDown) - Tunnel + P0
########################################################

        margin=AddMargin(MargeDown, Tunnel, RoCofStop_Time)
        #Envelope curves
        lower_envelope = Signal - margin + P0


        # Find index of closest time
        index = np.argmin(np.abs(t_DeltaP - (RoCofStop_Time-0.01))) # taking the value of P 10ms before RoCofStop_Time
        # Get value from the signal
        value_at_RoCofStop_Time = lower_envelope[index]
        print("value_at_RoCofStop_Time", str(value_at_RoCofStop_Time))
        mask = (t_DeltaP >= RoCofStop_Time) & (t_DeltaP <= End_Time)
        condition = mask & (lower_envelope < (value_at_RoCofStop_Time))
        lower_envelope = np.where(condition, value_at_RoCofStop_Time, lower_envelope)

       # lower_envelope = np.where(lower_envelope > (P0 - Tunnel), lower_envelope, P0 - Tunnel)

        # Signal_up_anal = Signal * (1 + MargeUp) + Tunnel + P0
        Signal_up_anal = upper_envelope

        # Signal_down_anal = Signal * (1 - MargeDown) - Tunnel + P0
        Signal_down_anal = lower_envelope

    Signal_up_anal = Cutsignal(Pmin_,Signal_up_anal,Pmax_)
    Signal_down_anal = Cutsignal(Pmin_, Signal_down_anal, Pmax_)
    print(type(Signal_up_anal),"the type of up_anal")
    return Signal_up_anal,Signal_down_anal

def AddMargin(MargeUp,Tunnel,RoCofStop_Time):
     initial_margin = MargeUp
     final_margin = Tunnel
     decay_rate = 3  # tune this to control how fast the margin narrows
     RoCofStop_Time = RoCofStop_Time
     #margin = final_margin + (initial_margin - final_margin) * np.exp(-decay_rate * t_DeltaP) + (0.1) * np.exp(-decay_rate * (t_DeltaP - 2))
     margin = initial_margin * ((RoCofStop_Time>=t_DeltaP) &(t_DeltaP>= Event_Time)) * np.exp(-decay_rate * t_DeltaP) + initial_margin * (t_DeltaP >= RoCofStop_Time) * np.exp(-decay_rate * (t_DeltaP - RoCofStop_Time)) + Tunnel
     return margin

#Variables needed to be fulfilled in order to implement the envelops



RoCoF = 0.01  # Rate of Change of Frequency (Hz/s)
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

Z_grid=1/SCR


'''
if EventOnZgrid == "up":
    Delta_ZGrid = EventOnZgrid_up(Z1,Z2)[0]
elif EventOnZgrid == "down":
    Delta_ZGrid = EventOnZgrid_down(Z1, Z2)[0]
else:
    print("Invalid mode")

'''
print("Final DeltaP",RoCoF*(2*H+D_damping*T_pll))


#second Order system

# Define the time vector for simulation
Start_Time = -1
End_Time = 2
RoCofStop_Time= 3
Event_Time=0

t0 = Start_Time
t1 = Event_Time
t2 = RoCofStop_Time
t3 = End_Time

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
Ppeak_array = []
Tunnel_array = []
Epsilon_array = []
P_up_anal_array = []
P_down_anal_array = []

D_array=[D_damping,D_damping*(1+Ratio_H_D_UP),D_damping*(1-Ratio_H_D_Down)]
H_array=[H,H*(1-Ratio_H_D_UP),H*(1+Ratio_H_D_Down)]
print("Set of D values to be considered:",D_array)
print("Set of H values to be considered", H_array)

#Retrieving the second order response and the Tunnel that will be used in the Margins

results = [GetDeltaP(D_array[i], H_array[i], Xtotal, P0) for i in range(len(D_array))]
DeltaP_array, Ppeak_array , Epsilon_array= map(np.array, zip(*results))
Tunnel_array = [GetTunnel(Ppeak_array[i]) for i in range(len(D_array))]


#Creating Envelops
MargeUp=0.2 # This is the Margin up used in DeltaP*(1+-MarginUp)+Tunnel
MargeDown=0.2 # This is the Margin down used in DeltaP*(1+-MargeDown)+Tunnel
DeltaP = DeltaP_array[0]
Tunnel = Tunnel_array[0]
epsilon = Epsilon_array[0]

print(GetEnvelops(MargeUp,MargeDown,DeltaP_array[0],Tunnel_array[0]))

results = [GetEnvelops(MargeUp,MargeDown,DeltaP_array[i],Tunnel_array[i]) for i in range(len(D_array))]
P_up_anal_array, P_down_anal_array = map(np.array, zip(*results))

#Theoretical Value
#P_PCC= P0+DeltaP
P_PCC=Cutsignal(Pmin_,P0+DeltaP,Pmax_)
#P_PCC = np.where(P_PCC < -1, -1, P_PCC)



#Envelop of 50% of Delta before t=10ms and after that it takes DeltaP
#P_50Prc = P0+ np.where(t_DeltaP < 0.1, DeltaP*0.5 , DeltaP)
#P_50Prc = P_50Prc+0.02
#P_50Prc=Cutsignal(Pmin_,P_50Prc,Pmax_)

#From different possibilities , we select the max value for Envelop UP and the MIN value for envelop DOWN
#P_up_finale = np.maximum(P_up_anal_array[0], P_up_anal_array[1],P_50Prc)

# Element-wise max across all


#P_up_finale = np.maximum.reduce(P_up_anal_array + [P_50Prc])
#P_down_finale = np.minimum.reduce(P_down_anal_array + [P_50Prc])


# Stack the arrays into a 2D NumPy array
stacked = np.vstack((P_up_anal_array ))
# Compute the element-wise max, ignoring NaNs
P_up_finale = np.nanmax(stacked, axis=0)

# Stack the arrays into a 2D NumPy array
stacked = np.vstack((P_down_anal_array ))
# Compute the element-wise max, ignoring NaNs
P_down_finale = np.nanmin(stacked, axis=0)

LocationFile= "P0="+ str(P0) +", RoCoF" + str(RoCoF) +  ", Epsilon= " + str(round(epsilon,3)) + ", ωd= " + ", D= " + str(D_damping) + ", H= " +str(H) + ", Xeff= " + str(Xeff)+".csv"
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

# Create the plot
plt.figure(figsize=(8, 5))  # Set figure size
plt.plot(t_DeltaP+shift_Time,P_PCC, label="Ppcc", color='black', linestyle='--')  # First plot
plt.plot(t_DeltaP+shift_Time,P_down_anal_array[0], label="Pdown_analytical", color='m', linestyle='--')  # First plot

plt.plot(t_DeltaP+shift_Time,P_up_anal_array[0], label="Pup_analytical", color='r', linestyle='--')  # First plot
plt.plot(t_DeltaP+shift_Time,P_down_anal_array[1], label="Pdown_analytical", color='m', linestyle=':')  # First plot
plt.plot(t_DeltaP+shift_Time,P_up_anal_array[1], label="Pup_analytical", color='r', linestyle=':')  # First plot
plt.plot(t_DeltaP+shift_Time,P_down_anal_array[2], label="Pdown_analytical", color='m', linestyle='-')  # First plot
plt.plot(t_DeltaP+shift_Time,P_up_anal_array[2], label="Pup_analytical", color='r', linestyle='-')  # First plot

#plt.plot(t_DeltaP+shift_Time,P_50Prc, label="P_50%", linewidth='3')  # First plot
plt.plot(t_DeltaP+shift_Time,P_down_finale, label="Pdown_final", linewidth=2)  # First plot
plt.plot(t_DeltaP+shift_Time,P_up_finale, label="Pup_final", linewidth=2)  # First plot


# Add vertical line at t = 10 seconds
plt.axvline(x=0.010, color='black', linestyle='--', label='t = 10ms')

# Add labels, title, and legend
plt.xlabel("sec")
plt.ylabel("P at PCC (pu)")
plt.title("RoCof= " +str(RoCoF) + "pu, SCR="+str(SCR) + ", Epsilon= " + str(round(epsilon,3)) + ", ωd= " + ", D= " + str(D_damping) + ", H= " +str(H) + "sec., Xeff= " + str(Xeff)+"pu")
plt.legend()  # Show legend

# Show the plot
plt.grid(True)  # Add grid for better visualization
plt.show()

