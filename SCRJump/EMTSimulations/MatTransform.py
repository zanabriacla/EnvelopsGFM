from scipy.io import loadmat
import pandas as pd
import matplotlib.pyplot as plt

# Load the .mat file
data = loadmat('/home/claudia/Documents/Scripts_GFM/DevPython/Git_Repo_GFM_Envelops/EnvelopsGFM/SCRJump/EMTSimulations/SCR_step_response.mat')

# View the keys (variables in the file)
print(data.keys())

# Access a variable
P_PCC = data['Ppos_GFM_ideal']
time =data['x1']
# Filter for time >= 3
start_time=3
mask = time >= start_time
time_filtered = time[mask]
P_PCC_filtered = P_PCC[mask]

# Shift time: subtract 3 from all values
time_shifted = time_filtered - start_time

# Save shifted data to CSV
df_shifted = pd.DataFrame({'Time': time_shifted, 'Signal': P_PCC_filtered})
df_shifted.to_csv('Ppos_GFM_ideal_P0=0.5.csv', index=False)

print(time)
plt.figure(figsize=(8, 5))  # Set figure size
plt.plot(time_shifted,P_PCC_filtered)
plt.show()
