import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


# Read the data
data_bj = pd.read_csv('../original/Beijing.txt', sep='\t')

# Replace ".sorted" with an empty string in the first column
data_bj.iloc[:, 0] = data_bj.iloc[:, 0].str.replace(".sorted", "", regex=False)

# Display the modified DataFrame to verify changes
print(data_bj.head())

# Select only the "Sample", "Mapped", and "Paired in sequencing" columns
data_bj_filtered = data_bj[['Sample', 'Total', 'Mapped', 'Paired in sequencing']]

# Display the filtered DataFrame to verify the selected columns
print(data_bj_filtered.head())

# Prepare the data for plotting
data_bj_converted = data_bj_filtered.astype({'Mapped': int, 'Paired in sequencing': int})
data_bj_converted['Mapped'] = data_bj_converted['Mapped'] / data_bj_converted['Total'] * 100
data_bj_converted['Paired in sequencing'] = data_bj_converted['Paired in sequencing'] / data_bj_converted['Total'] * 100

# Plot the data
plt.figure(figsize=(10, 6), dpi=300)
plt.plot(data_bj_converted['Mapped'], label='Mapped')
plt.plot(data_bj_converted['Paired in sequencing'], label='Paired in sequencing')
plt.xlabel('Sample')
plt.ylabel('Count/%')
plt.title('Mapped and Paired in sequencing reads')
plt.legend()
plt.savefig('../results/Beijing.png')
plt.show()

# Another Plot

# Calculate the ratio of "Mapped" to "Paired in sequencing" reads
data_bj_converted['Ratio1'] = data_bj_converted['Mapped'] / data_bj_converted['Paired in sequencing']

# Calculate the ratio of "Paired in sequencing" to "Total" reads
data_bj_converted['Ratio2'] = data_bj_converted['Paired in sequencing'] / 100

# Calculate the ratio of "Mapped" to "Total" reads
data_bj_converted['Ratio3'] = data_bj_converted['Mapped'] / 100

# Plot the ratio
plt.figure(figsize=(10, 6), dpi=300)
plt.plot(data_bj_converted['Ratio1'], label='Mapped/Paired in sequencing ratio')
plt.plot(data_bj_converted['Ratio2'], label='Paired in sequencing/Total ratio')
plt.plot(data_bj_converted['Ratio3'], label='Mapped/Total ratio')
plt.xlabel('Sample')
plt.ylabel('Ratio/%')
plt.title('Sequence read ratios')
plt.legend()
plt.savefig('../results/Beijing_ratio.png')
plt.show()
