import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

# Set the font properties
# Specify the path to your custom font file
font_path = './fonts/Roobert-Regular.ttf'  # Update this path to the actual font file location
# Create a FontProperties instance with your custom font
custom_font = FontProperties(fname=font_path, size=12)
custom_font_label = FontProperties(fname=font_path, size=10)

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
plt.xlabel('Sample', fontproperties=custom_font)
plt.ylabel('Count/%', fontproperties=custom_font)
plt.title('Mapped and Paired in sequencing reads', fontproperties=custom_font)
plt.legend(prop=custom_font_label)
for label in plt.gca().get_xticklabels():
    label.set_fontproperties(custom_font_label)
plt.savefig('../results/Beijing.png')
plt.show()

# Another Plot

# Calculate the ratio of "Mapped" to "Paired in sequencing" reads
data_bj_converted['Ratio1'] = data_bj_converted['Mapped'] / data_bj_converted['Paired in sequencing']

# Calculate the ratio of "Paired in sequencing" to "Total" reads
data_bj_converted['Ratio2'] = data_bj_converted['Paired in sequencing']

# Calculate the ratio of "Mapped" to "Total" reads
data_bj_converted['Ratio3'] = data_bj_converted['Mapped']

# Plot the ratio
plt.figure(figsize=(10, 6), dpi=300)
plt.plot(data_bj_converted['Ratio1'], label='Mapped/Paired in sequencing ratio')
plt.plot(data_bj_converted['Ratio2'], label='Paired in sequencing/Total ratio')
plt.plot(data_bj_converted['Ratio3'], label='Mapped/Total ratio')
plt.xlabel('Sample', fontproperties=custom_font)
plt.ylabel('Ratio/%', fontproperties=custom_font)
plt.title('Sequence read ratios', fontproperties=custom_font)
plt.legend(prop=custom_font_label)
for label in plt.gca().get_xticklabels():
    label.set_fontproperties(custom_font_label)
plt.savefig('../results/Beijing_ratio.png')
plt.show()
