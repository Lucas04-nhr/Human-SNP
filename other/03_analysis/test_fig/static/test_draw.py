import matplotlib.pyplot as plt
import pandas as pd

# Load the data
data = pd.read_csv("BJ001.static.csv")

# Sort the data by count disc
data_sorted = data.sort_values(by='Count', ascending=False)
data_sorted['Flag'] = data_sorted['Flag'].astype(str)

# List the top 5 flags
flag_head = data_sorted['Flag'].head()
print("Top 5 flags:" + "\n" + str(flag_head))


# Define the figure
plt.figure(figsize=(10, 10), dpi=300)

# Draw the histogram
plt.bar(data_sorted['Flag'], data_sorted['Count'])
plt.xticks(rotation=90)
plt.xlabel("Flag")
plt.ylabel("Count")
plt.title("Flag counts in BJ001")
plt.tight_layout()
plt.savefig("BJ001.histogram.pdf")