import pandas as pd
import matplotlib.pyplot as plt
# Read the data
data_bj = pd.read_csv('D:/pythonCODE/envTest/BeijingSorted.txt', sep='\t')

plt.plot(data_bj['Mapped/Paired in sequencing'], label='Paired in sequencing/Mapped ratio')
plt.xlabel('Sample')
plt.ylabel('Ratio')
plt.title('	Mapped/Paired in sequencing')
plt.savefig('D:/pythonCODE/envTest/Beijing_pm2.png')
plt.show()