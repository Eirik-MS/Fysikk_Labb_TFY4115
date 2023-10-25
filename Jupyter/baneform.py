#compare two polinolam like paths and plot them in the same plot
import csv
#open a text file an get colums 2 and 3
with open('Data/Vid3.txt', 'r') as f:
    reader = csv.reader(f, delimiter=',')
    print(reader)
    x = []
    y = []
    # removed first line
    next(reader)
    for row in reader:
        print(row)
        x.append(float(row[1]))
        y.append(float(row[2]))
        
#plot the data

import matplotlib.pyplot as plt

plt.plot(x,y)
plt.title('y(x)', fontsize=20)
plt.xlabel('$x$ (mm)',fontsize=20)
plt.ylabel('$y$(mm)',fontsize=20)
plt.grid()
plt.show()
