#compare two polinolam like paths and plot them in the same plot
import csv
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import CubicSpline

#open a text file an get colums 2 and 3
with open('Data/Vid3.txt', 'r') as f:
    reader = csv.reader(f, delimiter=',')
    print(reader)
    x_csv = []
    y_csv = []
    # removed first line
    next(reader)
    for row in reader:
        print(row)
        x_csv.append(float(row[1]))
        y_csv.append(float(row[2]))
        
x_csv = np.asarray(x_csv)*1000
y_csv = np.asarray(y_csv)*1000        
#plot the data
xmin = 0
xmax = 1401
dx = 1
x = np.arange(xmin,xmax,dx)

#Skruehøyder:
#yfast = np.zeros(8)
#yfast[0] = 300
#yfast[1] = yfast[0] - np.random.randint(40,60)
#yfast[2] = yfast[1] - np.random.randint(70,90)
#yfast[3] = yfast[2] + np.random.randint(-30,10)
#yfast[4] = yfast[3] + np.random.randint(30,70)
#yfast[5] = yfast[4] + np.random.randint(-20,20)
#yfast[6] = yfast[5] - np.random.randint(40,80)
#yfast[7] = yfast[6] + np.random.randint(-40,40)
#
#print(yfast)
Data = [300, 253, 182, 175, 205, 189, 123, 135]
h = 200
xfast=np.asarray([0,1,2,3,4,5,6,7])*h
cs = CubicSpline(xfast,Data,bc_type='natural')
y = cs(x)
dy = cs(x,1)
d2y = cs(x,2)

yfast = Data.copy()
baneform = plt.figure('y(x)',figsize=(12,6))
plt.plot(x,y,xfast,yfast,'*', label='Ideelle baneform')
plt.plot(x_csv,y_csv,'y', label='Målt baneform')
plt.title('Banens form', fontsize=20)
plt.xlabel('$x$ (mm)',fontsize=20)
plt.ylabel('$y(x)$ (mm)',fontsize=20)
plt.text(10,80,'Skruehøyder (mm):', fontsize=16)
plt.text(-40, 50, int(yfast[0]), fontsize=16)
plt.text(160, 50, int(yfast[1]), fontsize=16)
plt.text(360, 50, int(yfast[2]), fontsize=16)
plt.text(560, 50, int(yfast[3]), fontsize=16)
plt.text(760, 50, int(yfast[4]), fontsize=16)
plt.text(960, 50, int(yfast[5]), fontsize=16)
plt.text(1160, 50, int(yfast[6]), fontsize=16)
plt.text(1360, 50, int(yfast[7]), fontsize=16)
plt.ylim(0,300)
plt.xlim(-50,1450)
plt.legend(fontsize=16)
plt.grid()
#plt.show()
plt.savefig('BaneformRe.png', bbox_inches='tight', dpi=600)