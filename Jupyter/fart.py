import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import CubicSpline
import csv


#open a text file an get colums 2 and 3
with open('Data/Vid3.txt', 'r') as f:
    reader = csv.reader(f, delimiter=',')
    print(reader)
    x_csv = []
    y_csv = []
    # removed first line
    next(reader)
    next(reader)
    for row in reader:
        print(row)
        x_csv.append(float(row[1]))
        y_csv.append(float(row[3]))
        
x_csv = np.asarray(x_csv)*1000
y_csv = np.asarray(y_csv)*1000

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
#print(yfast)
Data = [300, 253, 182, 175, 205, 189, 123, 135]
h = 200
xfast=np.asarray([0,1,2,3,4,5,6,7])*h
cs = CubicSpline(xfast,Data,bc_type='natural')
y = cs(x)
dy = cs(x,1)
d2y = cs(x,2)

yfast = Data.copy()

#plt.savefig("Baneform", dpi = 600)
#Ta bort # hvis du ønsker å lagre grafen som pdf og/eller png.
#baneform.savefig("baneform.pdf", bbox_inches='tight')
#baneform.savefig("baneform.png", bbox_inches='tight')

y37 = y[400:1400]
y27 = y[200:1400]
y37min = np.min(y37)
y37max = np.max(y37)
y27min = np.min(y27)
y27max = np.max(y27)
K = d2y/(1+dy**2)**(1.5)
R = 1/(np.abs(K)+1E-8)  #unngår R = uendelig
Rmin = np.min(R)
beta = np.arctan(dy)
betadeg = beta*180/np.pi
startvinkel = betadeg[0]
maksvinkel = np.max(np.abs(betadeg))

print('Høyeste punkt etter 3.skrue (mm): %4.0f' %y37max)
print('Laveste punkt etter 2.skrue (mm): %4.0f' %y27min)
print('Starthelningsvinkel (grader): %4.1f' %startvinkel)
print('Maksimal helningsvinkel (grader): %4.1f' %maksvinkel)
print('Minste krumningsradius (mm): %4.0f' %Rmin)
print('Festepunkthøyder (mm):', yfast)

def finn_farten(y):
    c = 2/5
    y0 = 300
    g = 9810
    v = np.sqrt(2*g*(y0-y)/(1+c))
    
    return v

def plot_fart(y):
    v = finn_farten(y)
    fart = plt.figure('y(x)',figsize=(12,6))
    plt.plot(x, v, label='Ideell hastighet')
    plt.plot(x_csv, y_csv, 'y', label='Målt hastighet')
    plt.title('Fart vs posisjon', fontsize=20)
    plt.xlabel('$x$ (mm)',fontsize=20)
    plt.ylabel('$v$(mm/s)',fontsize=20)
    plt.grid()
    plt.legend()
    #plt.show()
    plt.savefig("fart.png", dpi = 600)

print(finn_farten(cs(1401))/1000)
plot_fart(y)