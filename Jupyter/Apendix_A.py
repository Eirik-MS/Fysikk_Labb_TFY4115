import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import CubicSpline

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
baneform = plt.figure('y(x)',figsize=(12,6))
plt.plot(x,y,xfast,yfast,'*')
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
plt.grid()
plt.show()
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
    plt.plot(x, v)
    plt.title('Fart vs posisjon', fontsize=20)
    plt.xlabel('$x$ (mm)',fontsize=20)
    plt.ylabel('$v$(mm/s)',fontsize=20)
    plt.grid()
    plt.savefig("fart.png", dpi = 600)
    
plot_fart(y)

#kinetisk energi er gitt som 1/2*mv**2
def plot_kinetic(y):
    #Total kinetisk energi K er summen av transelasjonsenergi mv**2/2 og rotasjonsenergi c*m*v**2/2
    #E = U = m*g*y_0
    #Konstanter:
    m = 0.031 #kg
    r = 0.011 #m
    c = 2/5
    Kinetic = ((1+c)/2)*m*finn_farten(y)**2
    Kinetic = Kinetic / 1000 # for å få energi i J
    kin = plt.figure('y(x)',figsize=(12,6))
    plt.plot(x, Kinetic)
    plt.title('Kinetic energi', fontsize=20)
    plt.xlabel('$x$ (mm)',fontsize=20)
    plt.ylabel('$E$(J)',fontsize=20)
    plt.grid()
    plt.show()


plot_kinetic(y)

def krumning_kappa():
    kappa = d2y/((1+(dy)**2)**(3/2))

    Krumning = plt.figure('y(x)',figsize=(12,6))


    plt.plot(x,kappa)
    plt.title('Banens krumning', fontsize=20)
    plt.xlabel('$x$ (mm)',fontsize=20)
    plt.ylabel('$K$(1/mm)',fontsize=20)
    plt.grid()
    plt.show()
    return kappa
    
krumning_kappa()

def plot_helningsvinkel():

    beta = np.arctan(dy)

    vinkel = plt.figure('y(x)',figsize=(12,6))
    plt.plot(x,beta)
    plt.title('Banens helningsvinkel', fontsize=20)
    plt.xlabel('$x$ (mm)',fontsize=20)
    plt.ylabel('$ß$(grader)',fontsize=20)
    plt.grid()
    plt.show()

plot_helningsvinkel()


def finn_tid():
    delta_v_x = np.zeros(len(y))
    delta_t = np.zeros(len(y))
    V = finn_farten(y)
    t = np.zeros(len(y))
    for i in range(1,1401):
        delta_v_x[i] = 1/2* (V[i]-V[i-1])* np.cos(beta[i])
        delta_t[i] = 2* delta_v_x[i]/(V[i]+V[i-1])* np.cos(beta[i])
        t[i] = t[i-1]+delta_t[i]
    return t[i]

finn_tid()
    
def sentripetalakselerasjon():
    a_ = (finn_farten(y)**2)*krumning_kappa()
    return a_



def normalKraft():
    g = 9810    
    m = 0.031
    N = m*(g*np.cos(beta)+sentripetalakselerasjon())
    return N
  

def plot_normalkraft():
    kraft = plt.figure('y(x)',figsize=(12,6))
    plt.plot(x,normalKraft())
    plt.title('Normalkraft', fontsize=20)
    plt.xlabel('$x$ (mm)',fontsize=20)
    plt.ylabel('$N$(mg)',fontsize=20)
    plt.grid()
    plt.show()

plot_normalkraft()
  
def friksjon():
    #Konstanter
    c = 2/5
    m = 0.031 #kg
    g = 9810 #mm/s**2
    return (c/(1+c)*m*g*np.sin(beta))/1000

def plot_friksjon():
    kraft = plt.figure('y(x)',figsize=(12,6))
    plt.plot(x,friksjon())
    plt.title('Friksjonskraft', fontsize=20)
    plt.xlabel('$x$ (mm)',fontsize=20)
    plt.ylabel('$f$(mg)',fontsize=20)
    plt.grid()
    plt.show()

plot_friksjon()


def plot_NdF():
    F = friksjon()
    N = normalKraft()
    NF = np.abs(F/N)
    kraft = plt.figure('y(x)',figsize=(12,6))
    plt.plot(x,NF)
    plt.title('Forholdet f/N', fontsize=20)
    plt.xlabel('$x$ (mm)',fontsize=20)
    plt.ylabel('$|f/N$|',fontsize=20)
    plt.grid()
    plt.show()

plot_NdF()

h_start = 0.287 #[m]

h_slutt = 0.155 #[m]

m_kule = 0.037 #[kg]

g = 9.81 #[m/s**2]

V_slutt = [1.391,1.377,1.35828 ,1.36367, 1.407, 1.361, 1.343, 1.345, 1.3435, 1.343, 1.339] #[m/s]

rulletid = [1.55, 1.54, 1.500, 1.5333, 1.500, 1.533, 1.500, 1.567, 1.533, 1.433] #[s]

def gjennomsnitt_rulletid():
    rulletid_sum = sum(rulletid)/len(rulletid)        
    return rulletid_sum


def gjennomsnitt_sluttfart():
    V_slutt_sum = sum(V_slutt)/len(V_slutt)        
    return V_slutt_sum


def tap_mekanisk_energi(): 
    Energi_sum = 0
    #finner gjennomsnittstap på de ti forsøkene
    
    for i in V_slutt:
        Energi_sum = h_start*g*m_kule - (i**2 * 0.5 * m_kule + m_kule*g*h_slutt)


        
    
    Energitap_gjennomsnitt = Energi_sum/10
    
    print(Energitap_gjennomsnitt)
    
    return Energitap_gjennomsnitt 

def tap_mekanisk_energi_liste():
    liste =[]
    for i in V_slutt:
        Energi_sum = h_start*g*m_kule - (i**2 * 0.5 * m_kule + m_kule*g*h_slutt)
        liste.append(Energi_sum)

    return liste 

def total_kinetisk():
    kin_sum = 0
    for i in V_slutt:
        m = 0.031
        r = 0.011
        c = 2/5
        Kinetic = ((1+c)/2)*m*i**2
        Kinetic = Kinetic / 1000
        kin_sum =+ Kinetic

    kin_gjennomsnitt = kin_sum/10

    print(kin_gjennomsnitt)
    return kin_gjennomsnitt

def finn_standaravik(data):
    sum1 = 0
    for value in data:
        sum1 = (value - tap_mekanisk_energi())**2
    delta_X = np.sqrt(1/(len(data)-1) * sum1)
    return delta_X

def standarfeil(data):
    standarfeil = finn_standaravik(data)/np.sqrt(len(data))
    return standarfeil

standarfeil(rulletid)
standarfeil(total_kinetisk())
standarfeil(V_slutt)

def view_diff():
    diff_list = []
    simulert_hatighet = yfast[-1]
    eksprementiell_hatighet = V_slutt
    for i in eksprementiell_hatighet:
        sum += simulert_hatighet - i
    diff_hastighet = sum / len(eksprementiell_hatighet)
    
    diff_list.append(diff_hastighet)
    
    simulert_rulletid = finn_tid()
    eksprementiell_rulletid = gjennomsnitt_rulletid(rulletid)
    
    diff_list.append(simulert_rulletid-eksprementiell_rulletid)
    
    
    return diff_list

print(f"Forskjell i hastighet: {view_diff()[0]} \nForskjell i rulletid {view_diff()[1]}")