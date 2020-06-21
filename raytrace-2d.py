import numpy as np
import matplotlib.pyplot as plt

#die linspaces sind sowas, wie die Intervalle, auf denen gewisse Geraden bze. Formen definiert sind

pi = np.pi
R1 = 100.0 #Radius der ersten Kugel, für die normale Verwendung wird R1 = 800 empfohlen, fuer sphärische Aberration R1=100 
R2 = R1 #Radius der zweiten Kugel
num = 500 #Anzahl der Punkte, die fuer jeden linspace erstellt werden
n0 = 1.0 #Brechungsindices
nD = 1.5
h = 70.7 # gleiche Hoehe, wie in der 3d Simulation
d1 = (R1-(R1**2-h**2)**0.5) # Jeweilige Durchmesser Haelften um gegebene Hoehe zu erreichen
d2 = (R2-(R2**2-h**2)**0.5)
d = d1+d2 #Gesamt Durchmesser der Linse

def Linsenschleifergleichung():
    D = (nD-1.0)*(1.0/R1-1.0/(-R2)+d*(nD-1.0)/(nD*R1*(-R2)))
    return 1.0/D

brennweite = Linsenschleifergleichung() # praktische Brennweite mit komischem Korrekturfaktor
print(Linsenschleifergleichung())


q1 = 2.0*brennweite+R1-d1 # Ursprünge der Kreise
q2 = 2.0*brennweite-R2+d2
weitung = 1.1 # Die linspaces werden um diesen Faktor verlaengert
dist = 4.0*brennweite # Entfernung von Objekt zu Schirm


def norm(x): # Norm/Betrag eines Vektors
    v = 0
    for i  in range(0, len(x)):
        v += x[i]**2
    return v**0.5

def sp(x, y): # Standardskalarprodukt zweier Vektoren
    v = 0
    for i  in range(0, len(x)):
        v += x[i]*y[i]
    return v

def vp(x, y): # Vektorprodukt
    v1 = x[1]*y[2]-x[2]*y[1]
    v2 = x[2]*y[0]-x[0]*y[2]
    v3 = x[0]*y[1]-x[1]*y[0]
    return np.array([v1, v2, v3])  

def angle2(x, y): # Winkel zwischen zwei Vektoren
    return np.arccos(np.abs(sp(x, y))/(norm(x)*norm(y)))

def angle(x, y):
    x2 = np.array([x[0], x[1], 0.0])
    y2 = np.array([y[0], y[1], 0.0])
    return np.arcsin(vp(x2, y2)[2]/(norm(x)*norm(y)))

def abc(a, b, c, index): #pq bzw. abc Formel. Der index gibt an welche der beiden Loesungen gewollt ist
    p = b/a
    q = c/a
    return -0.5*p + index * np.sqrt(0.25*p**2 - q)


def SchnittpunktKreis(a, b, q, r, index): # Bestimmt den Schnittpunkt zwischen einer Gerade a+k*b und einem Kreis (mit Mittelpunkt q) des Radiuses r um den Punkt z 
     C1 = sp(b, b)
     C2 = 2.0*b[0]*(a[0]-q) + 2.0*b[1]*a[1]
     C3 = (a[0]-q)**2 + a[1]**2 - r**2
     k = abc(C1, C2, C3, index) # Faktor k in der Geradengleichung a+k*b
     S = a+k*b # Schnittpunkt a+k*b
     L = S-np.array([q, 0.0]) # Lot, bzw Vektor von Zentrum zu Schnittpunkt
     L_ = L/norm(L)
     return np.array([S, L_])

def SchnittpunktGerade(a, b, c, d): # Gerade a+k*b trifft Gerade c+r*d
    Ergebnis = a-c
    Matrix = np.array(
        [[d[0], -b[0]],
        [d[1], -b[1]]])
    Loesung = np.linalg.solve(Matrix, Ergebnis)
    Schnittpunkt = a+b*Loesung[1]
    return Schnittpunkt


def Brechung(L, v, n1, n2): # Snelliussches Brechungsgesetz
    alpha = angle(-L, v)
    beta = np.arcsin((n1/n2)*np.sin(alpha)) #-np.arcsin... für Linse
    Rotationsmatrix = np.array([
         [np.cos(beta), -np.sin(beta)],
         [np.sin(beta), np.cos(beta)]])
    Rotiert = np.matmul(Rotationsmatrix, -L) #-L für Linse
    return Rotiert #np.array([np.cos(beta), -np.sin(beta)])


def Kreis(x, q, r): # um die Linse zu malen
    return np.sqrt(r**2-(x-q)**2)

def S(x, v, e): # Schnittpunkt einer Geraden x + k v mit der zur y-Achse parallelen Gerade durch x = e 
     k = (e-x[0])/v[0]
     return x+k*v


K1_ls = np.linspace(2.0*brennweite-d1, 2.0*brennweite, num) # linspaces fuer die Linsen
K2_ls = np.linspace(2.0*brennweite, 2.0*brennweite+d2, num) 
    
def Linse(pos, v): # Berechnet den Strahlengang in einer Linse
    K1_result = SchnittpunktKreis(pos, v, q1, R1, -1.0) #Schnittpunkt mit dem 1. Kreis, bestimmt die Position der 1. Brechung
    x = K1_result[0] # der Stuetzvektor der Gerade wird auf den Schnittpunkt gesetzt 
    x0alt = x[0] # speichert die alte Position
    x1alt = x[1]
    W1_ls = np.linspace(0.0, x[0], num) # linspace fuer den Weg vom Objekt zum 1. Schnittpunkt
    b1 = pos[1] # Y-Achsen Abschnitt
    m1 = v[1]/v[0] # Steigunmg der 1. Gerade
    def yW1(x):
        return b1+m1*x
    L = K1_result[1] # Lotvektor zum ersten Schnittpunkt
    v = Brechung(L, v, n0, nD) #Brechung beim Eintritt
    K2_result = SchnittpunktKreis(x, v, q2, R2, +1.0) #Schnittpunkt mit dem 2. Lreis, bestimmt die Position der 2. Brechung
    x = K2_result[0] # der Stuetzvektor der Gerade wird auf den Schnittpunkt gesetzt  
    W2_ls = np.linspace(x0alt, x[0], num) # linspace fuer den Weg vom 1. zum 2. Schnittpunkt
    m2 = v[1]/v[0] # Steigung der 2. Gerade
    b2 = x1alt-m2*x0alt # Y-Achsen Abschnitt
    def yW2(x):
        return b2+m2*x
    x0alt = x[0] # speichert die alte Position
    x1alt = x[1]
    L = K2_result[1] # Lotvektor zum zweiten Schnittpunkt
    v = Brechung(-L, v, nD, n0) #Brechung beim Austritt, hier wird -L verwendet, weil sich der Lichtstrahl in die gleiche Richtung, wie L bewegt
    x = S(x, v, dist) # der Stuetzvektor der Gerade wird auf den Schnittpunkt mit dem Schirm gesetzt
    W3_ls = np.linspace(x0alt, x[0]*weitung, num) # linspace vom 2. Schnittpunkt zum Schnittpunkt mit dem Schirm
    m3 = v[1]/v[0] # Steigung der 3. Gerade
    b3 = x1alt-m3*x0alt # Y-Achsen Abschnitt 
    def yW3(x):
        return b3+m3*x
    plt.plot(W1_ls, yW1(W1_ls), 'k') # plottet die Geraden
    plt.plot(W2_ls, yW2(W2_ls), 'k')
    plt.plot(W3_ls, yW3(W3_ls), 'k')    


mEbene = np.array([1.0, 1.0]) #Richtungsvektor der "Ebene", da die Simulation im2-dimensionalen ist, ist es eine Gerade
bEbene = np.array([0.0, 0.0]) #Stuetzvektor der "Ebene"
Xachse = np.linspace(0.0, 100, num)
Yachse = mEbene[1]/mEbene[0]*Xachse+bEbene[1]

def Ebene(pos, v): #Diese Methode berechnet die Brechung an einer Ebenen, ist für die Linsensimulation unnoetig und wird nicht ausgefuehrt
    Schnittpunkt = SchnittpunktGerade(pos, v, bEbene, mEbene)
    Steigungswinkel = np.arctan(mEbene[1]/mEbene[0])
    Normalenwinkel = Steigungswinkel + np.pi/2.0
    Normalenvektor = np.array([np.cos(Normalenwinkel), np.sin(Normalenwinkel)])
    Gebrochen = Brechung(Normalenvektor, v, n0, nD)
    Brechungswinkel = np.arccos(Gebrochen[0])
    W1_ls = np.linspace(pos[0], Schnittpunkt[0], num)
    W2_ls = np.linspace(Schnittpunkt[0], 2*Schnittpunkt[0], num)
    yW1 = (v[1]/v[0])*W1_ls+pos[1]
    yW2 = np.tan(Brechungswinkel)*W2_ls+Schnittpunkt[1]-np.tan(Brechungswinkel)*Schnittpunkt[0]
    plt.plot(W1_ls, yW1, label=str(pos[1])+' 1') # plottet die Geraden
    plt.plot(W2_ls, yW2, label=str(pos[1])+' 2')

pos = np.array([0.0, 50.0]) # Festlegung mehrerer Startpunkte
pos2 = np.array([0.0, 25.0])
pos3 = np.array([0.0, 5.0])
pos4 = np.array([0.0, -5.0])
pos5 = np.array([0.0, -25.0])
pos6 = np.array([0.0,  -50.0])

xAxZ = np.linspace(0, 4.5*brennweite, num) #
yAxZ = 0.0*xAxZ # optische Achse

v2 = np.array([1.0, 0.0])

# Berechnung von Strahlengaengen
Linse(pos, v2)
Linse(pos2, v2)
Linse(pos3, v2)
Linse(pos4, v2)
Linse(pos5, v2)
Linse(pos6, v2)
plt.plot(K1_ls, Kreis(K1_ls, q1, R1), 'k') # Linsen nach oben
plt.plot(K2_ls, Kreis(K2_ls, q2, R2), 'k')
plt.plot(K1_ls, -Kreis(K1_ls, q1, R1), 'k') # Linsen nach unten
plt.plot(K2_ls, -Kreis(K2_ls, q2, R2), 'k')
plt.plot(xAxZ, yAxZ, label='optische Achse')
xa = np.array([0.0, 1.0*brennweite, 2.0*brennweite, 3.0*brennweite, 4.0*brennweite])
plt.xticks(xa, ['0', '1f', '2f', '3f', '4f'])

#plt.axis('equal')
plt.legend()
plt.show()
