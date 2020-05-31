import numpy as np
import matplotlib.pyplot as plt
import cv2

filename = 'Fehlt.png' # Bild, welches verarbeitet wird

XSIZE = 50 # Groesse des Bildes
YSIZE = 50

AnzahlPixelSchirmX = 2000
AnzahlPixelSchirmY = AnzahlPixelSchirmX
PixelBreite = 0.05

gamma = np.pi/1000.0 # Oeffnungswinkel der Lichtkegel, die fuer jeden Pixel erzeugt werden, in RAD
direc = np.array([0.0, 0.0, 1.0]) # Haupt(bzw. Rotations)achse des Photonenkegels

R1 = 800.0 #Radius der ersten Kugel
R2 = R1 #Radius der zweiten Kugel

e1 = np.array([1.0, 0.0, 0.0]) # Richtungsvektoren der Schirmebene
e2 = np.array([0.0, 1.0, 0.0])
factor_num = 1 #Es werden factor_num*H_ij viele Photonen erzeugt


n0 = 1.0 #Brechungsindices
nD = 1.5

h = 2.0*((0.5*XSIZE)**2+(0.5*YSIZE)**2)**0.5 # Die Hoehe der Linse haengt von der Groesse des Bildes ab, um sphaerische Aberration zu vermeiden wird mit 2 multipliziert
d1 = (R1-(R1**2-h**2)**0.5) # Jeweilige Durchmesser Haelften um gegebene Hoehe zu erreichen
d2 = (R2-(R2**2-h**2)**0.5)
d = d1+d2 #Gesamt Durchmesser der Linse

Radius_Blende = 0.6*h #Radius der Blende

def Linsenschleifergleichung(): #Linsenschleifergleichung zur Berechnung der theoretischen Brennweite
    D = (nD-1.0)*(1.0/R1 - 1.0/(-R2) + d*(nD-1.0)/(nD*R1*(-R2)))
    return 1.0/D

brennweite = Linsenschleifergleichung()*0.3334 #0.3334 ist ein aus unerfindlichen Gruenden noetiger Korrekturfaktor

print(h)
K1_vec = np.array([0, 0, 2*brennweite+R1-d1]) #Vektor des Kugelmittelpunkts der ersten  Kugelmittelpunkts
K2_vec = np.array([0, 0, 2*brennweite-R2+d2]) #Vektor des Kugelmittelpunkts der zweiten Kugelmittelpunkts

Bild = [] #Array, das die Schirme, bei verschiedenen Abstaenden zur Linse speichert
e0 = [] #Array, welches die verschiedenen Abstaende der Schirme zur linse speichert, bzw. die Stuetzvektoren der Schirmebenen

#Bild.append(np.zeros([AnzahlPixelSchirmX, AnzahlPixelSchirmX]))
#e0.append(np.array([0.0, 0.0, 3.0*brennweite])) #Positioniert die Schirmebene in 3*f zur Objektebene, also in 1*f zur Linse
#Bild.append(np.zeros([AnzahlPixelSchirmX, AnzahlPixelSchirmX]))
#e0.append(np.array([0.0, 0.0, 3.5*brennweite])) #Positioniert die Schirmebene in 3.5*f zur Objektebene, also in 1.5*f zur Linse
Bild.append(np.zeros([AnzahlPixelSchirmX, AnzahlPixelSchirmX]))
e0.append(np.array([0.0, 0.0, 4.0*brennweite])) #Positioniert die Schirmebene in 4*f zur Objektebene, also in 2*f zur Linse

AS = 3 #Bei manchen Ausgaben, muss gerundet werden


def cut(x): #Rundung bzw. Wegschneiden der letzten Stellen
    return ""+str(int(x*10.0**AS)/10.0**AS)


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

def angle(x, y): # Winkel zwischen zwei Vektoren
    return np.arccos(np.abs(sp(x, y))/(norm(x)*norm(y)))

def abc(a, b, c, index): #pq bzw. abc Formel. Der index gibt an welche der beiden Loesungen gewollt ist
    p = b/a
    q = c/a
    return -0.5*p + index * np.sqrt(0.25*p**2 - q) #[-0.5*p, (0.25*p**2-q)**0.5]

def VektorErzeugungGleichmaessig(phi, n): # erzeugt einen Vektor im Winkelbereich beta <= phi um den Vektor n, entstehende Vektoren sind auf der Kegeloberflaeche gleichmaessig verteilt
    if (n[0]==n[1]==n[2]):
        ortho = np.array([np.sqrt(0.5), -np.sqrt(0.5), 0])
    else:
        ortho_ = np.array([n[1]-n[2], n[2]-n[0], n[0]-n[1]])
        ortho = ortho_*1.0/norm(ortho_)
        #ortho = np.array([n[1]-n[2], n[2]-n[0], n[0]-n[1]])
        #ortho = ortho * 1/np.sqrt(ortho[0]**2 + ortho[1]**2 + ortho[2]**2)
    if (phi < np.pi/2):
        sinMax = np.sin(phi)
    else:   
        sinMax = 1
    beta = phi * np.random.random_sample()
    y = np.random.random_sample()*sinMax
    while (y>np.sin(beta)):
       beta = phi * np.random.random_sample()
       y = np.random.random_sample()*sinMax 
    gekippt = np.cos(beta)*n + np.sin(beta)*ortho
    alpha = np.random.random_sample()*2.0*np.pi
    gedreht = Rot(n,  alpha, gekippt)
    return gedreht

def VekttorErzeugung(phi, n): # erzeugt einen Vektor im Winkelbereich beta <= phi um den Vektor n
    if (n[0]==n[1]==n[2]):
        ortho = np.array([np.sqrt(0.5), -np.sqrt(0.5), 0])
    else:     
        ortho = np.array([n[1]-n[2], n[2]-n[0], n[0]-n[1]])
        ortho = ortho * 1/np.sqrt(ortho[0]**2 + ortho[1]**2 + ortho[2]**2)
    
    phi *= np.random.random_sample()
    gekippt = np.cos(phi)*n + np.sin(phi)*ortho
    
    alpha = np.random.random_sample()*2.0*np.pi
    gedreht = Rot(n,  alpha, gekippt)
    
    return gedreht

def Rot(n, alpha, x): # x wird um den Winkel alpha um den Vektor n gedreht
    c = np.cos(alpha) # Um Rechenzeit zu sparen, werden die Ergebnisse trigonometrischer Funktionen gespeichert
    s = np.sin(alpha)
    R = np.array([
        [n[0]*n[0]*(1-c)+c,      n[0]*n[1]*(1-c)-n[2]*s, n[0]*n[2]*(1-c)+n[1]*s],
        [n[1]*n[0]*(1-c)+n[2]*s, n[1]*n[1]*(1-c)+c,      n[1]*n[2]*(1-c)-n[0]*s],
        [n[2]*n[0]*(1-c)-n[1]*s, n[2]*n[1]*(1-c)+n[0]*s, n[2]*n[2]*(1-c)+c]
        ])
    r = np.matmul(R, x)
    return r

def SchnittpunktKugel(a, b, z, r, index): # Bestimmt den Schnittpunkt zwischen einer Gerade a+k*b und einer Kugel des Radiuses r um den Punkt z 
    delta = a-z # Abstand zwischen Stuetzvektor und Zentrum der Kugel
    C1 = sp(b, b) # Skalarprodukt von b und b
    C2 = 2.0*sp(b, delta)
    C3 = sp(delta, delta)-r**2
    k = abc(C1, C2, C3, index) # Faktor k in der Geradengleichung a+k*b
    S = a+k*b # Schnittpunkt a+k*b
    L = S-z # Lot, bzw Vektor von Zentrum zu Schnittpunkt
    return np.array([S, L])

def SchnittpunktEbene(g0, g1, q): #Vereinfachte intersection Methode, bei der die Ebene orthogonal zur z-Achse ist und die z-Komponte q hat
    k = (q-g0[2])/g1[2]
    return g0+k*g1
  
  
def Brechung(L, x, n1, n2): # Berechnet den Einfallsvektor nach der Brechung an einer Ebene : Lot, Einfallsvektor, Brechzahlen
    z_ = -vp(-L, x)
    z = 1/norm(z_)*z_
    alpha = angle(-L, x)
    beta  = np.arcsin(n1*np.sin(alpha)/n2) #n1 sin alpha = n2 sin beta -> n1/n2 sin alpha = sin beta
    y_ = Rot(z, beta, -L)
    y = y_*1.0/norm(y_)
    return y
  

def Schirm(x, w): # Berechnet den Schnittpunkt eines Photons mit den Schirmebenen und addiert dann 1 zum Pixel auf, der getroffen wird
  global Bild
  for k in range(0, len(Bild)): # Geht durch alle Schirmebenen durch
      v = SchnittpunktEbene(x, w, e0[k][2])
      if(-PixelBreite*AnzahlPixelSchirmX*0.5 <= v[0] and v[0] <= PixelBreite*AnzahlPixelSchirmX*0.5 and -PixelBreite*AnzahlPixelSchirmY*0.5 <= v[1] and v[1] <= PixelBreite*AnzahlPixelSchirmY*0.5):
        x_ = v[0]+0.5*PixelBreite*AnzahlPixelSchirmX
        y_ = v[1]+0.5*PixelBreite*AnzahlPixelSchirmY
        i = int(x_/PixelBreite)-1
        j = int(y_/PixelBreite)-1
        Bild[k] [i] [j] += 1
      else: pass

def Linse(pos, v): # berechnet das Verhalten eines Lichtstrahls in einer Linse
    K1_result = SchnittpunktKugel(pos, v, K1_vec, R1, -1.0) #Schnittpunkt mit der 1. Kugel, bestimmt die Position der 1. Brechung
    x = K1_result[0]
    #x = pos
    if(x[0]**2+x[1]**2 <= Radius_Blende**2): #Wenn die Photonen zu weit von der z-Achse entfernt sind, werden sie von der Blende absorbiert
        L = K1_result[1]
        v = Brechung(L, v, n0, nD) #Brechung beim Eintritt
        K2_result = SchnittpunktKugel(x, v, K2_vec, R2, +1.0) #Schnittpunkt mit der 2. Kugel, bestimmt die Position der 2. Brechung
        x = K2_result[0]
        L = K2_result[1]
        v = Brechung(L, v, nD, n0) #Brechung beim Austritt
        Schirm(x, v)
        return [x, v]

def createall(filename): # Erzeugung der Photonen eines Bildes und Berechnung ihrer Strahlengaenge 
    img = cv2.imread(filename,0)
    #print(img)
    h = len(img[0])
    w = len(img)
    for i in range(w): # Durchgehen aller Pixel
        print(str(i/(1.0*w))) # Fortschrittsangabe
        for j in range(h):
            num = img[i, j]*factor_num # Anzahl an zu erzeugender Lichtstrahlen fuer den ij-ten Pixel
            pos = np.array([(i-w/2.0)*(1.0*XSIZE/w), (j-h/2.0)*(1.0*YSIZE/w), 0.0]) # Stuetzvektor der Lichtstrahlen, fuer den ij-ten Pixel
            for k in range(num):
                v = VektorErzeugungGleichmaessig(gamma, direc)
                Linse(pos, v)
                
                
                #counting_photons(pos, v)

#Ausfuehrung der benoetigten Methoden  
createall(filename)

#Speicherung der Bilddateien
for i in range(len(Bild)):
    plt.imshow(Bild[i])
    #plt.savefig('PL2-Rakete-1000-x'+str(factor_num)+'-OhneLinse-PB,'+str(PixelBreite)+'-APS,'+str(AnzahlPixelSchirmX)+'-d,'+str(e0[i][2]/brennweite-2.0)+'.pdf')
    plt.savefig('PL2-Rakete-1000-x'+str(factor_num)+'-R='+str(R1)+'-PB,'+str(PixelBreite)+'-APS,'+str(AnzahlPixelSchirmX)+'-d,'+str(e0[i][2]/brennweite-2.0)+'.pdf')


