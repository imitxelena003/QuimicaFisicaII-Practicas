from string import *
from numpy import *
import matplotlib
import matplotlib.pyplot as plt
    
def read_constants():
    global T
    global kb
    
    T=10.0
    kb=0.69503 #cm-1K
    
    return
    
def read_input():

    global Molekula
    global B
    global D
    global nue
    global Xe
    global alfae
    global JMax
    
    Molekula=input("Nombre de la molecula: ")
    B=input("Constante rotacional (cm-1): ")
    D=input("Distorsion centrifuga (cm-1): ")
    nue=input("Frecuencia vibracional (cm-1): ")
    Xe=input("Constante anarmonica: ")
    alfae=input("Constante de acoplamiento roto-vibracional (cm-1): ")
    JMax=input("Jmax: ")
    
    return
    
def FRot(J):

    FRot=B*J*(J+1) - D*(J*(J+1))**2

    return FRot

def GVib(v):

    GVib=nue*(v+0.5) - Xe*nue*(v+0.5)**2
    
    return GVib
        
def SRotVib(v,J):

    SRotVib = GVib(v) + FRot(J) - alfae*(v+0.5)*J*(J+1)
    
    return SRotVib
    
def NJ(J):

    NJ=(2.0*J+1.0)*exp(-1.0*FRot(J)/(kb*T))

    return NJ
    
def calc_qRot(JMax):
    
    qRot=0.0
    for J in range(0,JMax+1):
        qRot = qRot + (2.0*J+1.0)*exp(-1.0*FRot(J)/(kb*T))
        
    return qRot    
        
def calc_nu_RotVib_adarra(JMax,Adarra):

    # P Adarra hasten da J=1-n
    if Adarra == "P" :
       nu_Adarra=zeros([JMax],dtype=float64)
       #gogoratu range egingo du fortan do baino 1 gutxiago 
       for J in range(1,JMax+1):
             nu_Adarra[J-1]= SRotVib(1,J-1) - SRotVib(0,J)
             #print "nuP",J-1,nu_Adarra[J-1]
    elif Adarra == "Q":
       nu_Adarra=zeros([1],dtype=float64)
       nu_Adarra[0] = SRotVib(1,0) - SRotVib(0,0)
    elif Adarra == "R":
        # Q Adarra hasten da J=0-n
       nu_Adarra=zeros([JMax+1],dtype=float64)
       for J in range(0,JMax+1):
             nu_Adarra[J]= SRotVib(1,J+1) - SRotVib(0,J)
             #print "nuR",J,nu_Adarra[J]
       
    return nu_Adarra 

def calc_nu_RotVib(JMax):

    nu_RotVib=zeros([2*JMax+2],dtype=float64)
    
    nuP=calc_nu_RotVib_adarra(JMax,"P")
    nuQ=calc_nu_RotVib_adarra(JMax,"Q")
    nuR=calc_nu_RotVib_adarra(JMax,"R")
    
    for J in range(0,JMax):
        nu_RotVib[J]= nuP[JMax-J-1]
        print "nuP: ",J+1,nu_RotVib[J]
    nu_RotVib[JMax]=nuQ[0]
    for J in range(JMax+1,2*JMax+2):
        nu_RotVib[J]= nuR[J-JMax-1]
        print "nuR: ",J,nu_RotVib[J]

    return nu_RotVib

def generate_weights(JMax):

    qRot=calc_qRot(JMax)
    
    w_RotVib=zeros([2*JMax+2],dtype=float64)

#P Adarra
    for J in range(0,JMax):
        JIni=JMax-J
        w_RotVib[J]= NJ(JIni)/qRot
# Q Adarra
    w_RotVib[JMax]=0.0
# R Adarra
    for J in range(JMax+1,2*JMax+2):
        JIni=J-JMax-1
        w_RotVib[J]= NJ(JIni)/qRot     

    return w_RotVib
    
def create_spectra(nu,w):

#Originalmente se leia unal de y un f(i) pesos 
    sigma=1.0
    npuntos=10000
    wa=0.01
    numin=min(nu) -20.0
    numax=max(nu) +20.0
    
    rango=numax-numin
    paso=float(rango)/float(npuntos)
    
# ahora ponemos los pesos en 1
    xx=zeros([npuntos+1], dtype=float64)
    yy=zeros([npuntos+1], dtype=float64)
    
    f = open('./data', 'w')
    freq=numin
    for i in range(0,npuntos+1):
        freq=freq+paso*(1.0-wa/2.0)

        fw=0.0                   
        for k in range(0,size(nu)):

            fw=fw+w[k]/(1+4*(freq-nu[k])*(freq-nu[k])/(sigma*sigma))
        
        xx[i]=freq
        yy[i]=fw
 #      print freq,fw   

        f.write(str(freq) + " " + str(fw) + "\n")

    f.close()
    make_plot1(xx,yy)

    return
  
def make_plot(xx,yy):
    
    matplotlib.rcParams['axes.unicode_minus'] = False
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(xx[:], yy[:])
    ax.set_title("Espectro de rotacion-vibracion")
    plt.legend((kalkulua, 'Total message length'),
           'upper left', shadow=True, fancybox=True)
    plt.show()
    
    return
    
def make_plot1(xx,yy):
    from pylab import *

    font = {'family'     : 'serif',
        'color'      : 'r',
        'weight' : 'normal',
        'size'   : 12,
        }

    
    plot(xx, yy, 'b')
    title('Espectro de rotacion-vibracion de la molecula de '+ Molekula, font, size='large', color='r')
    text(xx[10], 0.75*max(yy),"B = " + str(B) + " cm"r'$^{-1}$'"\n" + \
                              "D = " + str(D) + " cm"r'$^{-1}$'"\n" + \
                              r'$\nu_e = $' + str(nue) + " cm"r'$^{-1}$'"\n" + \
                              r'$\chi_e = $' + str(Xe) + "\n" + \
                              r'$\alpha_e = $' + str(alfae) \
                              , color='k')


    xlabel('Frecuencia (cm-1)', font, style='italic')
    ylabel('Intensidad relativa', font)

    show()

    return

##### Programa Principal #######

read_constants()
read_input()

nu_RotVib=calc_nu_RotVib(JMax)
w_RotVib=generate_weights(JMax)
create_spectra(nu_RotVib,w_RotVib)
