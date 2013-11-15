# Using the magic encoding
# -*- coding: utf-8 -*-
from pylab import *
from uncertainties import  *
from scipy.optimize import curve_fit
def make_LaTeX_table(data,header, flip= 'false', onedim = 'false'):
    output = '\\begin{tabular}{'
    #Get dimensions
    if(onedim == 'true'):
        if(flip == 'false'):
        
            data = array([[i] for i in data])
        
        else:
            data = array([data])
        

    row_cnt, col_cnt = data.shape
    header_cnt = len(header)
    
    if(header_cnt == col_cnt and flip== 'false'):
        #Make Format
        output += '|'
        for i in range(col_cnt):
            output += 'c|'
        output += '}\n\\hline\n'+ header[0]
        for i in range (1,col_cnt):
            output += ' & ' + header[i]
        output += ' \\\\\n\\hline\n'
        for i in data:
            output += str(i[0])
            for j in range(1,col_cnt):
                output += ' & ' + str( i[j])
            output += '\\\\\n'
        output += '\\hline\n\\end{tabular}\n'
                            
        return output
    else:
        if(row_cnt == header_cnt):
            output += '|c|' + (col_cnt)*'c' + '|}\n\\hline\n'
            for i in range(row_cnt):
                output += header[i]
                for j in range(col_cnt):
                    output += ' & ' + str(data[i][j])
                output += '\\\\\n\\hline\n'
                
            output += '\\end{tabular}\n'
            return output
        else:
            return 'ERROR'

    
def err(data):
    mean = data.mean()
    N = len(data)
    err = 0
    for i in data:
        err += (i - mean)**2
    err = sqrt(err/((N-1)*N))
    return ufloat(mean,err)

def lin_reg(x,y):
    N = len(x)
    sumx = x.sum()
    sumy = y.sum()
    sumxx = (x*x).sum()
    sumxy = (x*y).sum()
    m = (sumxy -  sumx*sumy/N)/(sumxx- sumx**2/N)
    b = sumy/N - m*sumx/N
    
    sy = sqrt(((y - m*x - b)**2).sum()/(N-1))
    m_err = sy *sqrt(N/(N*sumxx - sumx**2))
    b_err= m_err * sqrt(sumxx/N)
    return ufloat(m,m_err), ufloat(b,b_err)


def T_quad(x,A,B,C):
    return A*x**2+B*x+C



#Messwerrte einlesen
t,p_a, p_b, T_2, T_1, P = loadtxt('Messwerte/mess', unpack = 'true')

t *= 60 # Messzeit in sec umwandeln
T_1 += 273.15 # Temperaturen in Kelvin umechenen
T_2 += 273.15 
p_a = (p_a+1)*1e5 
p_b = (p_b+1)*1e5 
data =array( [[t[i],p_a[i], p_b[i], T_1[i], T_2[i], P[i]] for i in range(0,33)])


print("Tabelle aller Messwerte")
print(make_LaTeX_table(data,['$t [s]$','$p_a [10^5Pa]$','$p_b [10^5Pa]$','$T_a [^\circ]$','$T_b [^\circ]$','P [W]']))

## Beide Datensetzeand a*t^2+b*t+c fitten und Fehler aus der Kovarianzmatrix berechnen
params, cov = curve_fit(T_quad,t,T_1)
A_1 = ufloat(params[0],sqrt(cov[0][0]))
B_1 = ufloat(params[1],sqrt(cov[1][1]))
C_1 = ufloat(params[2],sqrt(cov[2][2]))

params, cov = curve_fit(T_quad,t,T_2)
A_2 = ufloat(params[0],sqrt(cov[0][0]))
B_2 = ufloat(params[1],sqrt(cov[1][1]))
C_3 = ufloat(params[2],sqrt(cov[2][2]))

#Plot Data und Fit
x = linspace(0,33*60) # Werte fuer die Interpolationsfunktionen
#T1
plot(t,T_1,'ro', label= r'Messwerte $T_1$')
plot(x,T_quad(x,A_1.n,B_1.n,C_1.n),'r', label = u'Quadratischer Fit für $T_1$')

#T2 Das Gleoche in blau
plot(t,T_2,'bo', label= r'Messwerte $T_2$')
plot(x,T_quad(x,A_2.n,B_2.n,C_3.n),'b', label = u'Quadratischer Fit für $T_2$')

xlabel(u"Zeit [$s$]")
ylabel(u"Temperatur [$^\circ C$]")
legend(loc = 'upper left')


#Abeltungen definierent

dT1_dt = array([2*A_1 * t[i] + B_1 for i in [8,16,24,32]])
dT2_dt =array([ 2*A_2 * t[i] + B_2 for i in [8,16,24,32]])
T1_ = array([T_1[i] for i in [8,16,24,32]])
T2_ = array([T_2[i] for i in [8,16,24,32]])
pb_ = array([p_b[i] for i in [8,16,24,32]])
pa_ = array([p_a[i] for i in [8,16,24,32]])
t_ =60*  array([8,16,24,32])
P_ = [ P[i]for i in [8,16,24,32]]
#4 Zeiten fuer die Weitere Analyse auswaehlen und Ableitungen berechen
#und in Latex tabelle darstellen

data = array([t_, T1_,dT1_dt,T2_, dT2_dt ,P_])
header = [r"$t [s]$",r"$T_1 [^\cric C]$",r"$\frac{dT_1}{dt} [^\cric Cs^{-1}]$",r"$T_2 [^\cric C]$",r"$\frac{dT_1}{dt} [^\cric Cs^{-1}]$", r"$P [W]$"]
print(data)

print("Tabelle ausgewaehlter Messwerte mit Temperatur Ableitungen")
output = make_LaTeX_table(data.T, header)
print(output)

## Guetezaheln berechen fuer die Verschieden Zeitupunkte
rho_w = ufloat(0.99799,0) # kg/L Dichte von Wasser bei 20C Teubner phys. Praktikum
c_w =ufloat(4180,2)  # Spezifische Waerme Kapazitaet des wasser 20-40 C
V = ufloat(4,0.01)
R = 8.3143 # molare Gaskonstante
nu_real = ( V*rho_w*c_w + 750)* dT1_dt /P_
nu_id   =[round(T1_[i]/(T1_[i]-T2_[i]),2) for i in range(4)]

#Ausgabe der Tabelle

data = array([t_,T1_-T2_, nu_real, nu_id])
header = [r't [s]',r'$\Delta t [^\circ C]$',r'$\nu_{real}$',r'$nu_{id}$']
print("Bestimmung der Guetezahlen")
print(make_LaTeX_table(data.T,header))

# Erstellung des p-t Diagramms
close() # Alten plot schleissen

# Druecke und Temperaturen zusammenfuegen

p_ges = array(list(p_b[3:])+list(p_a[6:]))
T_ges = array(list(T_1[3:])+list(T_2[6:]))
plot(T_ges,p_ges, 'bx')

def exp_fit(x,A,B):
    return B*exp(-A/(R*x))
    
params , cov = curve_fit(exp_fit,T_ges,p_ges)

L = ufloat(params[0],sqrt(cov[0][0]))
print("Die Verdampfungswaereme des Arbeitsmediums betraegt : %s" % L)


x= linspace(270,330) # Temperaturbereich fuer die Ausgleichskurve
plot(x,exp_fit(x,params[0],params[1]))
show()
#Bestimmung des Massendurchsatz 

dm_dt = -(V*rho_w*c_w+750)*dT2_dt/L
data = array([t_, dm_dt])
header = ["$t [s]$",r"$\frac{dm}{dt} [\frac{kg}{s}]$"]
print("Massendursatz der WP:")
print(make_LaTeX_table(data.T, header))


#Bestimmung der mech. Kompressorleistung
#Zunaechst muss die Dichte des Arbeitsmediums mit Hilfe der idealen Gas Gl.
#bestimmt werden : PV/T =const


rho_km0 = 5.51 #kg/m^3 vom Blatt
T_0 = 271.15 # K  Standardbed. 0C
kappa = 1.14 # Adiabatenkoeff. vom Blatt
p_0 = 1e5 # 1 Bar Normaldruck
rho_km = T2_*pa_*rho_km0/(p_0*T_0)

N = ((pb_ * (pa_/pb_)**(1/kappa)-pa_)*dm_dt)/((kappa-1)*rho_km)

print(N)