from scipy import *
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


t,p_a, p_b, T_2, T_1, P = loadtxt('Messwerte/mess', unpack = 'true')

data =array( [[t[i],p_a[i], p_b[i], T_1[i], T_2[i], P[i]] for i in range(0,33)])
t *= 60
print(make_LaTeX_table(data,['$t [s]$','$p_a [10^5Pa]$','$p_b [10^5Pa]$','$T_a [^\circ]$','$T_b [^\circ]$','P [W]']))


def T_quad(x,A,B,C):
    return A*x**2+B*x+C
params, cov = curve_fit(T_quad,t,T_1)

A=[]
B=[]
C=[]
A.append(ufloat(params[0],sqrt(cov[0][0])))
B.append(ufloat(params[1],sqrt(cov[1][1])))
C.append(ufloat(params[2],sqrt(cov[2][2])))

params, cov = curve_fit(T_quad,t,T_2)
A.append(ufloat(params[0],sqrt(cov[0][0])))
B.append(ufloat(params[1],sqrt(cov[1][1])))
C.append(ufloat(params[2],sqrt(cov[2][2])))

dT_dt = [[0,0],[0,0],[0,0],[0,0],[0,0]]
for i in range(0,5):
    for j in range(0,2):
       dT_dt[i][j] = ( 2*60*(i+1)*A[j]+ B[j])

print(dT_dt)
x = linspace(0,t[-1])
print(make_LaTeX_table(array(dT_dt),['T1','T2']))

plot(t, T_2, 'bo' )
plot(t,T_1,'ro')


show()