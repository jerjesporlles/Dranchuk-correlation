########### Correlación de Dranchuk ########### 
import math

#Datos
GamaG = 0.7  #adim
P = 500      #psia
T = 150      #F
Zsup = 1     #adim
n = 1        #contador
Zrr = 1      #error de Z

#T y P Pseudocríticas
Tpc = 167 +316.67 * GamaG
Ppc = 702.5 - 50 * GamaG

print("Tpc=",Tpc,"[°F]     Ppc=",Ppc,"[psia]")

#Tpr y Ppr
Tpr = (T + 459.67)/ Tpc
Ppr = P / Ppc
DenR = 0.27 * Ppr / (Zsup * Tpr)

print("Tpr=",round(Tpr,5),"[adim]   Ppr=",round(Ppr,5),"[adim]", "Densidad R=",round(DenR,4))

#Constantes
A1 =  0.31506237
A2 = -1.0467099
A3 = -0.57832729
A4 =  0.53530771
A5 = -0.61232032
A6 = -0.10488813
A7 =  0.68157001
A8 =  0.68446549

#Dranchuk Cálculo
while (Zrr > 0.000005):   #indicamos Tolerancia
    Zparte1 = (A1 + A2/Tpr + A3/Tpr**3) * DenR
    Zparte2 = (A4 + A5/Tpr) * DenR**2
    Zparte3 = (A5 * A6 * DenR**5) / Tpr
    Zparte4 = (A7 * DenR**2 / Tpr**3) * (1 + A8 * DenR**2) * math.exp(-A8 * DenR**2)

    Z = 1 + Zparte1 + Zparte2 + Zparte3 + Zparte4

    Zrr = Zsup - Z
    Zsup = Z
    print("Nueva Z: ",round(Z,7),"El error es de:",round(Zrr,7),"Iteración: ",n)
    n += 1
    DenR = 0.27 * Ppr / (Zsup * Tpr)

print("Cálculo con Método de Dranchuk para factor de compresibilidad = ", round(Z,6))
