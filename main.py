from sympy import Symbol, diff, exp,Abs
from scipy.optimize import minimize
import scipy
import numpy as np
from scipy import integrate
from scipy.misc import derivative
import matplotlib.pyplot as plt
import csv
#sympy.init_printing(use_latex="mathjax")

def function ():

    print("Hello world")

function()



b = 3.14

c = 0.85

d = 0.45

h_bar = 6.582119 * 10**-16 #eV*s

m_c_squared = 939*10**6 #eV

c_squared =  8.98755179*10**16 #m**2/s**2

light = 299792458 #m/s

pi = 3.14159265359

A1 = -6.28

A2 = -6.28

h_c = 1970#A*eV

stala_hamilton = -7.76*10**-5#((h_c**2)/(2*m_c_squared))
#print(stala_hamilton)
R = Symbol("R")

B = Symbol("B")

C = Symbol("C")

D = Symbol("D")

x = Symbol("x")

C1 = Symbol("C1")

C2 = Symbol("C2")

C3 = Symbol("C3")

C4 = Symbol("C4")

X1 = Symbol("x1")

X2 = Symbol("x2")

#a1 = Symbol("a1")

#A1 = Symbol("A1")

#a2 = Symbol("a2")

#A2 = Symbol ("A2")

#k1 = Symbol("k1") #na razie przyjmujemy że to sqrt(2mE)/h

#k2 = Symbol("k2")

E = Symbol("E")

#alpha1 = Symbol("alpha1")#sqrt(mk/pi**2)

#alpha2 = Symbol("alpha2")#sqrt(mk/pi**2)

psi = Symbol("psi")





#narazie przyjmujemy ,że alfa1 alfa2 a1 a2 =1 , a polozenia studni jako 1 i 3




a1 = 1

a2 = 1

XWellA = 0.001

XWellB = 0.501

m_p = 1,67*10**-27
#######################################################

#Liczmy alfy

#omega = ((-A1*a1*c_squared/m_c_squared)**1/2)
#omega = (10**-29 / 10**-27)**1/2
#print(omega)
#k = (m_c_squared*omega/h_bar*light)**1/2

#alpha1 = (m_c_squared*k/(light**2 *pi**2))**1/4
#alpha2 = alpha1
#alpha1 = (m_c_squared*omega/h_c)**1/2
alpha1 = 0.7
alpha2 = alpha1
print(alpha1)

alpha2 = alpha1

def U(variable,potential_well): # U w zależności od studni i atomu

    if potential_well == "a":

        return ((alpha1**(1/2))/(pi**(1/4)))*np.exp((-1/2)*alpha1**2*(variable-XWellA)**2)

    else:

        return ((alpha2 ** (1 / 2)) / (pi ** (1 / 4))) * np.exp((-1 / 2) * alpha2 ** 2 * (variable - XWellB)**2)





def double_derivative(formula , variable): #do hamiltonianu

    outcome = ((-h_c**2)/(2*m_c_squared))*diff(formula,variable,variable)

    return outcome

internal_potential = B*(exp((-2*(R-C))/D) - 2*exp((-(R-C))/D)) #potencjał wewnętrzny

print(internal_potential)

internal_potential = internal_potential.subs(B,b) #podstawiam podane wartości

internal_potential = internal_potential.subs(C,c)

internal_potential = internal_potential.subs(D,d)

internal_potential = internal_potential.subs(R, Abs(X1-X2))

print(internal_potential)

def internal_pot(x,y):
    R=abs(x-y)
    #print(R)
    return 3.14*(np.exp((-2 * (R - 0.85)) / 0.45) - 2 * np.exp((-(R - 0.85)) / 0.45))
    #print('pot wew')
#print(scipy.integrate.dblquad(internal_pot,XWellA-3,XWellA+3,lambda x:XWellB-3,lambda x:XWellB+3))
#potencjal zewnetrzny
def ext_potential(variable):
    V_ext = A1*np.exp(-(a1*(variable-XWellA)**2)/2)+A2*np.exp(-(a2*(variable-XWellB)**2)/2)
    #V_ext = V_ext.subs(A1,-6.28)
    #V_ext = V_ext.subs(A2,-6.28)
    #V_ext = V_ext.subs(a1,1)
    #V_ext = V_ext.subs(a2,1)
    return V_ext
#od razu suma V_ext(x1) i V_ext(x2)
#print('ext pot calka')
#print(scipy.integrate.quad(ext_potential,XWellA-3,XWellA+3))




#hamiltonian





#wymnożone wszytkie U -- funkcja psi

#test = C1*U(X1,"a")*U(X2,"b")+C1,"a")*U(X1,"b")+C4*U(X2,"a")*U(X2,"b") #proszę sprawdzić czy dobrze

#psi --nie usuwaj mi tych tablic z psi

#psi_tab =[U(X1,"a")*U(X2,"b"),U(X2,"a")*U(X1,"b"),U(X1,"a")*U(X2,"a"),U(X1,"b")*U(X2,"b")]
def psi_w_a(x):
    return ((alpha1 ** (1 / 2)) / (pi ** (1 / 4)))* np.exp(-(alpha1**2)*((x - XWellA)**2) / 2)
def psi_w_b(x):
    return ((alpha2 ** (1 / 2)) / (pi ** (1 / 4)))* np.exp(-(alpha2**2)*((x - XWellB)**2) / 2)



def psiAB(x,y):
    return psi_w_a(x)*psi_w_b(y)
def psiBA(x,y):
    return psi_w_a(y)*psi_w_b(x)

def psiAA(x,y):
    return psi_w_a(x)*psi_w_a(y)

def psiBB(x,y):
    return psi_w_b(x)*psi_w_b(y)



def hamiltonianAA(x,y):
    return stala_hamilton*derivative(psi_w_a,x,n=2,dx=10**-6)*psi_w_a(y)+stala_hamilton*derivative(psi_w_a,y,n=2,dx=10**-6)*psi_w_a(x)+internal_pot(x,y)+ext_potential(x)+ext_potential(y)
def hamiltonianBB(x,y):
    return stala_hamilton*derivative(psi_w_b,x,n=2,dx=10**-6)*psi_w_b(y)+stala_hamilton*derivative(psi_w_b,y,n=2,dx=10**-6)*psi_w_b(x)+internal_pot(x,y)+ext_potential(x)+ext_potential(y)
def hamiltonianAB(x,y):
    return stala_hamilton*derivative(psi_w_a,x,n=2,dx=10**-6)*psi_w_b(y)+stala_hamilton*derivative(psi_w_b,y,n=2,dx=10**-6)*psi_w_a(x)+internal_pot(x,y)+ext_potential(x)+ext_potential(y)
def hamiltonianBA(x,y):
    return stala_hamilton*derivative(psi_w_b,x,n=2,dx=10**-6)*psi_w_a(y)+stala_hamilton*derivative(psi_w_a,y,n=2,dx=10**-6)*psi_w_b(x)+internal_pot(x,y)+ext_potential(x)+ext_potential(y)



#i =1 AB; i =2 BA; i = 3 AA ;i =4 BB
b = np.array([[0., 1.], [1., 1.]])
def S11 (x,y):
    return psiAB(x,y)*psiAB(x,y)
def S12 (x,y):
    return psiAB(x,y)*psiBA(x,y)
def S13 (x,y):
    return psiAB(x,y)*psiAA(x,y)
def S14 (x,y):
    return psiAB(x,y)*psiBB(x,y)
def S21 (x,y):
    return psiBA(x,y)*psiAB(x,y)
def S22 (x,y):
    return psiBA(x,y)*psiBA(x,y)
def S23 (x,y):
    return psiBA(x,y)*psiAA(x,y)
def S24 (x,y):
    return psiBA(x,y)*psiBB(x,y)
def S31 (x,y):
    return psiAA(x,y)*psiAB(x,y)
def S32 (x,y):
    return psiAA(x,y)*psiBA(x,y)
def S33 (x,y):
    return psiAA(x,y)*psiAA(x,y)
def S34 (x,y):
    return psiAA(x,y)*psiBB(x,y)
def S41 (x,y):
    return psiBB(x,y)*psiAB(x,y)
def S42 (x,y):
    return psiBB(x,y)*psiBA(x,y)
def S43 (x,y):
    return psiBB(x,y)*psiAA(x,y)
def S44 (x,y):
    return psiBB(x,y)*psiBB(x,y)

def h11 (x,y):
    return psiAB(x,y)*hamiltonianAB(x,y)
def h12 (x,y):
    return psiAB(x,y)*hamiltonianBA(x,y)
def h13 (x,y):
    return psiAB(x,y)*hamiltonianAA(x,y)
def h14 (x,y):
    return psiAB(x,y)*hamiltonianBB(x,y)
def h21 (x,y):
    return psiBA(x,y)*hamiltonianAB(x,y)
def h22 (x,y):
    return psiBA(x,y)*hamiltonianBA(x,y)
def h23 (x,y):
    return psiBA(x,y)*hamiltonianAA(x,y)
def h24 (x,y):
    return psiBA(x,y)*hamiltonianBB(x,y)
def h31 (x,y):
    return psiAA(x,y)*hamiltonianAB(x,y)
def h32 (x,y):
    return psiAA(x,y)*hamiltonianBA(x,y)
def h33 (x,y):
    return psiAA(x,y)*hamiltonianAA(x,y)
def h34 (x,y):
    return psiAA(x,y)*hamiltonianBB(x,y)
def h41 (x,y):
    return psiBB(x,y)*hamiltonianAB(x,y)
def h42 (x,y):
    return psiBB(x,y)*hamiltonianBA(x,y)
def h43 (x,y):
    return psiBB(x,y)*hamiltonianAA(x,y)
def h44 (x,y):
    return psiBB(x,y)*hamiltonianBB(x,y)
#Test potencjału zewnętrznego, wychodzi OK
'''
y = []
x = []
for i in np.arange(0,1000,1):
    x.append(i/500)
    y.append(ext_potential(i/500))

plt.figure(1)
plt.scatter(x,y)
plt.show()
'''
#Test potencjalu wewnetrzengo, jest OK
'''
y=[]
x=[]
for R in np.arange(0,1,0.001):
    x.append(R)
    y.append(internal_pot(0,R))
plt.figure(1)
plt.scatter(x,y)
plt.axhline(y=0)
plt.show()
'''

x_axis = []
C1_values = []
C2_values = []
C3_values = []
C4_values = []
E_values = []
for R in np.arange(0.5,3.5,0.1):

    print("liczę S11")
    w_S11 = scipy.integrate.dblquad(S11,XWellA-3,XWellA+3,lambda x:R-3, lambda x:R+3)
    print("liczę S12")
    w_S12 = scipy.integrate.dblquad(S12,XWellA-3,XWellA+3,lambda x:R-3, lambda x:R+3)
    print("liczę S13")
    w_S13 = scipy.integrate.dblquad(S13,XWellA-3,XWellA+3,lambda x:R-3, lambda x:R+3)
    print("liczę S14")
    w_S14 = scipy.integrate.dblquad(S14,XWellA-3,XWellA+3,lambda x:R-3, lambda x:R+3)
    print("liczę S21")
    w_S21 = scipy.integrate.dblquad(S21,XWellA-3,XWellA+3,lambda x:R-3, lambda x:R+3)
    print("liczę S22")
    w_S22 = scipy.integrate.dblquad(S22,XWellA-3,XWellA+3,lambda x:R-3, lambda x:R+3)
    print("liczę S23")
    w_S23 = scipy.integrate.dblquad(S23,XWellA-3,XWellA+3,lambda x:R-3, lambda x:R+3)
    print("liczę S24")
    w_S24 = scipy.integrate.dblquad(S24,XWellA-3,XWellA+3,lambda x:R-3, lambda x:R+3)
    print("liczę S31")
    w_S31 = scipy.integrate.dblquad(S31,XWellA-3,XWellA+3,lambda x:R-3, lambda x:R+3)
    print("liczę S32")
    w_S32 = scipy.integrate.dblquad(S32,XWellA-3,XWellA+3,lambda x:R-3, lambda x:R+3)
    print("liczę S33")
    w_S33 = scipy.integrate.dblquad(S33,XWellA-3,XWellA+3,lambda x:R-3, lambda x:R+3)
    print("liczę S34")
    w_S34 = scipy.integrate.dblquad(S34,XWellA-3,XWellA+3,lambda x:R-3, lambda x:R+3)
    print("liczę S41")
    w_S41 = scipy.integrate.dblquad(S41,XWellA-3,XWellA+3,lambda x:R-3, lambda x:R+3)
    print("liczę S42")
    w_S42 = scipy.integrate.dblquad(S42,XWellA-3,XWellA+3,lambda x:R-3, lambda x:R+3)
    print("liczę S43")
    w_S43 = scipy.integrate.dblquad(S43,XWellA-3,XWellA+3,lambda x:R-3, lambda x:R+3)
    print("liczę S44")
    w_S44 = scipy.integrate.dblquad(S44,XWellA-3,XWellA+3,lambda x:R-3, lambda x:R+3)

    print("liczę h11")
    w_h11 = scipy.integrate.dblquad(h11,XWellA-3,XWellA+3,lambda x:R-3, lambda x:R+3)
    print("liczę h12")
    w_h12 = scipy.integrate.dblquad(h12,XWellA-3,XWellA+3,lambda x:R-3, lambda x:R+3)
    print("liczę h13")
    w_h13 = scipy.integrate.dblquad(h13,XWellA-3,XWellA+3,lambda x:R-3, lambda x:R+3)
    print("liczę h14")
    w_h14 = scipy.integrate.dblquad(h14,XWellA-3,XWellA+3,lambda x:R-3, lambda x:R+3)
    print("liczę h21")
    w_h21 = scipy.integrate.dblquad(h21,XWellA-3,XWellA+3,lambda x:R-3, lambda x:R+3)
    print("liczę h22")
    w_h22 = scipy.integrate.dblquad(h22,XWellA-3,XWellA+3,lambda x:R-3, lambda x:R+3)
    print("liczę h23")
    w_h23 = scipy.integrate.dblquad(h23,XWellA-3,XWellA+3,lambda x:R-3, lambda x:R+3)
    print("liczę h24")
    w_h24 = scipy.integrate.dblquad(h24,XWellA-3,XWellA+3,lambda x:R-3, lambda x:R+3)
    print("liczę h31")
    w_h31 = scipy.integrate.dblquad(h31,XWellA-3,XWellA+3,lambda x:R-3, lambda x:R+3)
    print("liczę h32")
    w_h32 = scipy.integrate.dblquad(h32,XWellA-3,XWellA+3,lambda x:R-3, lambda x:R+3)
    print("liczę h33")
    w_h33 = scipy.integrate.dblquad(h33,XWellA-3,XWellA+3,lambda x:R-3, lambda x:R+3)
    print("liczę h34")
    w_h34 = scipy.integrate.dblquad(h34,XWellA-3,XWellA+3,lambda x:R-3, lambda x:R+3)
    print("liczę h41")
    w_h41 = scipy.integrate.dblquad(h41,XWellA-3,XWellA+3,lambda x:R-3, lambda x:R+3)
    print("liczę h42")
    w_h42 = scipy.integrate.dblquad(h42,XWellA-3,XWellA+3,lambda x:R-3, lambda x:R+3)
    print("liczę h43")
    w_h43 = scipy.integrate.dblquad(h43,XWellA-3,XWellA+3,lambda x:R-3, lambda x:R+3)
    print("liczę h44")
    w_h44 = scipy.integrate.dblquad(h44,XWellA-3,XWellA+3,lambda x:R-3, lambda x:R+3)

    matrixS = np.array([[w_S11[0],w_S12[0],w_S13[0],w_S14[0]],[w_S21[0],w_S22[0],w_S23[0],w_S24[0]],[w_S31[0],w_S32[0],w_S33[0],w_S34[0]],[w_S41[0],w_S42[0],w_S43[0],w_S44[0]]])
    #print(matrixS)
    matrixH = np.array([[w_h11[0],w_h12[0],w_h13[0],w_h14[0]],[w_h21[0],w_h22[0],w_h23[0],w_h24[0]],[w_h31[0],w_h32[0],w_h33[0],w_h34[0]],[w_h41[0],w_h42[0],w_h43[0],w_h44[0]]])
    #print(matrixH)


    #c_Vector = [C1,C2,C3,C4]
    E = []
    c_Vector = []

    E, c_Vector = scipy.linalg.eig(matrixH,matrixS)
    #print('energie')
    #print(E)
    #print('wektor C')
    #print(c_Vector)
    c_vals = []
    #sum = 0
    for i in range(0,4):
        c_vals.append(c_Vector[i][np.argmin(E)])
        if isinstance(c_vals[i], complex):
            c_vals[i]= np.absolute(c_vals[i])
        else:
            c_vals[i] = c_vals[i]**2
        #sum+=c_vals[i]**2
    #print(c_Vector[np.argmin(E)])
    #print(c_vals)
    #print(sum)

    print(c_Vector)
    x_axis.append(R/3)
    C1_values.append(c_vals[0])
    C2_values.append((c_vals[1]))
    C3_values.append((c_vals[2]))
    C4_values.append((c_vals[3]))
    E_values.append(min(E))
    print("C1")
    print(C1_values)
    print("C2")
    print(C2_values)
    print("C3")
    print(C3_values)
    print("C4")
    print(C4_values)
    print('E')
    print(E)
    print('Emin')
    print(min(E))


    #E[:] = []
    #c_Vector = []

plt.figure(1)
C1_scatter = plt.scatter(x_axis,C1_values,c='r',marker = 'o')
C2_scatter = plt.scatter(x_axis,C2_values,c='g',marker = ',')
C3_scatter = plt.scatter(x_axis,C3_values,c='b', marker = 'x')
C4_scatter = plt.scatter(x_axis,C4_values,c='y',marker = 'p')
plt.ylim(0,1)
plt.xlim(0.01,1.15)
plt.xlabel('d/D')
plt.ylabel('prawdopodobieństwo')
plt.legend((C1_scatter,C2_scatter,C3_scatter,C4_scatter),('|C1|^2','|C2|^2','|C3|^2','|C4|^2'))
plt.show()


with open('dane_metoda_wariacyjna2.csv', 'w') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerows(zip(x_axis,C1_values,C2_values,C3_values,C4_values,E_values))

    quit()




