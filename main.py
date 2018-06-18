from statistics import variance

import sympy
#from scipy.stats import alpha_gen
from mpmath import psi
from sympy import Symbol, integrate, Derivative, init_printing, symbols, solve, diff, exp,Abs, oo , Matrix,ones, simplify
from scipy.optimize import minimize
import scipy
import numpy as np
from scipy import  integrate

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

XWellA = 1

XWellB = 3

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
    return 3.14* (np.exp((-2 * (np.abs(x-y) - 0.85)) / 0.45) - 2 * np.exp((-np.abs(x-y) - 0.85)) / 0.45)
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
def psiAB(x,y):
    return ((alpha1**1/2)/pi**1/4)*np.exp(-alpha1**2*(x-XWellA)**2/2)*((alpha1**1/2)/pi**1/4)*np.exp(-alpha1**2*(y-XWellB)**2/2)
def psiBA(x,y):
    return ((alpha1 ** 1 / 2) / pi ** 1 / 4) * np.exp(-alpha1 ** 2 * (y - XWellA) ** 2 / 2) * (
                (alpha1 ** 1 / 2) / pi ** 1 / 4) * np.exp(-alpha1 ** 2 * (x - XWellB) ** 2 / 2)

def psiAA(x,y):
    return ((alpha1 ** 1 / 2) / pi ** 1 / 4) * np.exp(-alpha1 ** 2 * (x - XWellA) ** 2 / 2) * (
            (alpha1 ** 1 / 2) / pi ** 1 / 4) * np.exp(-alpha1 ** 2 * (y - XWellA) ** 2 / 2)

def psiBB(x,y):
    return ((alpha1 ** 1 / 2) / pi ** 1 / 4) * np.exp(-alpha1 ** 2 * (x - XWellB) ** 2 / 2) * (
                (alpha1 ** 1 / 2) / pi ** 1 / 4) * np.exp(-alpha1 ** 2 * (y - XWellB) ** 2 / 2)

def hamiltonianAA(x,y):
    return ((-h_c**2)/(2*m_c_squared))*diff(psiAA(x,y),x,x)+((-h_c**2)/(2*m_c_squared))*diff(psiAA,y,y)+(internal_potential(x,y)+ext_potential(x)+ext_potential(y))

def hamiltonianAB(x, y):
    return double_derivative(psiAB(x, y), x) + double_derivative(psiAB(x, y), y) + (
                    internal_potential(x, y) + ext_potential(x) + ext_potential(y))

def hamiltonianBA(x, y):
    return double_derivative(psiBA(x, y), x) + double_derivative(psiBA(x, y), y) + internal_potential(x, y) + ext_potential(x) + ext_potential(y)

def hamiltonianBB(x, y):
    return double_derivative(psiBB(x, y), x) + double_derivative(psiBB(x, y), y) + (internal_potential(x, y) + ext_potential(x) + ext_potential(y))


# Wektor C
#print('jeszcze raz internal')
#print(internal_potential)
#print('ext1')
#print(ext_potential(X1))
#print('druga pochodna')
#print(double_derivative(psi_tab[0],X1))

c_Vector = [C1,C2,C3,C4]
#print('hamiltionian')
#print(hamiltonian(psi_tab[0],psi_tab[1],X1,X2))
caleczka = scipy.integrate.dblquad(hamiltonianAA,0,1,lambda x:0, lambda x:1)
print(caleczka)

'''
H_Matrix = []
for i in range(0,3):
    new = []
    for j in range(0,3):
        print('cos sie robi')
        new.append(scipy.integrate.dblquad(lambda X1,X2: hamiltonian(psi_tab[i],psi_tab[j],X1,X2),(X1,XWellA-3,XWellA+3),(X2,XWellB-3,XWellB+3)))
    H_Matrix.append(new)
print(H_Matrix)

S_Matrix = []
for i in range(0,3):
    new2 = []
    for j in range(0,3):
        print('cos sie robi')
        new2.append(scipy.integrate.dblquad(lambda X1,X2: hamiltonian(psi_tab[i],psi_tab[j],X1,X2),(X1,XWellA-3,XWellA+3),(X2,XWellB-3,XWellB+3)))
    C_Matrix.append(new2)
print(C_Matrix)


'''
'''
print('Test')
print(U(X1,"a")*U(X2,"b"))
print()
test2 = integrate(psi_tab[0],(X1,XWellA-3,XWellA+3))*integrate(psi_tab[0],(X2,XWellB-3,XWellB+3))
print("wynik testu")
print(simplify(test2))
#h1 = integrate(psi_tab[0]*hamiltonian(psi_tab[0],X1,X2),X1)*integrate(psi_tab[0]*hamiltonian(psi_tab[0],X1,X2),X2)
#print(h1)
print(hamiltonian(psi_tab[0],psi_tab[1],X1,X2))

calka_ham =  integrate(hamiltonian(psi_tab[0],psi_tab[1],X1,X2),X1)*integrate(hamiltonian(psi_tab[0],psi_tab[1],X1,X2),X2)
print(calka_ham)











psi = test

print(test)

psi = psi.subs(x,1)

psi_squared = psi*psi

do_h_psi = hamiltonian(psi)

h_psi = psi * do_h_psi

print("podstawianie: ")

suma = []

for i in range(0,6,1):

    suma.append(0)

całka_h = 0

for i in range(0,6,1): #musi być int

    for j in range(0, 6, 1):

        suma[j] = suma[j] +((1/2)*(h_psi.subs([(X1,i),(X2,j)])+h_psi.subs([(X1,i+1),(X2,j)])))



print("całka z h_psi")

for i in range(0,5,1): #musi być int

    całka_h = całka_h + ((1/2)*(suma[i]+suma[i+1]))

#function_to_optimize = integrate(h_psi,(X1,-3/sympy.sqrt(alpha1),3/sympy.sqrt(alpha1)),(X2,-3/sympy.sqrt(alpha2),3/sympy.sqrt(alpha2)))/integrate(psi_squared,(X1,-3/sympy.sqrt(alpha1),3/sympy.sqrt(alpha1)),(X2,-3/sympy.sqrt(alpha2),3/sympy.sqrt(alpha2)))#FUNKCJA INTEGRATE NIE WSPÓŁPRACUJE



#print(function_to_optimize)

print(całka_h)





całka = 0

for i in range(0,6,1): #musi być int

    for j in range(0, 6, 1):

        suma[j] = suma[j] +((1/2)*(psi_squared.subs([(X1,i),(X2,j)])+psi_squared.subs([(X1,i+1),(X2,j)])))



print("całka z bez h:")

for i in range(0,5,1): #musi być int

    całka = całka + ((1/2)*(suma[i]+suma[i+1]))



print(całka)



def f(alfa):

    return całka_h.subs([(alpha1,alfa[0]),(alpha2,alfa[1])])/całka.subs([(alpha1,alfa[0]),(alpha2,alfa[1])])



# Stałe C1... chyba juz niepotrzebne...



#nie potrafię podstawić do funkcji fmin funkcji z dwoma zmiennymi
x0 = [3,2]
print("optymalizacja")
sol = minimize(f,x0)
print(sol)






#należy znaleźć minimum

'''