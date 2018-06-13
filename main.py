#import sympy
from sympy import Symbol, integrate, Derivative, init_printing, symbols, solve, diff, exp,Abs, oo , Matrix,ones
from scipy import fmin

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
c = 299792458 #m/s
pi = 3.14159265359
A1 = -6.14
A2 = -6.14

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
alpha1 = Symbol("alpha1")#sqrt(mk/pi**2)
alpha2 = Symbol("alpha2")#sqrt(mk/pi**2)
psi = Symbol("psi")


#narazie przyjmujemy ,że alfa1 alfa2 a1 a2 =1 , a polozenia studni jako 0 i 2


a1 = 1
a2 = 2
XWellA = 1
XWellB = 3
#######################################################

def U(variable,potential_well): # U w zależności od studni i atomu
    if potential_well == "a":
        return ((alpha1**(1/2))/(pi**(1/4)))*exp((-1/2)*alpha1**2*(variable-XWellA)**2)
    else:
        return ((alpha2 ** (1 / 2)) / (pi ** (1 / 4))) * exp((-1 / 2) * alpha2 ** 2 * (variable - XWellB)**2)


def double_derivative(formula , variable): #do hamiltonianu
    outcome = ((-h_bar*c_squared)/(2*m_c_squared))*diff(formula,variable,variable)
    return outcome





internal_potential = B*(exp((-2*(R-C))/D) - 2*exp((-(R-C))/D)) #potencjał wewnętrzny
print(internal_potential)
internal_potential = internal_potential.subs(B,b) #podstawiam podane wartości
internal_potential = internal_potential.subs(C,c)
internal_potential = internal_potential.subs(D,d)
internal_potential = internal_potential.subs(R, Abs(X1-X2))
print(internal_potential)


V_ext = A1*exp(-(a1*(x-X1)**2)/2)+A2*exp(-(a2*(x-X2)**2)/2)
V_ext = V_ext.subs(A1,-6.14)
V_ext = V_ext.subs(A2,-6.14)
V_ext = V_ext.subs(a1,1)
V_ext = V_ext.subs(a2,1)
#od razu suma V_ext(x1) i V_ext(x2)
#nie wiem co podstawić za a1,a2,A1,A2



#wymnożone wszytkie U -- funkcja psi
test = C1*U(X1,"a")*U(X2,"b")+C2*U(X2,"a")*U(X1,"b")+C3*U(X1,"a")*U(X1,"b")+C4*U(X2,"a")*U(X2,"b") #proszę sprawdzić czy dobrze




#hamiltonian
def hamiltonian(formula):
    return double_derivative(formula,X1)+double_derivative(formula,X2)+(internal_potential+V_ext)*formula



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
#function_to_optimize = integrate(h_psi,(X1,-oo,oo),(X2,-oo,oo))/integrate(psi_squared,(X1,-oo,oo),(X2,-oo,oo))#FUNKCJA INTEGRATE NIE WSPÓŁPRACUJE


print(całka_h)


całka = 0
for i in range(0,6,1): #musi być int
    for j in range(0, 6, 1):
        suma[j] = suma[j] +((1/2)*(psi_squared.subs([(X1,i),(X2,j)])+psi_squared.subs([(X1,i+1),(X2,j)])))

print("całka z bez h:")
for i in range(0,5,1): #musi być int
    całka = całka + ((1/2)*(suma[i]+suma[i+1]))

print(całka)

def f(alfa1,alfa2):
 return całka_h.subs([(alpha1,alfa1),(alpha2,alfa2)])/całka.subs([(alpha1,alfa1),(alpha2,alfa2)])


# Stałe C1... chyba juz niepotrzebne...

#nie potrafię podstawić do funkcji fmin funkcji z dwoma zmiennymi




#należy znaleźć minimum





