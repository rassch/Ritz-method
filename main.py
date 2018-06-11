#import sympy
from sympy import Symbol, integrate, Derivative, init_printing, symbols, solve, diff, exp,Abs, oo , Matrix

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
#alpha1 = Symbol("alpha1")#sqrt(mk/pi**2)
#alpha2 = Symbol("alpha2")#sqrt(mk/pi**2)
psi = Symbol("psi")
H11 = Symbol("H11")
H12 = Symbol("H12")
H13 = Symbol("H13")
H14 = Symbol("H14")

H21 = Symbol("H21")
H22 = Symbol("H22")
H23 = Symbol("H23")
H24 = Symbol("H24")

H31 = Symbol("H31")
H32 = Symbol("H32")
H33 = Symbol("H33")
H34 = Symbol("H34")

H41 = Symbol("H41")
H42 = Symbol("H42")
H43 = Symbol("H43")
H44 = Symbol("H44")

S11 = Symbol("S11")
S12 = Symbol("S12")
S13 = Symbol("S13")
S14 = Symbol("S14")

S21 = Symbol("S21")
S22 = Symbol("S22")
S23 = Symbol("S23")
S24 = Symbol("S24")

S31 = Symbol("S31")
S32 = Symbol("S32")
S33 = Symbol("S33")
S34 = Symbol("S34")

S41 = Symbol("S41")
S42 = Symbol("S42")
S43 = Symbol("S43")
S44 = Symbol("S44")


#narazie przyjmujemy ,że alfa1 alfa2 a1 a2 =1 , a polozenia studni jako 0 i 2

alpha1 = 1
alpha2 = 1
a1 = 1
a2 = 2
XWellA = 0
XWellB = 2
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
#print("funkcja psi:")
#print(test)
#deklaracja poszczegolnych psi
psi1 = U(X1,"a")*U(X2,"b")
psi2 = U(X2,"a")*U(X1,"b")
psi3 = U(X1,"a")*U(X1,"b")
psi4 = U(X2,"a")*U(X2,"b")

psi_tab = [psi1,psi2,psi3,psi4]

#print(psi_tab)

psi = test
psi = psi.subs(alpha1,1) #trzeba wyznaczyć k, żeby wstawić dobrą alfę
psi = psi.subs(alpha2,1) #inaczej nie policzy się całka
psi = psi.subs(x,1)
# iloczyn skalarny
psi_dot_product = integrate(psi1*psi2, (X1, -oo, oo), (X2, -oo, oo)) #test psi1 z psi2
print("iloczyn skalarny: ")
print (psi_dot_product)


#hamiltonian
def hamiltonian(formula):
    return double_derivative(formula,X1)+double_derivative(formula,X2)+internal_potential+V_ext

#iloczyn skalarny
def dot_product(i,j):
    return integrate(psi_tab[i]*psi_tab[j], (X1, -oo, oo), (X2, -oo, oo))
#deklaracja wektora stalych normujacych
C_VECTOR = Matrix([C1,C2,C3,C4])


#deklaracja maciezry S i H
#H_MATRIX = Matrix([[H11,H12,H13,H14],[H21,H22,H23,H24],[H31,H32,H33,H34],[H41,H42,H43,H44]])
#S_MATRIX = Matrix([[S11,S12,S13,S14],[S21,S22,S23,S24],[S31,S32,S33,S34],[S41,S42,S43,S44]])

H_MATRIX = Matrix([[1,1,1,1],[1,1,1,1],[1,1,1,1],[1,1,1,1]])#narazie sa puste
S_MATRIX = Matrix([[1,1,1,1],[1,1,1,1],[1,1,1,1],[1,1,1,1]])

#Obliczanie wspolczynnikow macierzy S

for i in (0,3,1):
    for j in (0,3,1):
        S_MATRIX[i,j]= S_MATRIX[i,j]*dot_product(i,j)

print(S_MATRIX)
