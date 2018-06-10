from sympy import Symbol, integrate, Derivative, init_printing, symbols, solve, diff, exp,Abs, oo


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

R = Symbol("R")
B = Symbol("B")
C = Symbol("C")
D = Symbol("D")
x = Symbol("x")
X1 = Symbol("x1")
X2 = Symbol("x2")
a1 = Symbol("a1")
A1 = Symbol("A1")
a2 = Symbol("a2")
A2 = Symbol ("A2")
k1 = Symbol("k1") #na razie przyjmujemy że to sqrt(2mE)/h
k2 = Symbol("k2")
E = Symbol("E")
alpha1 = Symbol("alpha1")#sqrt(mk/pi**2)
alpha2 = Symbol("alpha2")#sqrt(mk/pi**2)
psi = Symbol("psi")

def U(variable,potential_well): # U w zależności od studni i atomu
    if potential_well == "a":
        return ((alpha1**(1/2))/(pi**(1/4)))*exp((-1/2)*alpha1**2*(x-variable))
    else:
        return ((alpha2 ** (1 / 2)) / (pi ** (1 / 4))) * exp((-1 / 2) * alpha2 ** 2 * (x - variable))


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
#od razu suma V_ext(x1) i V_ext(x2)
#nie wiem co podstawić za a1,a2,A1,A2



#wymnożone wszytkie U -- funkcja psi
test = U(X1,"a")*U(X2,"b")+U(X2,"1")*U(X1,"b")+U(X1,"a")*U(X1,b)+U(X2,"a")*U(X2,"b") #proszę sprawdzić czy dobrze
print("funkcja psi:")
print(test)
psi = test
psi = psi.subs(alpha1,1) #trzeba wyznaczyć k, żeby wstawić dobrą alfę
psi = psi.subs(alpha2,1) #inaczej nie policzy się całka
psi = psi.subs(x,1)
# iloczyn skalarny
psi_dot_product = integrate(psi*psi, (X1, 0.5, 3), (X2, 0.5, 3)) #granice od 0.5 do 3 zgodnie z danymi
print("iloczyn skalarny: ")
print (psi_dot_product)


#hamiltonian
def hamiltonian(formula):
    return double_derivative(formula,X1)+double_derivative(formula,X2)+internal_potential+V_ext



#wierzba ladnie opisal
