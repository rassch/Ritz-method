from sympy import Symbol, Integral, Derivative, init_printing, symbols, solve, diff, exp,Abs


def function ():
    print("Hello world")
function()

b = 3.14
c = 0.85
d = 0.45
h_bar = 6.582119 * 10**-16 #eV*s
m_c_squared = 939*10**6 #eV


R = Symbol("R")
B = Symbol("B")
C = Symbol("C")
D = Symbol("D")
X1 = Symbol("x1")
X2 = Symbol("x2")

internal_potential = B*(exp((-2*(R-C))/D) - 2*exp((-(R-C))/D))
print(internal_potential)
internal_potential = internal_potential.subs(B,b)
internal_potential = internal_potential.subs(C,c)
internal_potential = internal_potential.subs(D,d)
internal_potential = internal_potential.subs(R, Abs(X1-X2))
print(internal_potential)

#wierzba ladnie opisal
