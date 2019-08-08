from __future__ import division
from sympy import *

init_printing()
x = symbols('x') # definir simbolos

def f(x): return exp(x) # funcao

print(f(0)) # f(0)
print(diff(f(x), x).subs(x, 0)) # f'(0)