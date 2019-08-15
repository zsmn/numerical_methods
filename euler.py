import math

u = [] # vetor com os valores 
def f(t, y): return (1 - t + (4 * y)) # y(t, y)

n = int(input()) # quantidade de pontos
h = float(input()) # valor do passo

u.append(1) # caso base
qt_it = math.ceil(n/h) # quantidade de iterações necessarias

t = 0 # t vai de 0 até n (0,n)

for i in range(1, qt_it+1):
    u.append(u[i-1] + h * f(t, u[i-1])) # metodo de euler
    t = t + h

print("Result = {}".format(u[qt_it]))