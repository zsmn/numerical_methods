import math
import matplotlib.pyplot as graphic

arry = [] # vetor com os resultados
arrt = [] # vetor com os passos
arrs = [] # vetor com os resultados exatos (solução da edo)

def f(t, y): return (1 - t + (4 * y)) # y(t, y)
def solve(t): return ((1/4)*t - (3/16) + (19/16)*math.exp(4*t)) # solucao da EDO

n = int(input()) # quantidade de pontos
h = float(input()) # valor do passo

arry.append(1) # caso base
arrs.append(1) # caso base da EDO

qt_it = math.ceil(n/h) # quantidade de iterações necessarias

t = 0 # t vai de 0 até n (0,n)
arrt.append(t) # iniciando com o primeiro passo

for i in range(1, qt_it+1):
    arry.append(arry[i-1] + h * f(t, arry[i-1])) # metodo de euler
    arrs.append(solve(t)) # solucao da EDO
    t = t + h # incrementando o passo
    arrt.append(t) # colocando o passo no array

# plotagem do grafico 
graphic.title('Euler Method')
graphic.xlabel("Passos")
graphic.ylabel("f(t, y)")

graphic.plot(arrt, arry, color = 'red') # cor vermelha pra o grafico da sol. de euler
graphic.plot(arrt, arrs, color = 'yellow') #cor amarela para o grafico da sol. da EDO

print("Result = {}".format(arry[qt_it]))
graphic.show()