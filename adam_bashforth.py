import math
#import matplotlib.pyplot as graphic

arry = [] # vetor com os resultados
arrt = [] # vetor com os passos
arrs = [] # vetor com os resultados exatos (solução da edo)

def f(t, y): return (1 - t + (4 * y)) # y(t, y)
def solve(t): return ((1/4)*t - (3/16) + (19/16)*math.exp(4*t)) # solucao da EDO

n = int(input()) # quantidade de pontos
h = float(input()) # valor do passo

arry.append(1) # caso base
arry.append(1.00501) # caso base 2

arrs.append(1) # caso base da EDO

qt_it = math.ceil(n/h) # quantidade de iterações necessarias

t = 0 # t vai de 0 até n (0,n)
arrt.append(t) # iniciando com o primeiro passo

ordem = 2
# ordem = 2 (provada em sala)

consts = [
    [1],
    [-1/2, 3/2],
    [5/12, -4/3, 23/12],
    [-3/8, 37/24, -59/24, 55/24],
    [251/720, -637/360, 109/30, -1387/360, 1901/720],
    [-95/288, 959/480, -3649/720, 4991/720, -2641/480, 4277/1440],
    [19087/60480, -5603/2520, 135713/20160, -10754/945, 235183/20160, -18637/2520, 198721/60480],
    [-5257/17280, 32863/13440, -115747/13440, 2102243/120960, -296053/13440, 242653/13440, -1152169/120960, 16083/4480]        
]

for i in range(ordem-1, qt_it+1):

    arrk = []
    aux = 1
    aux_t = t
    for j in range(0, ordem):
        arrk.append(f(aux_t, arry[(i-ordem)+aux]))
        aux = aux + 1 #pegar prx termo do arry
        aux_t = aux_t + h # proximo t

    sum = 0

    for j in range(0, len(arrk)):
        sum += arrk[j] * h * consts[ordem-1][j]

    arry.append(arry[len(arry)-1] + sum) # metodo de adam bashfort
    arrs.append(solve(t)) # solucao da EDO
    t = t + h # incrementando o passo
    arrt.append(t) # colocando o passo no array

# plotagem do grafico 
#graphic.title('Adam Multon Method')
#graphic.xlabel("Passos")
#graphic.ylabel("f(t, y)")

#graphic.plot(arrt, arry, color = 'red') # cor vermelha pra o grafico da sol. de adam multon
#graphic.plot(arrt, arrs, color = 'yellow') #cor amarela para o grafico da sol. da EDO

print("Result = {}".format(arry[len(arry)-1]))
print("Exact = {}".format(solve(n)))
#graphic.show()
