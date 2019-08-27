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
arrs.append(1) # caso base da EDO

qt_it = math.ceil(n/h) # quantidade de iterações necessarias

t = 0 # t vai de 0 até n (0,n)
arrt.append(t) # iniciando com o primeiro passo

ordem = 2
# ordem = 2 (provada em sala)

consts = [
    [1],
    [1/2, 1/2],
    [-1/12, 2/3, 5/12],
    [1/24, -5/24, 19/24, 3/8],
    [-19/720, 53/360, -11/30, 323/360, 251/720],
    [3/160, -173/1440, 241/720, -133/240, 1427/1440, 95/288],
    [-863/60480, 263/2520, -6737/20160, 586/945, -15487/20160, 2713/2520, 19087/60480],
    [275/24192, -11351/120960, 1537/4480, -88547/120960, 123133/120960, -4511/4480, 139849/120960, 5257/17280]        
]

for i in range(ordem-1, qt_it+1):

    arrk = []
    aux = 1
    brincadeira = t
    for j in range(0, ordem-1):
        arrk.append(f(brincadeira, arry[(i-ordem)+aux]))
        aux = aux + 1
        brincadeira = brincadeira + h

    arrk.append(f(t, arry[i-1] + arrk[len(arrk) - 1] * h))

    sum = 0

    for j in range(0, len(arrk)):
        sum += arrk[j] * h * consts[ordem-1][j]

    arry.append(arry[i-1] + sum) # metodo de adam multon
    arrs.append(solve(t)) # solucao da EDO
    t = t + h # incrementando o passo
    arrt.append(t) # colocando o passo no array

# plotagem do grafico 
#graphic.title('Adam Multon Method')
#graphic.xlabel("Passos")
#graphic.ylabel("f(t, y)")

#graphic.plot(arrt, arry, color = 'red') # cor vermelha pra o grafico da sol. de adam multon
#graphic.plot(arrt, arrs, color = 'yellow') #cor amarela para o grafico da sol. da EDO

print("Result = {}".format(arry[qt_it]))
print("Exact = {}".format(solve(n)))
#graphic.show()
