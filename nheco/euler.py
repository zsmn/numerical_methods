u = [] # vetor com os valores
def f(x): return x # funcao f(x)

num_it = int(input()) # iteracao que se deseja 
tam = float(input()) # tamanho do intervalo (quanto menor, maior a precisao)

u.append(1) # caso base

qt_iteration = int((num_it/tam)+1) # quantia de iteracoes necessarias para chegar no resultado

for i in range(0, qt_iteration):
    u.append(u[i] + tam*f(u[i])) # metodo de euler

print("Result = {}".format(u[qt_iteration-1])) # resultado final