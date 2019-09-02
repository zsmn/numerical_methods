import math
import sympy as sym
from sympy.parsing.sympy_parser import parse_expr

pts_y = [] # resultados
pts_t = [] # passos
y, t = sym.symbols('y t')

consts_multon = [
    [1],
    [1/2, 1/2],
    [-1/12, 2/3, 5/12],
    [1/24, -5/24, 19/24, 3/8],
    [-19/720, 53/360, -11/30, 323/360, 251/720],
    [3/160, -173/1440, 241/720, -133/240, 1427/1440, 95/288],
    [-863/60480, 263/2520, -6737/20160, 586/945, -15487/20160, 2713/2520, 19087/60480],
    [275/24192, -11351/120960, 1537/4480, -88547/120960, 123133/120960, -4511/4480, 139849/120960, 5257/17280]        
]

consts_bashforth = [
    [1],
    [-1/2, 3/2],
    [5/12, -4/3, 23/12],
    [-3/8, 37/24, -59/24, 55/24],
    [251/720, -637/360, 109/30, -1387/360, 1901/720],
    [-95/288, 959/480, -3649/720, 4991/720, -2641/480, 4277/1440],
    [19087/60480, -5603/2520, 135713/20160, -10754/945, 235183/20160, -18637/2520, 198721/60480],
    [-5257/17280, 32863/13440, -115747/13440, 2102243/120960, -296053/13440, 242653/13440, -1152169/120960, 16083/4480]        
]

def euler(funct, y0, t0, qt_it, h):
    f = parse_expr(funct)
    arrt = []
    arry = []
    
    arrt.append(float(t0))
    arry.append(float(y0))

    for i in range(1, int(qt_it)+1):
        arry.append(arry[i-1] + float(h) * f.subs(t, arrt[i-1]).subs(y, arry[i-1]))
        arrt.append(arrt[i-1] + float(h))
    
    return arry, arrt

def reverse_euler(funct, y0, t0, qt_it, h):
    f = parse_expr(funct)
    arrt = []
    arry = []
    
    arrt.append(float(t0))
    arry.append(float(y0))

    for i in range(0, int(qt_it)):
        yaux = arry[i] + float(h) * f.subs(t, arrt[i]).subs(y, arry[i])
        taux = arrt[i] + float(h)

        arry.append(arry[i] + float(h) * f.subs(t, taux).subs(y, yaux))
        arrt.append(arrt[i] + float(h))

    return arry, arrt

def aprimorated_euler(funct, y0, t0, qt_it, h):
    f = parse_expr(funct)
    arrt = []
    arry = []
    
    arrt.append(float(t0))
    arry.append(float(y0))

    for i in range(1, int(qt_it)+1):
        k1 = f.subs(t, arrt[i-1]).subs(y, arry[i-1])
        k2 = f.subs(t, arrt[i-1]+float(h)).subs(y, arry[i-1] + float(h) * k1)
        arry.append(arry[i-1] + (float(h)/2) * (k1 + k2))
        arrt.append(arrt[i-1] + float(h))
    
    return arry, arrt

def runge_kutta(funct, y0, t0, qt_it, h):
    f = parse_expr(funct)
    arrt = []
    arry = []
    
    arrt.append(float(t0))
    arry.append(float(y0))

    for i in range(1, int(qt_it)+1):
        k1 = f.subs(t, arrt[i-1]).subs(y, arry[i-1])
        k2 = f.subs(t, arrt[i-1] + 0.5 * float(h)).subs(y, arry[i-1] + 0.5 * float(h) * k1)
        k3 = f.subs(t, arrt[i-1] + 0.5 * float(h)).subs(y, arry[i-1] + 0.5 * float(h) * k2)
        k4 = f.subs(t, arrt[i-1] + float(h)).subs(y, arry[i-1] + float(h) * k3)
        arry.append(arry[i-1] + (float(h)/6) * (k1 + (2*k2) + (2*k3) + k4))
        arrt.append(arrt[i-1] + float(h))

    return arry, arrt

def adam_multon(funct, vety, vett, qt_it, h, ordem):
    f = parse_expr(funct)
    arrt = vett
    arry = vety
    for i in range(int(ordem)-1, int(qt_it)+(int(ordem)-1)):
        arrk = []
        aux = 1
        aux_t = float(arrt[i-1])

        for j in range(0, int(ordem)-1):
            arrk.append(f.subs(t, aux_t).subs(y, float(arry[(i-int(ordem))+aux])))
            aux = aux + 1
            aux_t = aux_t + float(h)

        arrk.append(f.subs(t, arrt[i-1]).subs(y, float(arry[i-1]) + arrk[len(arrk) - 1] * float(h)))
        sum = 0

        for j in range(0, len(arrk)):
            sum += arrk[j] * float(h) * consts_multon[int(ordem)-1][j]

        arry.append(float(arry[i-1]) + sum)
        arrt.append(float(arrt[i-1]) + float(h))

    return arry, arrt

def adam_bashforth(funct, vety, vett, qt_it, h, ordem):
    f = parse_expr(funct)
    arrt = vett
    arry = vety
    for i in range(int(ordem)-1, int(qt_it)+(int(ordem)-1)):
        arrk = []
        aux = 1
        aux_t = float(arrt[i-1])

        for j in range(0, int(ordem)):
            arrk.append(f.subs(t, aux_t).subs(y, float(arry[(i-int(ordem))+aux])))
            aux = aux + 1
            aux_t = aux_t + float(h)

        sum = 0

        for j in range(0, len(arrk)):
            sum += arrk[j] * float(h) * consts_bashforth[int(ordem)-1][j]

        arry.append(float(arry[len(arry)-1]) + sum) 
        arrt.append(float(arrt[i-1]) + float(h))
    
    return arry, arrt

def main():
    arq = open('input.txt', 'r')

    for linha in arq:
        param = linha.split()
        if(param[0] == "euler"):
            pts_y, pts_t = euler(param[5], param[2], param[1], param[4], param[3])
            print(pts_y[int(param[4])])
        elif(param[0] == "euler_inverso"):
            pts_y, pts_t = reverse_euler(param[5], param[2], param[1], param[4], param[3])
            print(pts_y[int(param[4])])
        elif(param[0] == "euler_aprimorado"):
            pts_y, pts_t = aprimorated_euler(param[5], param[2], param[1], param[4], param[3])
            print(pts_y[int(param[4])])
        elif(param[0] == "runge_kutta"):
            pts_y, pts_t = runge_kutta(param[5], param[2], param[1], param[4], param[3])
            print(pts_y[int(param[4])])
        elif("adam_multon" in param[0]):
            if("by_euler_inverso" in param[0]):
                #print("eh com euler inverso")
                pts_y, pts_t = reverse_euler(param[len(param) - 2], param[2], param[1], int(param[len(param) - 1]) - 2, param[3])
                pts_y, pts_t = adam_multon(param[len(param) - 2], pts_y, pts_t, param[len(param) - 3], param[len(param) - 4], param[len(param) - 1])
                print(pts_y[len(pts_y) - 1])
            elif("by_euler_aprimorado" in param[0]):
                #print("eh com euler aprimorado")
                pts_y, pts_t = aprimorated_euler(param[len(param) - 2], param[2], param[1], int(param[len(param) - 1]) - 2, param[3])
                pts_y, pts_t = adam_multon(param[len(param) - 2], pts_y, pts_t, param[len(param) - 3], param[len(param) - 4], param[len(param) - 1])
                print(pts_y[len(pts_y) - 1])
            elif("by_runge_kutta" in param[0]):
                #print("eh com runge kutta")
                pts_y, pts_t = runge_kutta(param[len(param) - 2], param[2], param[1], int(param[len(param) - 1]) - 2, param[3])
                pts_y, pts_t = adam_multon(param[len(param) - 2], pts_y, pts_t, param[len(param) - 3], param[len(param) - 4], param[len(param) - 1])
                print(pts_y[len(pts_y) - 1])
            elif("by_euler" in param[0]):
                #print("eh com euler normal")
                pts_y, pts_t = euler(param[len(param) - 2], param[2], param[1], int(param[len(param) - 1]) - 2, param[3])
                pts_y, pts_t = adam_multon(param[len(param) - 2], pts_y, pts_t, param[len(param) - 3], param[len(param) - 4], param[len(param) - 1])
                print(pts_y[len(pts_y) - 1])
            else:
                #print("eh so o multon")
                pts_y = []
                pts_t = []
                h_aux = 0
                for i in range(int(param[len(param) - 1]) - 1):
                    pts_y.append(float(param[i+2]))
                    pts_t.append(float(param[1]) + h_aux)
                    #print(i, " ", pts_y[i], " ", pts_t[i])
                    h_aux = h_aux + float(param[len(param) - 4])
                
                pts_y, pts_t = adam_multon(param[len(param) - 2], pts_y, pts_t, param[len(param) - 3], param[len(param) - 4], param[len(param) - 1])
                print(pts_y[len(pts_y) - 1])
        elif("adam_bashforth" in param[0]):
            if("by_euler_inverso" in param[0]):
                #print("bash eh com euler inverso")
                pts_y, pts_t = reverse_euler(param[len(param) - 2], param[2], param[1], int(param[len(param) - 1]) - 1, param[3])
                pts_y, pts_t = adam_bashforth(param[len(param) - 2], pts_y, pts_t, param[len(param) - 3], param[len(param) - 4], param[len(param) - 1])
                print(pts_y[len(pts_y) - 1])
            elif("by_euler_aprimorado" in param[0]):
                #print("bash eh com euler aprimorado")
                pts_y, pts_t = aprimorated_euler(param[len(param) - 2], param[2], param[1], int(param[len(param) - 1]) - 1, param[3])
                pts_y, pts_t = adam_bashforth(param[len(param) - 2], pts_y, pts_t, param[len(param) - 3], param[len(param) - 4], param[len(param) - 1])
                print(pts_y[len(pts_y) - 1])
            elif("by_runge_kutta" in param[0]):
                #print("bash eh com runge kutta")
                pts_y, pts_t = runge_kutta(param[len(param) - 2], param[2], param[1], int(param[len(param) - 1]) - 1, param[3])
                pts_y, pts_t = adam_bashforth(param[len(param) - 2], pts_y, pts_t, param[len(param) - 3], param[len(param) - 4], param[len(param) - 1])
                print(pts_y[len(pts_y) - 1])
            elif("by_euler" in param[0]):
                #print("bash eh com euler normal")
                pts_y, pts_t = euler(param[len(param) - 2], param[2], param[1], int(param[len(param) - 1]) - 1, param[3])
                pts_y, pts_t = adam_bashforth(param[len(param) - 2], pts_y, pts_t, param[len(param) - 3], param[len(param) - 4], param[len(param) - 1])
                print(pts_y[len(pts_y) - 1])
            else:
                #print("eh so o bashforth")
                pts_y = []
                pts_t = []
                h_aux = 0
                for i in range(int(param[len(param) - 1])):
                    pts_y.append(float(param[i+2]))
                    pts_t.append(float(param[1]) + h_aux)
                    #print(i, " ", pts_y[i], " ", pts_t[i])
                    h_aux = h_aux + float(param[len(param) - 4])
                
                pts_y, pts_t = adam_bashforth(param[len(param) - 2], pts_y, pts_t, param[len(param) - 3], param[len(param) - 4], param[len(param) - 1])
                print(pts_y[len(pts_y) - 1])
main()