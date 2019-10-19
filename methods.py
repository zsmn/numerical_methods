import math
import sympy as sym
import matplotlib.pyplot as graphic
from sympy.parsing.sympy_parser import parse_expr

pts_y = [] # resultados
pts_t = [] # passos
y, t = sym.symbols('y t')

consts_bashforth = [
    [1.0],
	[3.0/2.0,-1.0/2.0],
	[23.0/12.0,-4.0/3.0,5.0/12.0],
	[55.0/24.0,-59.0/24.0,37.0/24.0,-3.0/8.0],
	[1901.0/720.0,-1387.0/360.0,109.0/30.0,-637.0/360.0,251.0/720.0],
	[4277.0/1440.0,-2641.0/480.0,4991.0/720.0,-3649.0/720.0,959.0/480.0,-95.0/288.0],
	[198721.0/60480.0,-18637.0/2520.0,235183.0/20160.0,-10754.0/945.0,135713.0/20160.0,-5603.0/2520.0,19087.0/60480.0],
	[16083.0/4480.0,-1152169.0/120960.0,242653.0/13440.0,-296053.0/13440.0,2102243.0/120960.0,-115747.0/13440.0,32863.0/13440.0,-5257.0/17280.0]
]

consts_multon = [
	[1.0],
	[1.0/2.0,1.0/2.0],
	[5.0/12.0,2.0/3.0,-1.0/12.0],
	[3.0/8.0,19.0/24.0,-5.0/24.0,1.0/24.0],
	[251.0/720.0,323.0/360.0,-11.0/30.0,53.0/360.0,-19.0/720.0],
	[95.0/288.0,1427.0/1440.0,-133.0/240.0,241.0/720.0,-173.0/1440.0,3.0/160.0],
	[19087.0/60480.0,2713.0/2520.0,-15487.0/20160.0,586.0/945.0,-6737.0/20160.0,263.0/2520.0,-863.0/60480.0],
	[5257.0/17280.0,139849.0/120960.0,-4511.0/4480.0,123133.0/120960.0,-88547.0/120960.0,1537.0/4480.0,-11351.0/120960.0,275.0/24192.0]
]

def plotGraphic(method, arrt, arry, color):
    graphic.title(method)
    graphic.xlabel("Passos")
    graphic.ylabel("f(t, y)")

    graphic.plot(arrt, arry, color)

    graphic.show()

    return

def euler(funct, y0, t0, qt_it, h):
    print('Euler')
    print('y({}) = {}'.format(t0, y0))
    print('it = {}'.format(qt_it))
    print('h = {}'.format(h))

    f = parse_expr(funct)
    arrt = []
    arry = []

    arrt.append(float(t0))
    arry.append(float(y0))

    for i in range(1, int(qt_it)+2):
        print(i-1, ' ', arry[i-1])
        arry.append(arry[i-1] + float(h) * f.subs(t, arrt[i-1]).subs(y, arry[i-1]))
        arrt.append(arrt[i-1] + float(h))

    return arry, arrt

def reverse_euler(funct, y0, t0, qt_it, h):
    print('Reverse Euler')
    print('y({}) = {}'.format(t0, y0))
    print('it = {}'.format(qt_it))
    print('h = {}'.format(h))

    f = parse_expr(funct)
    arrt = []
    arry = []
    
    arrt.append(float(t0))
    arry.append(float(y0))

    for i in range(0, int(qt_it) + 1):
        print(i, ' ', arry[i])
        yaux = arry[i] + float(h) * f.subs(t, arrt[i]).subs(y, arry[i])
        taux = arrt[i] + float(h)

        arry.append(arry[i] + float(h) * f.subs(t, taux).subs(y, yaux))
        arrt.append(arrt[i] + float(h))

    return arry, arrt

def aprimorated_euler(funct, y0, t0, qt_it, h):
    print('Aprimorated Euler')
    print('y({}) = {}'.format(t0, y0))
    print('it = {}'.format(qt_it))
    print('h = {}'.format(h))
    
    f = parse_expr(funct)
    arrt = []
    arry = []
    
    arrt.append(float(t0))
    arry.append(float(y0))

    for i in range(1, int(qt_it)+2):
        print(i - 1, ' ', arry[i])
        k1 = f.subs(t, arrt[i-1]).subs(y, arry[i-1])
        k2 = f.subs(t, arrt[i-1]+float(h)).subs(y, arry[i-1] + float(h) * k1)
        arry.append(arry[i-1] + (float(h)/2) * (k1 + k2))
        arrt.append(arrt[i-1] + float(h))
    
    return arry, arrt

def runge_kutta(funct, y0, t0, qt_it, h):
    print('Runge Kutta')
    print('y({}) = {}'.format(t0, y0))
    print('it = {}'.format(qt_it))
    print('h = {}'.format(h))

    f = parse_expr(funct)
    arrt = []
    arry = []
    
    arrt.append(float(t0))
    arry.append(float(y0))

    for i in range(1, int(qt_it)+2):
        print(i - 1, ' ', arry[i-1])
        k1 = f.subs(t, arrt[i-1]).subs(y, arry[i-1])
        k2 = f.subs(t, arrt[i-1] + 0.5 * float(h)).subs(y, arry[i-1] + 0.5 * float(h) * k1)
        k3 = f.subs(t, arrt[i-1] + 0.5 * float(h)).subs(y, arry[i-1] + 0.5 * float(h) * k2)
        k4 = f.subs(t, arrt[i-1] + float(h)).subs(y, arry[i-1] + float(h) * k3)
        arry.append(arry[i-1] + (float(h)/6) * (k1 + (2*k2) + (2*k3) + k4))
        arrt.append(arrt[i-1] + float(h))

    return arry, arrt

def bashfort_dynamic_p(funct, vety, tf, h, ordem):
    aux = 0.0
    yf = vety[len(vety) - 1]

    f = parse_expr(funct)

    for j in range(len(consts_bashforth[int(ordem) - 1])):
        aux += float(h) * consts_bashforth[int(ordem) - 1][j] * f.subs(t, float(tf) - (float(h) * j)).subs(y, float(vety[len(vety) - j - 1]))
    
    yf = yf + aux

    return yf

def adam_multon(funct, vety, vett, qt_it, h, ordem, complemento = ''):
    print('Adam Multon', complemento)
    print('it = {}'.format(qt_it))
    print('h = {}'.format(h))

    f = parse_expr(funct)
    arrt = vett
    arry = vety

    for i in range(len(arry)):
        print(i, ' ', arry[i])

    for i in range(int(ordem), int(qt_it) + 1, 1):
        aux = 0.0
        pd_term = float(h) * consts_multon[int(ordem) - 1][0] * f.subs(t, float(arrt[i - 1]) + float(h)).subs(y, bashfort_dynamic_p(funct, arry, arrt[len(arrt) - 1], h, ordem))
        aux += pd_term

        for j in range(1, len(consts_multon[int(ordem) - 1]), 1):
            aux += float(h) * consts_multon[int(ordem) - 1][j] * f.subs(t, arrt[i - 1] - (float(h) * (j - 1))).subs(y, arry[len(arry) - (j - 1) - 1])

        arry.append(float(arry[i - 1]) + aux)
        arrt.append(float(arrt[i - 1]) + float(h))


        print(i, ' ', arry[i])

    return arry, arrt

def adam_bashforth(funct, vety, vett, qt_it, h, ordem, complemento = ''):
    print('Adam Bashforth', complemento)
    print('it = {}'.format(qt_it))
    print('h = {}'.format(h))

    f = parse_expr(funct)
    arrt = vett
    arry = vety

    for i in range(len(arry)):
        print(i, ' ', arry[i])

    for i in range(int(ordem), int(qt_it) + 1, 1):
        aux = 0.0

        for j in range(0, len(consts_bashforth[int(ordem) - 1]), 1):
            aux += float(h) * consts_bashforth[int(ordem) - 1][j] * f.subs(t, arrt[i - 1] - (float(h) * j)).subs(y, arry[len(arry) - j - 1])

        arry.append(float(arry[i - 1]) + aux)
        arrt.append(float(arrt[i - 1]) + float(h))


        print(i, ' ', arry[i])

    return arry, arrt

## falta formula inversa!

def main():
    arq = open('input.txt', 'r')

    # euler, euler-inverso, euler-aprimorado, runge-kutta:
    # metodo y0 t0 qt_passos func

    # adam_multon_by_euler, by_euler_inverso, by_runge_kutta, by_euler_aprimorado
    # metodo y0 t0 qt_passos func ordem

    # adam_multon
    # metodo y0 y1 ... yn - 1 t0 qt_passos func ordem

    # adam_bashforth_by_euler, by_euler_inverso, by_runge_kutta, by_euler_aprimorado
    # metodo y0 t0 qt_passos func ordem

    # adam_bashforth
    # metodo y0 y1 ... yn t0 qt_passos func ordem

    # formula inversa
    # ?
     

    for linha in arq:
        param = linha.split()
        if(param[0] == "euler"):
            pts_y, pts_t = euler(param[5], param[1], param[2], param[4], param[3])
            plotGraphic("Euler", pts_t, pts_y, 'red')
        elif(param[0] == "euler_inverso"):
            pts_y, pts_t = reverse_euler(param[5], param[1], param[2], param[4], param[3])
            plotGraphic("Euler Inverso", pts_t, pts_y, 'yellow')
        elif(param[0] == "euler_aprimorado"):
            pts_y, pts_t = aprimorated_euler(param[5], param[1], param[2], param[4], param[3])
            plotGraphic("Euler Aprimorado", pts_t, pts_y, 'blue')
        elif(param[0] == "runge_kutta"):
            pts_y, pts_t = runge_kutta(param[5], param[1], param[2], param[4], param[3])
            plotGraphic("Runge Kutta", pts_t, pts_y, 'black')
        elif("adam_multon" in param[0]):
            if("by_euler_inverso" in param[0]):
                pts_y, pts_t = reverse_euler(param[len(param) - 2], param[1], param[2], int(param[len(param) - 1]) - 2, param[3])
                pts_y, pts_t = adam_multon(param[len(param) - 2], pts_y, pts_t, param[len(param) - 3], param[len(param) - 4], param[len(param) - 1], ' by Euler Inverso')
                plotGraphic("Adam Multon by Euler Inverso", pts_t, pts_y, 'cyan')
            elif("by_euler_aprimorado" in param[0]):
                pts_y, pts_t = aprimorated_euler(param[len(param) - 2], param[1], param[2], int(param[len(param) - 1]) - 2, param[3])
                pts_y, pts_t = adam_multon(param[len(param) - 2], pts_y, pts_t, param[len(param) - 3], param[len(param) - 4], param[len(param) - 1], ' by Euler Aprimorado')
                plotGraphic("Adam Multon by Euler Aprimorado", pts_t, pts_y, 'cyan')
            elif("by_runge_kutta" in param[0]):
                pts_y, pts_t = runge_kutta(param[len(param) - 2], param[1], param[2], int(param[len(param) - 1]) - 2, param[3])
                pts_y, pts_t = adam_multon(param[len(param) - 2], pts_y, pts_t, param[len(param) - 3], param[len(param) - 4], param[len(param) - 1], ' by Runge Kutta')
                plotGraphic("Adam Multon by Runge Kutta", pts_t, pts_y, 'cyan')
            elif("by_euler" in param[0]):
                pts_y, pts_t = euler(param[len(param) - 2], param[1], param[2], int(param[len(param) - 1]) - 2, param[3])
                pts_y, pts_t = adam_multon(param[len(param) - 2], pts_y, pts_t, param[len(param) - 3], param[len(param) - 4], param[len(param) - 1], ' by Euler')
                plotGraphic("Adam Multon by Euler", pts_t, pts_y, 'cyan')
            else:
                pts_y = []
                pts_t = []
                h_aux = 0
                for i in range(int(param[len(param) - 1]) - 1):
                    pts_y.append(float(param[i+1]))
                    pts_t.append(float(param[len(param) - 4]) + h_aux)
                    h_aux = h_aux + float(param[len(param) - 4])
                
                pts_y, pts_t = adam_multon(param[len(param) - 2], pts_y, pts_t, param[len(param) - 3], param[len(param) - 4], param[len(param) - 1])
                plotGraphic("Adam Multon", pts_t, pts_y, 'cyan')
        elif("adam_bashforth" in param[0]):
            if("by_euler_inverso" in param[0]):
                pts_y, pts_t = reverse_euler(param[len(param) - 2], param[1], param[2], int(param[len(param) - 1]) - 2, param[3])
                pts_y, pts_t = adam_bashforth(param[len(param) - 2], pts_y, pts_t, param[len(param) - 3], param[len(param) - 4], param[len(param) - 1], ' by Euler Inverso')
                plotGraphic("Adam Bashforth by Euler Inverso", pts_t, pts_y, 'brown')
            elif("by_euler_aprimorado" in param[0]):
                pts_y, pts_t = aprimorated_euler(param[len(param) - 2], param[1], param[2], int(param[len(param) - 1]) - 2, param[3])
                pts_y, pts_t = adam_bashforth(param[len(param) - 2], pts_y, pts_t, param[len(param) - 3], param[len(param) - 4], param[len(param) - 1], ' by Euler Aprimorado')
                plotGraphic("Adam Bashforth by Euler Aprimorado", pts_t, pts_y, 'brown')
            elif("by_runge_kutta" in param[0]):
                pts_y, pts_t = runge_kutta(param[len(param) - 2], param[1], param[2], int(param[len(param) - 1]) - 2, param[3])
                pts_y, pts_t = adam_bashforth(param[len(param) - 2], pts_y, pts_t, param[len(param) - 3], param[len(param) - 4], param[len(param) - 1], ' by Runge Kutta')
                plotGraphic("Adam Bashforth by Runge Kutta", pts_t, pts_y, 'brown')
            elif("by_euler" in param[0]):
                pts_y, pts_t = euler(param[len(param) - 2], param[1], param[2], int(param[len(param) - 1]) - 2, param[3])
                pts_y, pts_t = adam_bashforth(param[len(param) - 2], pts_y, pts_t, param[len(param) - 3], param[len(param) - 4], param[len(param) - 1], ' by Euler')
                plotGraphic("Adam Bashforth by Euler", pts_t, pts_y, 'brown')
            else:
                pts_y = []
                pts_t = []
                h_aux = 0
                for i in range(int(param[len(param) - 1])):
                    pts_y.append(float(param[i+1]))
                    pts_t.append(float(param[len(param) - 4]) + h_aux)
                    h_aux = h_aux + float(param[len(param) - 4])
                
                pts_y, pts_t = adam_bashforth(param[len(param) - 2], pts_y, pts_t, param[len(param) - 3], param[len(param) - 4], param[len(param) - 1])
                plotGraphic("Adam Bashforth", pts_t, pts_y, 'brown')

main()