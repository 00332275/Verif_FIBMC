# -*- coding: utf-8 -*-
"""
Programa para verificação de seções poligonais de concreto armado submetidas 
à flexão composta normal segundo o FIB model code (2010).

@author: GabrielMachado
Last-Modified: 02/08/2022
"""

import numpy as np
import pandas as pd
# Para resolver o sistema de eq. não lineares
from scipy.optimize import fsolve
import xlsxwriter
from datetime import datetime
import matplotlib.pyplot as plt


# Entrada de dados
#

with open('dados.xlsx', 'rb') as target:
    sheet =  pd.read_excel(target, sheet_name='Planilha1')
    data  =  sheet.values
#
nc = np.int_( data[0,0] )         # coluna A: nc
nc1 = nc+1
xc = np.zeros(nc1)
yc = np.zeros(nc1)
xc[0:nc] = data[0:nc,1]           # coluna B: xc
xc[-1] = data[0,1]
yc[0:nc] = data[0:nc,2]           # coluna C: yc
yc[-1] = data[0,2]
As = np.single(data[0,3])         # coluna D: As
ns = np.int_(data[0,4])           # coluna E: ns
xs = np.zeros(ns)
ys = np.zeros(ns)
xs = data[0:ns,5]                 # coluna F: xs
ys = data[0:ns,6]                 # coluna G: ys
rj = data[0:ns,7]                 # coluna H: rj
fck = np.single(data[0,8])        # coluna I: fck
fyk = np.single(data[0,9])        # coluna J: fyk
Es = np.single(data[0,10])        # coluna K: Es
gamac = np.single(data[0,11])     # coluna L: gamac
gamas = np.single(data[0,12])     # coluna M: gamas
Nad = np.single(data[0,13])       # coluna N: Nad
Maxd = np.single(data[0,14])      # coluna O: Maxd
#
# Final da entrada de dados
#
################################################
# cálculo das const. para diag. stress-strain
#
cc  = 200000/Es
fck = fck*cc
fcm = fck + 8

    # determinação do k
k_i = pd.read_csv('kvalues.csv')
xk  = k_i['X'].values
yk  = k_i['Y'].values
k4  = np.poly1d(np.polyfit(xk, yk, 4))
k   = np.round(k4(fck),2)

    # determinação de epsc limite
if fck<50:
    epsc_lim = 3.5/1000
else:
    epsc_lim = np.round(-0.01*fck + 3.9,1)/1000

    # determinação de epsc1
epsc1_i = pd.read_csv('epsc1values.csv')
xeps  = epsc1_i['X'].values
yeps  = epsc1_i['Y'].values
epsc1_6  = np.poly1d(np.polyfit(xeps, yeps, 6))
epsc1   = np.round(epsc1_6(fck),1)/1000

    # ajuste de sigma_c
epsc = np.zeros(1001)
for i in range(1001):
    if i==1000:
        continue
    else:
        epsc[i+1] = epsc[i] + epsc_lim/1000
eta = epsc/epsc1
sigmac_i = fcm*(k*eta - eta**2)/(1+(k-2)*eta)
ajust4 = np.poly1d(np.polyfit(epsc, sigmac_i, 4))
a0 = ajust4[0]
a1 = ajust4[1]
a2 = ajust4[2]
a3 = ajust4[3]
a4 = ajust4[4]

# plot gráficos
sigmac4 = ajust4(epsc)
fig, ax = plt.subplots()
line_down, = ax.plot(epsc, sigmac_i, 'b', label='Eq. FIB model code')
line_up, = ax.plot(epsc, sigmac4, 'r--', label='Polinômio 4a ordem')
ax.legend(handles=[line_up, line_down])
ax.set_xlabel('concrete strain εc<0')
ax.set_ylabel('concrete stress σc<0[MPa]');

#
xmax = np.amax(xc, axis=0)
xmin = np.amin(xc, axis=0)
ymax = np.amax(yc, axis=0)
ymin = np.amin(yc, axis=0)

# cálculo das propriedades geométricas da seção
Lx = xmax-xmin   #b
Ly = ymax-ymin   #h
Ac = 0
Sx = 0
Jx = 0
Sy = 0
Jy = 0
Jxy = 0
for i in range(nc):
    dx = xc[i+1]-xc[i]
    dy = yc[i+1]-yc[i]
    Ac = Ac + (xc[i]+dx/2)*dy
    Sx = Sx + (xc[i]*(yc[i]+dy/2) + dx*(yc[i]/2+dy/3))*dy
    Jx = Jx + (xc[i]*(yc[i]*(dy+yc[i])+dy*dy/3) + \
               dx*(yc[i]*(yc[i]/2+dy/1.5) + dy*dy/4))*dy
    Sy = Sy + (xc[i]*(xc[i]+dx)+dx*dx/3)*dy/2
    Jy = Jy + (dx**3/4+xc[i]*(dx*dx+xc[i]*(1.5*dx+xc[i])))*dy/3
    Jxy = Jxy + (xc[i]*(xc[i]*(yc[i]+dy/2) + dx*(yc[i]+dy/1.5)) + \
        dx*dx*(yc[i]/3+dy/4))*dy/2
xg   = Sy/Ac
yg   = Sx/Ac
Sxg  = Sx - yg*Ac
Syg  = Sy - xg*Ac
Jxg  = Jx - Ac*yg**2
Jyg  = Jy - Ac*xg**2
Jxyg = Jxy - Ac*xg*yg

# coordenadas no sistema local xg,yg
xc = xc - xg
yc = yc - yg
xs = xs - xg
ys = ys - yg

###################################################

#############################
"""
Função esfor() retorna os esforços resistentes MRxd e NRd, calculados
a partir das funções poli(), reg1(), reg2() - integração da área de 
concreto comprimida -  e aco() - diagrama tensão-deformação das armaduras.
Variáveis:
    X - altura da linha neutra
    Ly - altura da seção
    xc, yc - vértices da seção em relação ao cg.
    ys - coordenadas y das armaduras em relação ao cg.
    epscu, epsc2 - def. específicas do diagrama parábola-retângulo 
            (NBR 6118:2014)
    xg, yg - coord. do centroide
    YS, YI - alturas acima e abaixo da LN
    d - altura útil da seção
    epsS, epsI - def. específicas nos pontos extremos da seção
    b, c - variáveis de ajuste usadas em D0, D1 e D2
    a1, a2 - coeficientes do ajuste por regressão do diagr. parábola/ret
    fcd - resistência à compressão do concreto em projeto
    fyk - resistência à tração característica do aço
    gamas - coeficiente de segurança aço
    rj() - porcentagem de armaduras das barras
    Es - módulo de def. longitudinal do aço
"""
def esfor(As, a1, a2, Es, epscu, epsc2, fcd, fyk, gamas, Ly, rj, X, xc, xg, yc, yg, ys):
    YS = np.max(yc)
    YI = np.min(ys)
    d  = YS-YI
    # calculo de epsilon S e epsilon I
    if X>(-1E50) and X<=((epscu*d)/(10/1000+epscu)):
        epsS = -10/1000*X/(d-X)
        epsI = 10/1000
    elif X>=((epscu*d)/(10/1000+epscu)) and X<=(d):
        epsS = -epscu
        epsI = epscu*(d-X)/X
    elif X>=(d) and X<=(Ly):
        epsS = -epscu
        epsI = 0
    else:
        epsS = -epsc2*(X/(X-Ly*((epscu-epsc2)/epscu)))
        epsI = -epsc2*((X-Ly)/(X-Ly*((epscu-epsc2)/epscu)))

    b = (epsS-epsI)/(YS-YI)
    c = epsS - b*YS

    xx1, yy1, xx2, yy2 = poli(xc, yc, c, b, epsc2)
    Mrx1, Nr1 = reg1(a1, a2, c, b, fcd, xx1, yy1)
    Mrx2, Nr2 = reg2(fcd, xx2, yy2)
    
    fyd = fyk/gamas
    Mrxas, Nras = aco(As, b, c, Es, fyd, rj, ys)
    
    Mrxd = Mrx1 + Mrx2 + Mrxas
    Nrd = Nr1 + Nr2 + Nras
    
    return (Mrxd, Nrd, epsS, epsI)
#############################

#############################
"""
Função poli() retorna as variáveis xx1,yy1,xx2,yy2 que determinam
as poligonais das regiões 1 e 2 da área comprimida da seção.
Variáveis:
    xc, yc - vértices da seção em relação ao cg.
    xx1, yy1, xx2, yy2 - vértices das poligonais das reg. 1 e 2
    c, b - variáveis de ajuste usadas em D0, D1 e D2
    epsc2 - def. específica c2 do diagrama parábola-retângulo 
            (NBR 6118:2014)
    y01 - ordenada limite entre as regiões 0 e 1
    y12 - ordenada limite entre as regiões 1 e 2

Last-Modified: 17/07/2022
Status: Ok - Testado para seções retangulares e T.
"""
def poli(xc, yc, c, b, epsc2):
    y01 = (-c/b)
    y12 = (-epsc2-c)/b

    xx1 = np.copy(xc)
    yy1 = np.copy(yc)

    for i in range(np.size(yy1)):
        if yy1[i]<y01:
            yy1[i] = np.nan
            xx1[i] = np.nan
        elif np.isnan(yy1[i-1])==True:
            yy1[i-1] = y01
            xx1[i-1] = (y01-yc[i-1])/(yc[i]-yc[i-1])*(xc[i]-xc[i-1])+xc[i-1]
        elif i<(np.size(yy1)-1):
            if (yy1[i+1]<y01) and (yy1[i]>y01):
                yy1[i+1] = y01
                xx1[i+1] = (y01-yc[i])/(yc[i+1]-yc[i])*(xc[i+1]-xc[i])+xc[i]
        else:
            continue
    
    for i in range(np.size(yy1)):
        if i<(np.size(yy1)-1):
            if yy1[i]<y12 and yy1[i+1]>y12:
                yy1[i+1] = y12
                xx1[i+1] = (y12-yc[i])/(yc[i+1]-yc[i])*(xc[i+1]-xc[i])+xc[i]
            elif yy1[i]>y12 and yy1[i+1]<y12:
                yy1[i] = y12
                xx1[i] = (y12-yc[i+1])/(yc[i]-yc[i+1])*(xc[i]-xc[i+1])+xc[i+1]
            elif yy1[i]>y12:
                yy1[i] = np.nan
                xx1[i] = np.nan
            else:
                continue
    yy1=yy1[np.logical_not(np.isnan(yy1))]
    xx1=xx1[np.logical_not(np.isnan(xx1))]

    xx2 = np.copy(xc)
    yy2 = np.copy(yc)

    for i in range(np.size(yy2)):
        if yy2[i]<y12:
            yy2[i] = np.nan
            xx2[i] = np.nan
        elif np.isnan(yy2[i-1])==True:
            yy2[i-1] = y12
            xx2[i-1] = (y12-yc[i-1])/(yc[i]-yc[i-1])*(xc[i]-xc[i-1])+xc[i-1]
        elif i<(np.size(yy2)-1):
            if (yy2[i+1]<y12) and (yy2[i]>y12):
                yy2[i+1] = y12
                xx2[i+1] = (y12-yc[i])/(yc[i+1]-yc[i])*(xc[i+1]-xc[i])+xc[i]
        else:
            continue
    yy2=yy2[np.logical_not(np.isnan(yy2))]
    xx2=xx2[np.logical_not(np.isnan(xx2))]
    
    return (xx1, yy1, xx2, yy2)
#############################

#############################
"""
Função reg1() integra a região 1 de compressão do concreto, retornando momento
e esforço normal resistente.
Variáveis:
    D0, D1, D2 - coef. geométricos da seção p/ calc. das tensões no concreto
    scd = sigma cd, tensão de cálculo do concreto
    Mrx1, Nr1 - esforços resistentes da reg. 1
    G00, G01, G02, G03 - polinômios de integração
"""
def reg1(a1, a2, c, b, fcd, xx1, yy1):
    D0   = a1*c + a2*c**2
    D1   = a1*b + 2*a2*b*c
    D2   = a2*b**2
    scd  = 0.85*fcd
    Mrx1 = 0
    Nr1  = 0
    for i in range(np.size(xx1)-1):
        x1 = xx1[i]
        y1 = yy1[i]
        x2 = xx1[i+1]
        y2 = yy1[i+1]
        dx = x2 - x1
        dy = y2 - y1
        dy1  = dy/2
        if dy==0:
            continue
        else:
            dy2  = dy*dy
            dy3  = dy2*dy
            G00  = (x1 + dx/2)*dy
            G01  = (x1*(y1+dy1)+dx*(y1/2+dy/3))*dy
            G02  = (x1*(y1*(dy+y1)+dy2/3) + dx*(y1*(y1/2+dy/1.5)+dy2/4))*dy
            G03  = (x1*(y1*(dy2+y1*(1.5*dy+y1))+dy3/4) + 
                    dx*(y1*(0.75*dy2+y1*(dy+y1/2))+dy3/5))*dy
            Mrx1 = Mrx1 + scd*(D0*G01 + D1*G02 + D2*G03)
            Nr1  = Nr1 + scd*(D0*G00 + D1*G01 + D2*G02)
    return (Mrx1, Nr1)
#############################
#############################
"""
Função reg2() integra a região 2 de compressão do concreto, retornando momento
e esforço normal resistente.
Variáveis:
    scd = sigma cd, tensão de cálculo do concreto
    Mrx2, Nr2 - esforços resistentes da reg. 2
    G00, G01, G02, G03 - polinômios de integração
"""
def reg2(fcd, xx2, yy2):
    scd  = 0.85*fcd
    Mrx2 = 0
    Nr2  = 0
    for i in range(np.size(xx2)-1):
        x1 = xx2[i]
        y1 = yy2[i]
        x2 = xx2[i+1]
        y2 = yy2[i+1]
        dx = x2 - x1
        dy = y2 - y1
        dy1  = dy/2 
        if dy==0:
            continue
        else:
            G00  = (x1 + dx/2)*dy
            G01  = (x1*(y1+dy1)+dx*(y1/2+dy/3))*dy
            Mrx2 = Mrx2 - scd*G01
            Nr2  = Nr2 - scd*G00
    return (Mrx2, Nr2)
#############################
#############################
"""
Função aco() retorna a parcela resistente das armaduras da seção.
Variáveis:
    As - área total de aço
    b,c - variáveis de ajuste
    Nrasi - força normal na barra i
    Nras, Mrxas - esforços resistentes do aço
    sig - tensão na armadura
    epsb - deformação específica da barra
    fyd - resistência de projeto do aço
    Es - módulo elástico do aço
    rj - porcentagem de área de aço das armaduras
    ys - ordenadas das armaduras em relação ao cg
"""
def aco(As, b, c, Es, fyd, rj, ys):
    Nras  = 0
    Mrxas = 0
    for i in range(np.size(ys)):
        epsb = b*ys[i] + c
        epsyd = fyd/Es
        if (np.abs(epsb)<=epsyd):
            sig = Es*epsb
        else:
            sig = np.sign(epsb)*fyd
        Nrasi = rj[i]*As*sig
        Nras  = Nras + Nrasi
        Mrxas = Mrxas + Nrasi*ys[i]
    return (Mrxas, Nras)
#############################
"""
Função nlsistema() define o sistema de equações não lineares a ser solucionado, 
a partir da chamada das demais funções
"""

def nlsistema(var, *var_aux):
    (X, lamb) = var
    const, epscu, epsc2, rj, xc, xg, yc, yg, ys = var_aux
    As    = const[0]
    a1    = const[1]
    a2    = const[2]
    Es    = const[3]
    fcd   = const[4]
    fyk   = const[5]
    gamas = const[6]
    Ly    = const[7]
    Maxd  = const[8]
    Nad   = const[9]
    Mrxd, Nrd, epsS, epsI = esfor(As, a1, a2, Es, epscu, epsc2, fcd, fyk,
                                  gamas, Ly, rj, X, xc, xg, yc, yg, ys)
    f = lamb*(Mrxd) - Maxd
    g = lamb*(Nrd) - Nad
    return[f,g]
"""
const = np.array([As, a1, a2, Es, fcd, fyk, gamas, Ly, Maxd, Nad])
lamb_i = 1 #lambda inicial
X_i = 0.5*Ly #altura da LN inicial
s0 = np.array([X_i, lamb_i])
var_aux = (const, epscu, epsc2, rj, xc, xg, yc, yg, ys)
#X, lamb  = fsolve(nlsistema, s0, var_aux)
#Mrxd, Nrd, epsS, epsI = esfor(As, a1, a2, Es, epscu, epsc2, fcd, fyk,
#                              gamas, Ly, rj, X, xc, xg, yc, yg, ys)
#FS = 1/lamb


# Exportar resultados para Excel
out = xlsxwriter.Workbook('saída.xlsx')
worksheet = out.add_worksheet()
lista = (
    ['ESFORÇOS SOLICITANTES:',''],
    ['NAd:', Nad],
    ['Maxd:', Maxd],
    ['',''],
    ['ESFORÇOS RESISTENTES:', ''],
    ['NRd:', Nrd],
    ['MRxd:', Mrxd],
    ['',''],
    ['RESULTADOS:', ''],
    ['Altura da LN:', X],
    ['Fator de segurança:', FS],
    ['Deformação na fibra superior: ', epsS],
    ['Deformação na fibra inferior: ', epsI],
)
worksheet.set_column('B:B', 50)
cell_format1 = out.add_format({'bold': True, 'font_size': 18, 'top': True, 'bottom': True})
cell_format2 = out.add_format({'bold': True, 'italic': True ,'font_size': 12, 'top': True, 'bottom': True})
cell_format3 = out.add_format({'font_size': 12, 'top': True, 'bottom': True})
worksheet.write(1, 1, 'RELATÓRIO DE VERIFICAÇÃO DA SEÇÃO:', cell_format1)
worksheet.write(1, 2, '', cell_format1)
row=3
col=1

for i, j in (lista):
    worksheet.write(row, col, i, cell_format2)
    worksheet.write(row, col+1, j, cell_format3)
    row += 1

now = datetime.now()
dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
worksheet.write(row+1, col, dt_string, cell_format3)
worksheet.write(row+1, col+1, '', cell_format3)
out.close()
#

print('X=', X, '\nFS=', FS, '\nMrxd=', Mrxd, '\nNrd=', Nrd)
"""