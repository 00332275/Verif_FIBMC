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

with open('arezoumandi-rca-4.xlsx', 'rb') as target:
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
Nad = np.single(data[0,13])       # coluna N: Nad
Maxd = np.single(data[0,14])      # coluna O: Maxd
#
# Final da entrada de dados
#
################################################
# cálculo das const. para diag. de tensão de compressão
#
cc   = 1e5/(10**np.floor(np.log10(Es)))
fck  = fck*cc
#fcm = (fck/1.4*0.85)/cc
fcm  = (fck + 8)/cc

    # determinação do k
k_i = pd.read_csv('kvalues.csv')
xk  = k_i['X'].values
yk  = k_i['Y'].values
k4  = np.poly1d(np.polyfit(xk, yk, 4))
k   = np.round(k4(fck),2)

    # determinação de epsc limite
if fck<50:
    epsc_lim = -3.5/1000
else:
    epsc_lim = np.round(0.01*fck - 3.9,1)/1000

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
        epsc[i+1] = (epsc[i] + epsc_lim/1000)
eta = epsc/epsc1
sigmac_i = -(k*eta - eta**2)/(1+(k-2)*eta)   #sigmac_i/fcm
ajust4 = np.poly1d(np.polyfit(epsc, sigmac_i, 4))
a0 = ajust4[0]
a1 = ajust4[1]
a2 = ajust4[2]
a3 = ajust4[3]
a4 = ajust4[4]

"""
# plot gráficos
sigmac4 = ajust4(epsc)
fig, ax = plt.subplots()
line_down, = ax.plot(-epsc, -sigmac_i, 'b', label='Eq. FIB model code')
line_up, = ax.plot(-epsc, -sigmac4, 'r--', label='Polinômio 4a ordem')
ax.legend(handles=[line_up, line_down])
ax.set_xlabel('concrete strain εc<0')
ax.set_ylabel('concrete stress σc<0[MPa]');
"""

################################################
# cálculo das const. para diag. de tensão de tração
if fck>50:
    fctm = (2.12*np.log(1 + 0.1*(fck + 8)))/cc
else:
    fctm = (0.3*(fck)**(2/3))/cc

alfE = 1.0 # granito e gnaisse (1.2=basalto, 0.9=calcário, 0.7=arenito)
Eci  = (21.5 * 10**3 * alfE * ((fck + 8)/10)**(1/3))/cc
b1   = Eci

epsct0 = 0.9*fctm/Eci
c1   = 0.1*fctm / (0.15/1000 - epsct0)
c0   = fctm * (0.9 - (0.1*epsct0)/(0.15/1000 - epsct0))


################################################
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
    rj() - porcentagem de armaduras das barras
    Es - módulo de def. longitudinal do aço
"""
def esfor(As, a1, a2, a3, a4, b1, c0, c1, Es, epscu, epsc2, epsct0, fcm, fyk,
          Ly, rj, X, xc, xg, yc, yg, ys):
    YS = np.max(yc)
    rjj= np.copy(rj)
    for i in range(np.size(rjj)):
        if ys[i] > 0:
            rjj[i] = 0
        else:
            continue
    rjj2 = rjj/(np.sum(rjj))
    YI = np.sum(rjj2*ys)
    d  = YS-YI
    # calculo de epsilon S e epsilon I
    # domínios 1 e 2
    if X>(-1E50) and X<=((epscu*d)/(50/(1000)+epscu)):
        epsS = -50/(1000)*X/(d-X)
        epsI = 50/(1000)
    # domínios 3 e 4
    elif X>=((epscu*d)/(50/(1000)+epscu)) and X<=(d):
        epsS = -epscu
        epsI = epscu*(d-X)/X
    # domínio 4a
    elif X>=(d) and X<=(Ly):
        epsS = -epscu
        epsI = 0
    # domínio 5
    else:
        epsS = -epsc2*(X/(X-Ly*((epscu-epsc2)/epscu)))
        epsI = -epsc2*((X-Ly)/(X-Ly*((epscu-epsc2)/epscu)))

    b = (epsS-epsI)/d
    c = epsS - b*YS

    xx1, yy1 = poli(xc, yc, c, b)
    Mrx1, Nr1 = reg1(a1, a2, a3, a4, c, b, fcm, xx1, yy1)
    
    # Param. abertura de fissuras
    kk   = np.argmin(ys)
    As1  = rj[kk]*As
    phi1 = np.sqrt(As1*4/np.pi)
    cob  = np.abs(np.min(yc) - np.min(ys)) - phi1/2
    Ast  = As*(np.sum(rjj))
    if 2.5*(Ly-d)<(Ly-X)/3:
        Ac_ef = 2.5*(Ly-d)
    else:
        Ac_ef = (Ly-X)/3
    rhos_ef = Ast/Ac_ef
    lsmax = 1.0*cob + 0.25*phi1/(1.8*rhos_ef)
    Gf  = (73 * (fcm*200000/Es)**0.18)/1e5 #kN/cm
    w1  = Gf/(fctm)
    k1  = 0.8*fctm/(0.15/1000 - w1/lsmax)
    k0  = fctm - 0.15/1000*k1
    kk1 = - 0.2*fctm/(4*w1/lsmax)
    kk0 = - 5*kk1*w1/lsmax
    #
    if epsI>=0:
        y01 = -c/b
        y12 = (epsct0 - c)/b
        y23 = ((0.15/1000) - c)/b
        y34 = ((w1/lsmax) - c)/b
        y44 = ((5*w1/lsmax) - c)/b
        xxt1, yyt1 = polit(xc, yc, y01, y12)
        xxt2, yyt2 = polit(xc, yc, y12, y23)
        xxt3, yyt3 = polit(xc, yc, y23, y34)
        xxt4, yyt4 = polit(xc, yc, y34, y44)
    else:
        xxt1 = np.zeros(np.size(xc))
        yyt1 = xxt2 = yyt2 = xxt3 = yyt3 = xxt4 = yyt4 = xxt1
        
    Mrxt1, Nrt1 = regt(b, c, b1, 0, xxt1, yyt1)
    Mrxt2, Nrt2 = regt(b, c, c1, c0, xxt2, yyt2)
    Mrxt3, Nrt3 = regt(b, c, k1, k0, xxt3, yyt3)
    Mrxt4, Nrt4 = regt(b, c, kk1, kk0, xxt4, yyt4)
    
    Mrxas, Nras = aco(As, b, c, Es, fyk, rj, ys)
    
    Mrxd = Mrx1 + (Mrxt1 + Mrxt2 + Mrxt3 + Mrxt4) + Mrxas
    Nrd = Nr1 + (Nrt1 + Nrt2 + Nrt3 + Nrt4) + Nras
        
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
def poli(xc, yc, c, b):
    y01 = (-c/b)
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
    
    yy1=yy1[np.logical_not(np.isnan(yy1))]
    xx1=xx1[np.logical_not(np.isnan(xx1))]
    return (xx1, yy1)
#############################

#############################
"""
Função polit() retorna as variáveis xxt1,yyt1 que determinam
as poligonais da região correspondente à área tracionada da seção.
Variáveis:
    xc, yc - vértices da seção em relação ao cg.
    xxt1, yyt1 - vértices das poligonais das região
    c, b - variáveis de ajuste usadas em D0, D1 e D2
    epsct0 - def. específica limite do 1o trecho
    y01 - ordenada limite entre as regiões 0 e 1
    y12 - ordenada limite entre as regiões 1 e 2
    y22 - ordenada limite da região 2 (epsctlim = 0.15/1000)

Last-Modified: 28/08/2022
Status: 
"""
def polit(xc, yc, y01, y12):
    xxt1 = np.copy(xc)
    yyt1 = np.copy(yc)

    if y01<np.min(yc):
        yyt1 = xxt1 = np.array([])
    else:
        for i in range(np.size(yyt1)):
            if yyt1[i]>y01:
                yyt1[i] = np.nan
                xxt1[i] = np.nan
            elif np.isnan(yyt1[i-1])==True:
                yyt1[i-1] = y01
                xxt1[i-1] = (y01-yc[i-1])/(yc[i]-yc[i-1])*(xc[i]-xc[i-1])+xc[i-1]
            elif i<(np.size(yyt1)-1):
                if (yyt1[i+1]>y01) and (yyt1[i]<y01):
                    yyt1[i+1] = y01
                    xxt1[i+1] = (y01-yc[i])/(yc[i+1]-yc[i])*(xc[i+1]-xc[i])+xc[i]
            else:
                continue
        for i in range(np.size(yyt1)):
            if yyt1[i]<y12:
                yyt1[i] = np.nan
                xxt1[i] = np.nan
            elif np.isnan(yyt1[i-1])==True:
                yyt1[i-1] = y12
                xxt1[i-1] = (y12-yc[i-1])/(yc[i]-yc[i-1])*(xc[i]-xc[i-1])+xc[i-1]
            elif i<(np.size(yyt1)-1):
                if (yyt1[i+1]<y12) and (yyt1[i]>y12):
                    yyt1[i+1] = y12
                    xxt1[i+1] = (y12-yc[i])/(yc[i+1]-yc[i])*(xc[i+1]-xc[i])+xc[i]
            else:
                continue
    yyt1=yyt1[np.logical_not(np.isnan(yyt1))]
    xxt1=xxt1[np.logical_not(np.isnan(xxt1))]
    
    return (xxt1, yyt1)
#############################

#############################
"""
Função reg1() integra a região de compressão do concreto, retornando momento
e esforço normal resistente.
Variáveis:
    D0, D1,..., D4 - coef. geométricos da seção p/ calc. das tensões no concreto
    scd = sigma cd, tensão de cálculo do concreto
    Mrx1, Nr1 - esforços resistentes da reg. comprimida
    G00, G01,..., G05 - polinômios de integração
"""
def reg1(a1, a2, a3, a4, c, b, fcm, xx1, yy1):
    D0   = a1*c + a2*c**2 + a3*c**3 + a4*c**4
    D1   = a1*b + 2*a2*b*c + 3*a3*b*c**2 + 4*a4*b*c**3
    D2   = a2*b**2 + 3*a3*b**2*c + 6*a4*b**2*c**2
    D3   = a3*b**3 + 4*a4*b**3*c
    D4   = a4*b**4 
    scd  = fcm
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
            dy4  = dy3*dy
            dy5  = dy4*dy
            dy6  = dy5*dy
            G00  = (x1 + dx/2)*dy
            G01  = (x1*(y1+dy1)+dx*(y1/2+dy/3))*dy
            G02  = (x1*(y1*(dy+y1)+dy2/3) + dx*(y1*(y1/2+dy/1.5)+dy2/4))*dy
            G03  = (x1*(y1*(dy2+y1*(1.5*dy+y1))+dy3/4) + 
                    dx*(y1*(0.75*dy2+y1*(dy+y1/2))+dy3/5))*dy
            G04  = (x1*y1**4*dy + dx*dy5/6 + dy4*(x1*dy+4*y1*dx)/5 +
                    dy3*(2*x1*y1*dy + 3*y1**2*dx)/2 + dy2*(6*x1*y1**2*dy + 4*y1**3*dx)/3 +
                    dy*(4*x1*y1**3*dy + y1**4*dx)/2)
            G05  = (x1*y1**5*dy + dx*dy6/7 + dy5*(x1*dy + 5*y1*dx)/6 +
                    dy4*(x1*y1*dy + 2*y1**2*dx) + dy3*(5*x1*y1**2*dy + 5*y1**3*dx)/2 +
                    dy2*(10*x1*y1**3*dy + 5*y1**4*dx)/3 + dy*(5*x1*y1**4*dy + y1**5*dx)/2)
            Mrx1 = Mrx1 + scd*(D0*G01 + D1*G02 + D2*G03 + D3*G04 + D4*G05)
            Nr1  = Nr1 + scd*(D0*G00 + D1*G01 + D2*G02 + D3*G03 + D4*G04)
    return (Mrx1, Nr1)
#############################
#############################
"""
Função regt() integra as regiões de tração/fissuração do concreto, definidas
por funções lineares, retornando momento e esforço normal resistente.
Variáveis:
    Mrxt, Nrt - esforços resistentes da regiao
    G00, G01, G02 - polinômios de integração
"""
def regt(b, c, k1, k0, xt, yt):
    D0    = k1*c + k0
    D1    = k1*b
    Mrxt = 0
    Nrt  = 0
    for i in range(np.size(xt)-1):
        x1 = xt[i]
        y1 = yt[i]
        x2 = xt[i+1]
        y2 = yt[i+1]
        dx = x2 - x1
        dy = y2 - y1
        dy1  = dy/2 
        if dy==0:
            continue
        else:
            dy2  = dy*dy
            G00  = (x1 + dx/2)*dy
            G01  = (x1*(y1+dy1)+dx*(y1/2+dy/3))*dy
            G02  = (x1*(y1*(dy+y1)+dy2/3) + dx*(y1*(y1/2+dy/1.5)+dy2/4))*dy
            Mrxt = Mrxt + D0*G01 + D1*G02
            Nrt  = Nrt + D0*G00 + D1*G01
    return (Mrxt, Nrt)
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
def aco(As, b, c, Es, fyk, rj, ys):
    Nras  = 0
    Mrxas = 0
    epsuk = 5/100
    epsyk = fyk/Es
    a1 = (0.08*fyk)/(epsuk - epsyk)
    a0 = fyk - a1*epsyk
    for i in range(np.size(ys)):
        epsb = b*ys[i] + c
        if (np.abs(epsb)<=epsyk):
            sig = Es*epsb
        else:
            sig = a1*epsb + a0
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
    a3    = const[3]
    a4    = const[4]
    b1    = const[5]
    c0    = const[6]
    c1    = const[7]
    Es    = const[8]
    fcm   = const[9]
    fyk   = const[10]
    Ly    = const[11]
    Maxd  = const[12]
    Nad   = const[13]
    epsct0= const[14]
    #print(X,lamb)
    Mrxd, Nrd, epsS, epsI = esfor(As, a1, a2, a3, a4, b1, c0, c1, Es, epscu, epsc2, epsct0, fcm, fyk,
                                  Ly, rj, X, xc, xg, yc, yg, ys)
    f = lamb*(Mrxd) - Maxd
    g = lamb*(Nrd) - Nad
    return[f,g]

const = np.array([As, a1, a2, a3, a4, b1, c0, c1, Es, fcm, fyk, Ly, Maxd, Nad, epsct0])
lamb_i = 1 #lambda inicial
X_i = 0.5*Ly #altura da LN inicial
s0 = np.array([X_i, lamb_i])
var_aux = (const, -epsc_lim, -epsc1, rj, xc, xg, yc, yg, ys)
X, lamb  = fsolve(nlsistema, s0, var_aux)
Mrxd, Nrd, epsS, epsI = esfor(As, a1, a2, a3, a4, b1, c0, c1, Es, -epsc_lim, -epsc1, epsct0, fcm, fyk,
                              Ly, rj, X, xc, xg, yc, yg, ys)
FS = 1/lamb


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
