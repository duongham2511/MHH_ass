# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.

import math

vent_lee = 0
vent_wind = 0

CO2_out = 412.89
cap_CO2_Air = 3.8
cap_CO2_Top = 0.4

P_Blow = 0
n_HeatCO2=0.057
U_Blow=0

A_Flr = 1.4 * 10**4
U_Pad = 0
phi_Pad = 0

O_ExtCO2=7.2*10**4
U_ExtCO2=0

t_air = 20
t_top = 20
t_out = 18
p_air = 1
p_top = 1
g = 9.81
k_th_scr = 0.05*10**(-3)
u_th_scr = 0

eta_Roof_Thr = 0.9
eta_Side_Thr = eta_Roof_Thr

eta_Roof = 1
eta_Side = 0

zeta_InsScr = 1
v_wind = 2.4
c_leakage = 10**-4

C_d = 0.75
C_w = 0.09
h_sideroof = 0
u_roof = (vent_lee + vent_wind) / 2
u_side = 0
a_roof = 0.3 * A_Flr
a_side = 0

U_VentForced = 0
phi_VentForced = 0

M_CH20 = 30*10**(-3)
CBuf = 0
CMaxBuf = 20*10**3

theta=0.7
alpha=0.385
PAR_Can = 100
LAI = 2
J_MAX_Leaf = 210
E_j = 37 * 10**3
T_CanK = 293.2
T_K = 298.15
R = 8.314
S = 710
H = 22 * 10**4

def dx(CO2_Air,CO2_Top):
    dot_CO2_Air = (MC_BlowAir() + MC_ExtAir() + MC_PadAir(CO2_out,CO2_Air) - MC_AirTop(CO2_Air,CO2_Top) - MC_AirOut(CO2_Air,CO2_out) - MC_AirCan(CO2_Air))/cap_CO2_Air
    dot_CO2_Top = (MC_AirTop(CO2_Air,CO2_Top) - MC_TopOut(CO2_out,CO2_Top)) / cap_CO2_Top
    return dot_CO2_Air,dot_CO2_Top

#equation 3
def MC_BlowAir():
    return (n_HeatCO2*U_Blow*P_Blow)/A_Flr

#equation 4
def MC_ExtAir():
    return (O_ExtCO2*U_ExtCO2)/A_Flr

#equation 5
def MC_PadAir(CO2_out,CO2_Air):
    return U_Pad *phi_Pad/A_Flr*(CO2_out-CO2_Air)

#equation 15
def MC_TopOut(CO2_out,CO2_top):
    return fVentRoof()*(CO2_top-CO2_out)

#equation 6
def MC_AirTop(co2_air,co2_top):
    return F_Th_Scr()*(co2_air - co2_top)

#equation 9
def MC_AirOut(co2_air, co2_out):
    return (f_VentSide() + f_VentForced()) * (co2_air - co2_out)

def MC_AirCan(CO2_Air):
    return M_CH20 * hCBuf() * (P(CO2_Air) - R_func(CO2_Air))

#equation 7
def F_Th_Scr():
    return u_th_scr * k_th_scr * math.pow(math.fabs(t_air - t_top),2/3) + (1-u_th_scr) * math.sqrt(g*(1-u_th_scr)*math.fabs(p_air-p_top)/(2 * (p_air+p_top)/2))

#equation 10
def f_vent_roof_side(u_roof,u_side,a_roof,a_side):
    t_mean_air = (t_air + t_out)/2
    a = u_roof**2 * u_side**2 * a_roof**2 * a_side**2 / (u_roof**2 * a_roof**2 + u_side**2 * a_side**2) if u_roof**2 * a_roof**2 + u_side**2 * a_side**2 !=0 else 0
    return C_d/A_Flr * math.sqrt(2 * g * h_sideroof * a*(t_air - t_out)/t_mean_air + ((u_roof*a_roof + u_side*a_side)/2)**2 * C_w * v_wind**2)

#equation 11
def eta_InsScr():
    return zeta_InsScr * (2 - zeta_InsScr)

#equation 12
def f_leakage():
    if (v_wind < 0.25):
        return 0.25 * c_leakage
    else:
        return v_wind * c_leakage

#equation 14
def f_VentForced():
    return (eta_InsScr() * U_VentForced * phi_VentForced) / A_Flr


# cong thuc 17
def fVentRoof2():
    t_mean_air = (t_air + t_out)/2
    fVentRoof2 = (C_d*u_roof*a_roof)/(2*A_Flr)*(g*h_sideroof*(t_air-t_out)/(2*t_mean_air)+C_w*pow(v_wind,2))**1/2
    return fVentRoof2

#equation 16
def fVentRoof():
    if(eta_Roof>=eta_Roof_Thr):
        return eta_InsScr()*fVentRoof2() + 0.5*f_leakage()
    else:
        return eta_InsScr()*(u_th_scr *fVentRoof2()+(1-u_th_scr)*f_vent_roof_side(u_roof,u_side,a_roof,a_side)*eta_Side)+0.5*f_leakage()
    
#Giai thuat Euler
#Theo ly thuyet
# y'(t) = f(t,y(t)) <=> CO2_Air'(t) = f(t,CO2_Air(t)) <=> CO2_Air' = dx
# => CO2_Air(t0+h) = CO2_Air(t0) + h*f(t0,CO2_Air(t0)) = CO2_Air(t0) + h*dx
def Euler(CO2_Air_t0,CO2_Top_t0,h,t=dx):
    CO2_Air = CO2_Air_t0 + h*t(CO2_Air_t0,CO2_Top_t0)[0]
    CO2_Top = CO2_Top_t0 + h*t(CO2_Air_t0,CO2_Top_t0)[1]
    return [CO2_Air,CO2_Top] 


#Quan
#equation 13
def f_VentSide():
    if eta_Side >= eta_Side_Thr:
        return eta_InsScr()*f_vent_roof_side(u_roof,u_side,0,a_side) + 0.5*f_leakage()
    else:
        f2VentSide1 = f_vent_roof_side(u_roof, u_side, a_roof, a_side)
        f2VentSide2 = f_vent_roof_side(u_roof, u_side, 0, a_side)

        return eta_InsScr()*(u_th_scr*f2VentSide2 + (1 - u_th_scr)*f2VentSide1*eta_Side) + 0.5*f_leakage()
#Runge-Kutta bac 4
#gia su dx()[0] la co2_air', dx()[1] la co2_top'
def rk4(co2_air_t0, co2_top_t0, h, t, g = dx):
    #buoc nhay: (t + h - t) / h
    n = 1
    #gia tri tai (t + h) cua co2_air
    co2_air = co2_air_t0
    #gia tri tai (t + h) cua co2_top
    co2_top = co2_top_t0
    for i in range(1, n + 1):
        k1 = h * g[0](t, co2_air)
        k2 = h * g[0](t + 0.5 * h, co2_air + 0.5 * k1)
        k3 = h * g[0](t + 0.5 * h, co2_air + 0.5 * k2)
        k4 = h * g[0](t + h, co2_air + k3)
        co2_air = co2_air + (1.0 / 6.0)*(k1 + 2 * k2 + 2 * k3 + k4)
        t = t + h
    for i in range(1, n + 1):
        k1 = h * g[1](t, co2_top)
        k2 = h * g[1](t + 0.5 * h, co2_top + 0.5 * k1)
        k3 = h * g[1](t + 0.5 * h, co2_top + 0.5 * k2)
        k4 = h * g[1](t + h, co2_top + k3)
        co2_top = co2_top + (1.0 / 6.0)*(k1 + 2 * k2 + 2 * k3 + k4)
        t = t + h
    return [co2_air, co2_top]

#cong thuc 19
def hCBuf():
    if CBuf > CMaxBuf:
        hCBuf=0
    else:
        hCBuf=1
    return hCBuf

#Tri (P va J)
def P(CO2_Air):
    return (J()*(CO2_Stom(CO2_Air)-Gamma()))/(4*(CO2_Stom(CO2_Air)+2*Gamma()))

def R_func(CO2_Air):
    return P(CO2_Air) * Gamma() / CO2_Stom(CO2_Air)
#Global variable
def J():
    return (J_POT() + alpha * PAR_Can - math.sqrt(pow(J_POT() + alpha * PAR_Can, 2)- 4 * 0.7 * J_POT() * alpha * PAR_Can))/(2*theta)

#J_MAX_25_CAN
def J_MAX_25_CAN():
    return LAI * J_MAX_Leaf
#J_POT
def J_POT():
    return J_MAX_25_CAN() * math.exp(E_j * (T_CanK - T_K) / (R * T_CanK * T_K)) * (1 + math.exp((S * T_K - H) / (R * T_K))) / (1 + math.exp((S * T_CanK - H) / (R * T_CanK)))

c_gamma = 1.7
nCo2Air_Stom = 0.67
#CO2_Stom
def CO2_Stom(CO2_Air):
    return nCo2Air_Stom * CO2_Air
#Gamma
def Gamma(T_Can = t_air):
    return (J_MAX_Leaf * c_gamma * T_Can) / J_MAX_25_CAN() + 20 * c_gamma * (1 - J_MAX_Leaf / J_MAX_25_CAN())

# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    print(Euler(484,484,5))

