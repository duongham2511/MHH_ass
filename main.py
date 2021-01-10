# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.

import math

vent_lee = 0
vent_wind = 0

CO2_out = 412.89
cap_CO2_Air = 3.8
cap_CO2_Top = 0.4

CO2_dosage = 0
P_Blow = 0
n_HeatCO2=0.057
U_Blow= 0

A_Flr = 1.4 * 10**4
U_Pad = 0
phi_Pad = 0

CO2_dosage = 0.14
O_ExtCO2=7.2*10**4
U_ExtCO2= 0

t_air = 19.9
t_top = 19.9
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
c_gamma = 1.7
nCo2Air_Stom = 0.67

def dx(t,CO2_Air,CO2_Top):
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
    CO2_Air = CO2_Air_t0 + h*t(0,CO2_Air_t0,CO2_Top_t0)[0]
    CO2_Top = CO2_Top_t0 + h*t(0,CO2_Air_t0,CO2_Top_t0)[1]
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
def rk4(co2_air_t0, co2_top_t0, h, g = dx):
    #buoc nhay
    n = 1

    #gia tri tai (t + h) cua co2_air
    co2_air = co2_air_t0
    #gia tri tai (t + h) cua co2_top
    co2_top = co2_top_t0
    for i in range(1, n + 1):
        k1 = g(0,co2_air, co2_top)[0]
        l1 = g(0,co2_air, co2_top)[1]
        k2 = g(0,co2_air + 0.5 *h*k1, co2_top + 0.5 *h* l1)[0]
        l2 = g(0,co2_air + 0.5 *h*k1, co2_top + 0.5 *h* l1)[1]
        k3 = g(0,co2_air + 0.5 * h * k2, co2_top + 0.5 * h * l2)[0]
        l3 = g(0,co2_air + 0.5 * h * k2, co2_top + 0.5 * h * l2)[1]
        k4 = g(0,co2_air + h * k3, co2_top + l3 * h)[0]
        l4 = g(0,co2_air + h * k3, co2_top + l3 * h)[1]
        co2_air = co2_air + (h / 6.0)*(k1 + 2 * k2 + 2 * k3 + k4)
        co2_top = co2_top + (h / 6.0)*(l1 + 2 * l2 + 2 * l3 + l4)
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
#CO2_Stom
def CO2_Stom(CO2_Air):
    return nCo2Air_Stom * CO2_Air
#Gamma
def Gamma(T_Can = t_air):
    return (J_MAX_Leaf * c_gamma * T_Can) / J_MAX_25_CAN() + 20 * c_gamma * (1 - J_MAX_Leaf / J_MAX_25_CAN())


# Câu 5 -------------------------------------------------------------------------------------------------------------------------------------
#Global variable
denta_H = 2.45*10**6
gamma = 65.8
eta_HeatVap = 4.43*10**-8
M_Water = 18
U_Fog = 0
phi_Fog = 0
VP_Out = 1.5
h_Air = 3.8
VP_Can = VP_ThScr = 1 #
HEC_AirThSrc = 0 # Do u_th_scr = 0
c_p_air = 10**3
rb = 275
rs = 82
A_Cov = 1.8*10**4
t_cov_in = 17 #
VP_TopCov = 1.2 #
VP_in = 1 #
c_HEC_in = 1.86
h_Top = 3 #
def MV_BlowAir():
    return eta_HeatVap*U_Blow*P_Blow/A_Flr
def MV_PadAir():
    return p_air*U_Pad*phi_Pad/A_Flr
def MV_Airout_Pad(VP_Air):
    return (U_Pad*phi_Pad/A_Flr)*(M_Water/R)*(VP_Air/(t_air+273.15))
def MV_FogAir():
    return U_Fog*phi_Fog/A_Flr
def MV_AirTop(VP_Air,VP_Top):
    return ((M_Water*F_Th_Scr())/R)*((VP_Air/(t_air+273.15))-VP_Top/(t_top+273.15))
def MV_AirOut(VP_Air):
    return ((M_Water*(f_VentSide()+f_VentForced()))/R)*((VP_Air/(t_air+273.15))-VP_Out/(t_out+273.15))
def cap_VP_Air():
    return (M_Water*h_Air)/(R*(t_air+273.15))
def MV_AirThSrc(VP_Air):
    if VP_Air < VP_ThScr:
        return 0
    else:
        return 6.4*10**-9*HEC_AirThSrc*(VP_Air-VP_ThScr)
def VEC_CanAir():
    return (2*p_air*c_p_air*LAI)/(denta_H*gamma*(rb+rs))
def MV_CanAir(VP_Air):
    return VEC_CanAir()*(VP_Can-VP_Air)
def MV_TopOut(VP_Top):
    return (M_Water*fVentRoof()/R)*((VP_Top/(t_top+273.15))-(VP_Out/(t_out+273.15)))
def HEC_TopCov_in():
    return c_HEC_in*(t_top-t_cov_in)**0.33*A_Cov/A_Flr
def MV_TopCov_in():
    if(VP_TopCov<VP_in):
        return 0
    else:
        return 6.4*10**-9*HEC_TopCov_in()*(VP_TopCov-VP_in)
def cap_VP_Top():
    return (M_Water*h_Top)/(R*(t_top+273.15))
def dx2(VP_Air,VP_Top):
    a = (MV_CanAir(VP_Air) + MV_PadAir() + MV_FogAir() + MV_BlowAir() - MV_AirThSrc(VP_Air) - MV_AirTop(VP_Air,VP_Top) - MV_AirOut(VP_Air) - MV_Airout_Pad(VP_Air))/cap_VP_Air()
    b = (MV_AirTop(VP_Air,VP_Top) - MV_TopCov_in() - MV_TopOut(VP_Top))/cap_VP_Top()
    return [a,b]

# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    print(dx(0,484,484))
    print(Euler(484,484,5*60))
    print(rk4(484,484,5*60))

