# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.

import math


def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press Ctrl+F8 to toggle the breakpoint.


#equation 7
def F_Th_Scr(t_air,t_top,p_air,p_top,u_th_scr,k_th_scr = 0.03*10**(-3),g = 9.81):
    return u_th_scr * k_th_scr * math.pow(math.fabs(t_air - t_top),2/3) + (1-u_th_scr) * math.sqrt(g*(1-u_th_scr)*math.fabs(p_air-p_top)/(2 * (p_air+p_top)/2))

#equation 10
def f_vent_roof_side(c_d,c_w,a_flr,u_roof,u_side,a_roof,a_side,h_sideroof,g,t_air,t_out,v_wind):
    t_mean_air = (t_air + t_out)/2
    return c_d/a_flr * math.sqrt(u_roof**2 * u_side**2 * a_roof**2 * a_side**2 / (u_roof**2 * a_roof**2 + u_side**2 * a_side**2) * 2*g*h_sideroof*(t_air - t_out)/t_mean_air + ((u_roof*a_roof + u_side*a_side)/2)**2 * c_w * v_wind**2)


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    test1 = f_th_scr(0.8,0.2,9.8,37,38,0.9,1.0)
    print(test1)

# See PyCharm help at https://www.jetbrains.com/help/pycharm/

#Tri
#equation 3
def MC_BlowAir(P_Blow,n_HeatCO2=0.057,U_Blow=0.5,A_Flr=1.3*10**4):
    return (n_HeatCO2*U_Blow*P_Blow)/A_Flr

#equation 4
def MC_ExtAir(O_ExtCO2=7.2*10**4,U_ExtCO2=0.5,A_Flr=1.3*10**4):
    return (O_ExtCO2*U_ExtCO2)/A_Flr

#equation 5
def MC_PadAir(CO2_out,CO2_Air,U_Pad,phi_Pad,A_Flr):
    return u_Pad*phi_Pad/A_Flr*(CO2_out-CO2_Air)

#equation 15
def MC_TopOut(fVentRoof,CO2_out,CO2_top):
    return fVentRoof*(CO2_top-CO2_out)

#equation 16
def fVentRoof(eta_InsScr,fVentRoof2,f_leakage,eta_Roof,fVentRoofSide,eta_Side,U_ThScr=0.3,eta_Roof_Thr=0.9):
    if(eta_Roof>=eta_Roof_Thr):
        return eta_InsScr*fVentRoof2+0.5*f_leakage
    else:
        return eta_InsScr*(U_ThScr*fVentRoof2+(1-U_ThScr)*fVentRoofSide*eta_Side)+0.5*f_leakage
    
#Giai thuat Euler
#hàm dx mẫu lan luot chua gia tri cua CO2_Air' va CO2_Top'
def dx():
    return [3,6]
#Theo ly thuyet
# y'(t) = f(t,y(t)) <=> CO2_Air'(t) = f(t,CO2_Air(t)) <=> CO2_Air' = dx
# => CO2_Air(t0+h) = CO2_Air(t0) + h*f(t0,CO2_Air(t0)) = CO2_Air(t0) + h*dx
def Euler(CO2_Air_t0,CO2_Top_t0,h,t=dx()):
    CO2_Air = CO2_Air_t0 + h*t[0]
    CO2_Top = CO2_Top_t0 + h*t[1]
    return [CO2_Air,CO2_Top] 


#Hien
#equation 11
def eta_InsScr(zeta_InsScr = 1):
    return zeta_InsScr * (2 - zeta_InsScr)

#equation 12
def f_leakage(v_wind = 2.4, c_leakage = 10**-4):
    if (v_wind < 0.25):
        return 0.25 * c_leakage
    else:
        return v_wind * c_leakage

#equation 14
def f_VentForced(phi_VentForced, U_VentForced = 0.5, A_Flr = 1.3 * 10**4):
    return (eta_InsScr(1) * U_VentForced * phi_VentForced) / A_Flr

# cong thuc 17
def fVentRoof2(Cd,URoof,ARoof,AFlr,g,hRoof,TAir,TOut,TMeanAir,Cw,vWind):
    fVentRoof2 = (Cd*URoof*ARoof)/(2*AFlr)*(g*hRoof*(TAir-TOut)/(2*TMeanAir)+Cw*pow(vWind,2))**1/2
    return fVentRoof2

#cong thuc 24
def fT(e,Hd,R,T0,Hd,S,T):
    fT= (1+pow(e,-Hd/R*(1/T0-1/(Hd/S))))/(1+pow(e,-Hd/R*(1/T-1/(Hd/S))))
    return fT

#cong thuc 27
def L(L0,K,e,LAI,m):
    L= L0*(1-(K*pow(e,-K*LAI))/(1-m))
    return L

#Quan
#equation 6
def MC_AirTop(F_Th_Scr,co2_air,co2_top):
    return F_Th_Scr*(co2_air - co2_top)
#equation 13
def f_VentSide(eta_InsScr,eta_Side, eta_Side_Thr, f_leakage, U_ThScr, c_d,c_w,a_flr,u_roof,u_side,a_roof,a_side,h_sideroof,g,t_air,t_out,v_wind):
    if eta_Side >= eta_Side_Thr:
        return eta_InsScr*f_vent_roof_side(c_d,c_w,a_flr,u_roof,u_side,0,a_side,h_sideroof,g,t_air,t_out,v_wind) + 0.5*f_leakage
    else:
        f2VentSide1 = f_vent_roof_side(c_d, c_w, a_flr, u_roof, u_side, a_roof, a_side, h_sideroof, g, t_air, t_out, v_wind)
        f2VentSide2 = f_vent_roof_side(c_d, c_w, a_flr, u_roof, u_side, 0, a_side, h_sideroof, g, t_air, t_out, v_wind)

        return eta_InsScr*(U_ThScr*f2VentSide2 + (1 - U_ThScr)*f2VentSide1*eta_Side) + 0.5*f_leakage
#equation 9
def MC_AirOut(f_VentSide, f_VentForced, co2_air, co2_out):
    return (f_VentSide + f_VentForced) * (co2_air - co2_out)

#cong thuc 19
def hCBuf(CBuf,CMaxBuf):
    if CBuf > CMaxBuf:
    hCBuf = 0
    else:
    hCBuf=1
    return hCBuf
