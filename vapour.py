import math
import copy
import pandas
import numpy
import matplotlib
import matplotlib.pyplot as plt

class Environment():
    #environmental factor
    g = 9.81
    Md = 0.0289654
    Mv = 0.018016
    rho = 101325
    M_Water = 18

    #outside parameters
    #VP_Out = 700 #(Pa) #1500
    Rh_Out = 0.88
    T_Out = 11.4 + 273.15
    v_Wind = 1

    #inside parameters
    delta_T = 1
    T_Air = 18.8 + 273.15
    T_Top = T_Air + 1
    Rh = 0.817
    VP_Out = math.exp(77.3450 + 0.0057 * T_Out - 7235 / T_Out) / T_Out**8.2 * Rh
    VP_Air = math.exp(77.3450 + 0.0057 * T_Air - 7235 / T_Air) / T_Air**8.2 * Rh
    VP_Top = math.exp(77.3450 + 0.0057 * T_Top - 7235 / T_Top) / T_Air**8.2 * Rh
    T_Can = T_Air + 1
    VP_Can = (610.78) * math.exp((T_Can - 273.15) / (T_Can - 273.15 + 238.3) * 17.2694) * Rh_Out

    #greenhouse control parameters
    Vent_Lee = float(0)/100
    Vent_Wind = float(0)/100
    U_Side = 0
    U_Roof = (Vent_Lee + Vent_Wind)/2
    eta_Side = 0
    eta_Roof = 1
    U_Blow = 0 
    U_Pad = 0
    U_Th_Scr = 95/100
    U_Vent_Forced = 0
    U_Fog = 0

    #greenhouse design parameters
    A_Flr = 1.4 * 10**4
    A_Roof = 0.3 * A_Flr
    A_Side = 0
    C_d = 0.75
    C_w = 0.09
    h_SideRoof = 0
    zeta_Ins_Scr = 1
    c_leakage = 10**(-4)
    eta_Roof_Thr = 0.9
    eta_Side_Thr = 0.9

    #height of greenhouse (m)
    h_Air = 3.8
    h_Top = 0.4

    #direct air heater parameters
    eta_HeatVap = 4.43*10**-8
    P_Blow = 0

    #pad and fan system parameters:
    phi_Pad = 0

    #fog system parameters
    phi_Fog = 0

    #thermal screen parameters
    K_Th_Scr = 0.05 * 10**(-3)

    #forced cooling system parameters
    phi_Vent_Forced = 0

    #photosynthesis related
    LAI = 2
    R = 8.314
    denta_H = 2.45*10**6
    gamma = 65.8
    c_p_air = 10**3
    rb = 275
    rs = 82

    #pad and fan
    eta_Pad = 0
    x_Pad = 0
    x_Out = 0

    #Cov,in
    A_Cov = 1.8*10**4
    T_Cov_in = T_Air + 1
    VP_Cov_in = (610.78) * math.exp((T_Cov_in - 273.15) / (T_Cov_in - 273.15 + 238.3) * 17.2694) * Rh_Out
    c_HEC_in = 1.86

    #ThScr
    T_ThScr = T_Air + 1 
    VP_ThScr = (610.78) * math.exp((T_ThScr - 273.15) / (T_ThScr - 273.15 + 238.3) * 17.2694) * Rh_Out

    #Mech
    VP_Mech = 0
    HEC_AirMech = 0

def dx(env1: Environment):
    a = (MV_CanAir(env1) + MV_PadAir(env1) + MV_FogAir(env1) + MV_BlowAir(env1) - MV_AirThSrc(env1) - MV_AirTop(env1) - MV_AirOut(env1) - MV_Airout_Pad(env1) - MV_AirMech(env1))/cap_VP_Air(env1)
    b = (MV_AirTop(env1) - MV_TopCov_in(env1) - MV_TopOut(env1))/cap_VP_Top(env1)
    return [a,b]

def MV_CanAir(env1: Environment):
    return VEC_CanAir(env1) * (env1.VP_Can - env1.VP_Air)

def MV_PadAir(env1: Environment):
    return rho_Air(env1) * env1.U_Pad * env1.phi_Pad / env1.A_Flr * (env1.eta_Pad * (env1.x_Pad - env1.x_Out) + env1.x_Out)

def MV_FogAir(env1: Environment):
    return env1.U_Fog * env1.phi_Fog / env1.A_Flr

def MV_BlowAir(env1: Environment):
    return env1.eta_HeatVap * env1.U_Blow * env1.P_Blow / env1.A_Flr

def MV_AirThSrc(env1: Environment):
    if env1.VP_Air < env1.VP_ThScr:
        return 0
    else:
        return 6.4*10**-9 * HEC_AirThSrc(env1) * (env1.VP_Air - env1.VP_ThScr)

def MV_AirTop(env1: Environment):
    return ((env1.M_Water * f_Th_Scr(env1)) / env1.R) * (env1.VP_Air / env1.T_Air - env1.VP_Top / env1.T_Top)

def MV_AirOut(env1: Environment):
    return ((env1.M_Water * (f_Vent_Side(env1) + f_Vent_Forced(env1))) / env1.R) * (env1.VP_Air / env1.T_Air - env1.VP_Out / env1.T_Out)

def MV_Airout_Pad(env1: Environment):
    return (env1.U_Pad * env1.phi_Pad / env1.A_Flr) * (env1.M_Water / env1.R) * (env1.VP_Air / env1.T_Air)

def MV_AirMech(env1: Environment):
    if (env1.VP_Air < env1.VP_Mech):
        return 0
    else:
        return 6.4*10**-9*env1.HEC_AirMech*(env1.VP_Air-env1.VP_Mech)

def MV_TopCov_in(env1: Environment):
    if(env1.VP_Top < env1.VP_Cov_in):
        return 0
    else:
        return 6.4*10**-9*HEC_TopCov_in(env1)*(env1.VP_Top-env1.VP_Cov_in)

def MV_TopOut(env1: Environment):
    return (env1.M_Water * f_Vent_Roof(env1) / env1.R) * (env1.VP_Top / env1.T_Top - env1.VP_Out / env1.T_Out)

#Calaculate capacity
def cap_VP_Air(env1: Environment):
    return (env1.M_Water * env1.h_Air) / (env1.R * env1.T_Air)

def cap_VP_Top(env1: Environment):
    return (env1.M_Water * env1.h_Top) / (env1.R * env1.T_Top)

#calculate fluxes
def f_Th_Scr(env1: Environment):
    rho_Mean_Air = (rho_Air(env1) + rho_Top(env1))/2
    return env1.U_Th_Scr * env1.K_Th_Scr * math.fabs(env1.T_Air - env1.T_Top)**(2/3) + (1-env1.U_Th_Scr) * math.sqrt(env1.g * (1- env1.U_Th_Scr) / (2* rho_Mean_Air) * math.fabs(rho_Air(env1) - rho_Top(env1)))

def f_Vent_Side(env1: Environment):
    if env1.eta_Side >= env1.eta_Side_Thr:
        return eta_Ins_Scr(env1) * diff2_f_Vent_Side(env1) + 0.5 * f_leakage(env1)
    else:
        return eta_Ins_Scr(env1) * (env1.U_Th_Scr * diff2_f_Vent_Side(env1) + (1 - env1.U_Th_Scr) * diff2_f_Vent_Roof_Side(env1) * env1.eta_Side) + 0.5 * f_leakage(env1)

def f_Vent_Forced(env1: Environment):
    return eta_Ins_Scr(env1) * env1.U_Vent_Forced * env1.phi_Vent_Forced / env1.A_Flr

def f_Vent_Roof(env1: Environment):
    if env1.eta_Roof >= env1.eta_Roof_Thr:
        return eta_Ins_Scr(env1) * diff2_f_Vent_Roof(env1) + 0.5 * f_leakage(env1)
    else:
        return eta_Ins_Scr(env1) * (env1.U_Th_Scr * diff2_f_Vent_Roof(env1) + (1 - env1.U_Th_Scr) * diff2_f_Vent_Roof_Side(env1) * env1.eta_Side) + 0.5 * f_leakage(env1)


#Other subfunctions
def VEC_CanAir(env1: Environment):
    return (2* rho_Air(env1) * env1.c_p_air * env1.LAI) / (env1.denta_H * env1.gamma * (env1.rb + env1.rs))

def HEC_AirThSrc(env1: Environment):
    return 1.7 * env1.U_Th_Scr * math.fabs(env1.T_Air - env1.T_ThScr)**0.33

def HEC_TopCov_in(env1: Environment):
    return env1.c_HEC_in * (env1.T_Top-env1.T_Cov_in)**0.33 * env1.A_Cov / env1.A_Flr

def diff2_f_Vent_Roof_Side(env1: Environment):
    T_Mean_Air = (env1.T_Air + env1.T_Out)/2
    A = env1.U_Roof**2 * env1.U_Side**2 * env1.A_Roof**2 * env1.A_Side**2 / (env1.U_Roof**2 * env1.A_Roof**2 + env1.U_Side**2 * env1.A_Side**2) if (env1.U_Roof**2 * env1.A_Roof**2 + env1.U_Side**2 * env1.A_Side**2) != 0 else 0
    return env1.C_d/env1.A_Flr * math.sqrt(2 * env1.g * env1.h_SideRoof * A *(env1.T_Air - env1.T_Out)/T_Mean_Air + ((env1.U_Roof*env1.A_Roof + env1.U_Side*env1.A_Side)/2)**2 * env1.C_w * env1.v_Wind**2)

def diff2_f_Vent_Side(env1: Environment):
    return env1.C_d * env1.U_Side * env1.A_Side * env1.v_Wind / (2 * env1.A_Flr) * math.sqrt(env1.C_w)

def diff2_f_Vent_Roof(env1: Environment):
    T_Mean_Air = (env1.T_Air + env1.T_Out)/2
    return (env1.C_d * env1.U_Roof * env1.A_Roof)/(2 * env1.A_Flr)*math.sqrt(env1.g * env1.h_SideRoof * (env1.T_Air-env1.T_Out)/(2*T_Mean_Air)+env1.C_w*pow(env1.v_Wind,2))

def eta_Ins_Scr(env1: Environment):
    return env1.zeta_Ins_Scr * (2- env1.zeta_Ins_Scr)

def f_leakage(env1: Environment):
    return env1.v_Wind * env1.c_leakage if env1.v_Wind >= 0.25 else 0.25 * env1.c_leakage


def rho(Md,Mv,p, R,T,Rh):
    return (p* Md + Rh * 610.78 * 10**(7.5 * (T-273.15) / ((T-273.15) + 237.3)) * (Mv - Md)) / (R*(T-273.15))

def rho_Air(env1: Environment):
    return rho(env1.Md,env1.Mv,env1.rho,env1.R,env1.T_Air,env1.Rh)

def rho_Top(env1: Environment):
    return rho(env1.Md,env1.Mv,env1.rho,env1.R,env1.T_Top,env1.Rh)


#Euler
def Euler(env1: Environment,h,t=dx):
    VP_Air = env1.VP_Air + h*t(env1)[0]
    VP_Top = env1.VP_Top + h*t(env1)[1]
    return round(VP_Air,2),round(VP_Top,2)

#rk4
def rk4(env1: Environment,h,g = dx):
    #buoc nhay
    #gia tri tai (t + h) cua VPair
    #gia tri tai (t + h) cua VPtop
    k1 = round(g(env1)[0],2)
    l1 = round(g(env1)[1],2)
    env2 = copy.deepcopy(env1)
    env2.VP_Air = env1.VP_Air + 0.5 *h*k1
    env2.VP_Top = env1.VP_Top + 0.5 *h*l1
    k2 = round(g(env2)[0],2)
    l2 = round(g(env2)[1],2)
    env2.VP_Air = env1.VP_Air + 0.5 *h*k2
    env2.VP_Top = env1.VP_Top + 0.5 *h*l2
    k3 = round(g(env2)[0],2)
    l3 = round(g(env2)[1],2)
    env2.VP_Air = env1.VP_Air + h*k3
    env2.VP_Top = env1.VP_Top + h*l3
    k4 = round(g(env2)[0],2)
    l4 = round(g(env2)[1],2)
    VPair = env1.VP_Air + (h / 6.0)*(k1 + 2 * k2 + 2 * k3 + k4)
    VPtop = env1.VP_Top + (h / 6.0)*(l1 + 2 * l2 + 2 * l3 + l4)
    return round(VPair,2), round(VPtop,2)

def update(env1: Environment, row):
    if not pandas.isnull(row["Temp"]):
        env1.T_Air = row["Temp"] + 273.15
    if not pandas.isnull(row["Rh"]):
        env1.Rh = row["Rh"] /100
    if not pandas.isnull(row["VentLee"]):    
        env1.Vent_Lee = row["VentLee"] /100
    if not pandas.isnull(row["VentWind"]):
        env1.Vent_Wind = row["VentWind"] /100
    if not pandas.isnull(row["EnergyCurtain"]):
        env1.U_Th_Scr = row["EnergyCurtain"] /100
    if not pandas.isnull(row["OutsideTemp"]):
        env1.T_Out = row["OutsideTemp"] + 273.15
    if not pandas.isnull(row["OutsideRh"]):
        env1.Rh_Out = row["OutsideRh"] / 100
    if not pandas.isnull(row["WindSpeed"]):
        env1.v_Wind = row["WindSpeed"]
    env1.U_Roof = (env1.Vent_Lee + env1.Vent_Wind)/2
    env1.VP_Out= math.exp(77.3450 + 0.0057 * env1.T_Out - 7235 / env1.T_Out) / env1.T_Out**8.2 * env1.Rh
    env1.T_Can = env1.T_Air + 1
    env1.VP_Can = (610.78) * math.exp((env1.T_Can - 273.15) / (env1.T_Can - 273.15 + 238.3) * 17.2694) * env1.Rh_Out
    env1.T_Top = env1.T_Air + 1
    env1.T_Cov_in = env1.T_Air + 1
    env1.VP_Cov_in = (610.78) * math.exp((env1.T_Cov_in - 273.15) / (env1.T_Cov_in - 273.15 + 238.3) * 17.2694) * env1.Rh_Out
    env1.T_ThScr = env1.T_Air + 1 
    env1.VP_ThScr = (610.78) * math.exp((env1.T_ThScr - 273.15) / (env1.T_ThScr - 273.15 + 238.3) * 17.2694) * env1.Rh_Out

def run_Euler(env1: Environment, data_set,comparison_table,output_name):
    for i in range(10):         #for each 5 minutes
        row = data_set.iloc[i]
        #print("New environment:",row)
        update(env1,row)
        if i ==0:
            env1.VP_Air = math.exp(77.3450 + 0.0057 * env1.T_Air - 7235 / env1.T_Air) / env1.T_Air**8.2 * env1.Rh
            env1.VP_Top = math.exp(77.3450 + 0.0057 * env1.T_Top - 7235 / env1.T_Top) / env1.T_Air**8.2 * env1.Rh
        vp_air_real = math.exp(77.3450 + 0.0057 * env1.T_Air - 7235 / env1.T_Air) / env1.T_Air**8.2 * env1.Rh
        new_row = {"Timestamp":row["Timestamp"],"VP_Air_real":vp_air_real,"VP_Air_calc":env1.VP_Air}
        comparison_table = comparison_table.append(new_row, ignore_index = True)
        for j in range(5*60):    #for each seconds
            #print("Step:",i*5*60+j)
            #print("dx:",dx(env1))
            prediction = Euler(env1,1)
            #print("result =",prediction)
            env1.VP_Air = prediction[0]
            env1.VP_Top = prediction[1]
    print(comparison_table)
    comparison_table.set_index("Timestamp").plot(figsize=(10,5), grid=True, title= "euler_plot (Pa)")
    plt.show()
    #comparison_table.to_excel(output_name)

def run_rk4(env1: Environment, data_set,comparison_table,output_name):
    for i in range(10):         #for each 5 minutes
        row = data_set.iloc[i]
        #print("New environment:",row)
        update(env1,row)
        if i ==0:
            env1.VP_Air = math.exp(77.3450 + 0.0057 * env1.T_Air - 7235 / env1.T_Air) / env1.T_Air**8.2 * env1.Rh
            env1.VP_Top = math.exp(77.3450 + 0.0057 * env1.T_Top - 7235 / env1.T_Top) / env1.T_Air**8.2 * env1.Rh
        vp_air_real = math.exp(77.3450 + 0.0057 * env1.T_Air - 7235 / env1.T_Air) / env1.T_Air**8.2 * env1.Rh
        new_row = {"Timestamp":row["Timestamp"],"VP_Air_real":vp_air_real,"VP_Air_calc":env1.VP_Air}
        comparison_table = comparison_table.append(new_row, ignore_index = True)
        for j in range(5*60):    #for each seconds
            #print("Step:",i*5*60+j)
            #print("dx:",dx(env1))
            prediction = rk4(env1,1)
            #print("result =",prediction)
            env1.VP_Air = prediction[0]
            env1.VP_Top = prediction[1]
    print(comparison_table)
    comparison_table.set_index("Timestamp").plot(figsize=(10,5), grid=True, title= "rk4_plot (Pa)")
    plt.show()
    #comparison_table.to_excel(output_name)


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    data_set = pandas.read_csv("environment.csv",usecols=["Timestamp","Temp","Rh","VentLee","VentWind","EnergyCurtain","OutsideTemp", "OutsideRh","WindSpeed"])
    data_set = data_set.iloc[0:10,]
    data_set = data_set.reset_index(drop = True)
    print(data_set)
    comparison_table = pandas.DataFrame(columns = ["Timestamp","VP_Air_real","VP_Air_calc"])
    env1 = Environment()
    print(dx(env1))
    #run_Euler(env1,data_set,comparison_table,"output_Euler.xlsx")
    run_rk4(env1,data_set,comparison_table,"output_rk4.xlsx")
