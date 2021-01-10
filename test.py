import math
import copy
import pandas
import numpy
import matplotlib
import matplotlib.pyplot as plt

class Environment:

    #environmental factor
    g = 9.81
    M_CH20 = 30 * 10**(-3)
    T_25 = 298.15
    Md = 0.0289654
    Mv = 0.018016
    rho = 101325

    #conversion factor
    eta_mg_ppm = 0.554 #mg/m^3 to ppm

    #outside parameters
    CO2_Out = 668 #mg/m^3
    T_Out = 11.6 + 273.15
    v_Wind = 1

    #inside parameters
    delta_T = 1
    CO2_Air = 1296.01
    CO2_Top = 1296.01
    T_Air = 18.9 + 273.15
    T_Top = T_Air + 1
    T_Can = T_Top
    Rh = 81.6 /100

    #greenhouse control parameters
    Vent_Lee = float(0)/100
    Vent_Wind = float(0)/100
    U_Side = 0
    U_Roof = (Vent_Lee + Vent_Wind)/2
    eta_Side = 0
    eta_Roof = 1
    U_Blow = 0 
    U_Ext_CO2 = 0
    U_Pad = 0
    U_Th_Scr = 95/100
    U_Vent_Forced = 0

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

    #capacity of greenhouse (m)
    cap_CO2_Air = 3.8
    cap_CO2_Top = 0.4

    #direct air heater parameters
    eta_Heat_CO2 = 0.057
    P_Blow = 0

    #external CO2 source parameters
    phi_Ext_CO2 = 7.2 * 10**4

    #pad and fan system parameters:
    phi_Pad = 0

    #thermal screen parameters
    K_Th_Scr = 0.05 * 10**(-3)

    #forced cooling system parameters
    phi_Vent_Forced = 0

    #photosynthesis related
    c_Buf = 0
    c_Max_Buf = 20 * 10**3
    eta_CO2_Stom = 0.67
    LAI = 2
    J_MAX_25_Leaf = 210
    c_gamma = 1.7
    alpha = 0.385
    PAR_Can = 100
    Theta = 0.7
    E_j = 37 * 10**3
    R = 8.314
    S = 710
    H = 22 * 10**4


def dx(env1: Environment):
    dot_CO2_Air = round((MC_BlowAir(env1) + MC_ExtAir(env1) + MC_PadAir(env1) - MC_AirCan(env1) -MC_AirTop(env1) - MC_AirOut(env1))/float(env1.cap_CO2_Air),2)
    dot_CO2_Top = round((MC_AirTop(env1) - MC_TopOut(env1))/float(env1.cap_CO2_Top),2)
    return dot_CO2_Air,dot_CO2_Top

def MC_BlowAir(env1: Environment):
    return env1.eta_Heat_CO2 * env1.U_Blow * env1.P_Blow / env1.A_Flr

def MC_ExtAir(env1: Environment):
    return env1.U_Ext_CO2 * env1.phi_Ext_CO2 / env1.A_Flr

def MC_PadAir(env1: Environment):
    return env1.U_Pad * env1.phi_Pad / env1.A_Flr * (env1.CO2_Out - env1.CO2_Air)

def MC_AirCan(env1: Environment):
    return env1.M_CH20 * h_c_Buf(env1) * (P(env1) * 0.99)

def MC_AirTop(env1: Environment):
    return f_Th_Scr(env1) * (env1.CO2_Air - env1.CO2_Top)

def MC_AirOut(env1: Environment):
    return (f_Vent_Side(env1) + f_Vent_Forced(env1)) * (env1.CO2_Air - env1.CO2_Out)

def MC_TopOut(env1: Environment):
    return f_Vent_Roof(env1) * (env1.CO2_Top - env1.CO2_Out)

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

#calculate photosynthesis variable
def h_c_Buf(env1: Environment):
    return 0 if env1.c_Buf > env1.c_Max_Buf else 1

def P(env1: Environment):
    return J(env1) * (CO2_Stom(env1) - Gamma(env1)) / (4 * (CO2_Stom(env1) - 2 * Gamma(env1)))

def J(env1: Environment):
    return (J_POT(env1) + env1.alpha * env1.PAR_Can - math.sqrt((J_POT(env1) + env1.alpha * env1.PAR_Can)**2 - 4 * env1.Theta * J_POT(env1) * env1.alpha * env1.PAR_Can)) / (2 * env1.Theta)

def J_POT(env1: Environment):
    return J_MAX_25_Can(env1) * math.exp(env1.E_j * (env1.T_Can - env1.T_25) / (env1.R * env1.T_Can * env1.T_25)) * (1 + math.exp((env1.S * env1.T_25 - env1.H) / (env1.R * env1.T_25))) / (1 + math.exp((env1.S * env1.T_Can - env1.H) / (env1.R * env1.T_Can)))

def CO2_Stom(env1: Environment):
    return env1.eta_CO2_Stom * env1.CO2_Air

def Gamma(env1: Environment):
    return env1.J_MAX_25_Leaf / J_MAX_25_Can(env1) * env1.c_gamma * env1.T_Can + 25 * env1.c_gamma * (1 - env1.J_MAX_25_Leaf/J_MAX_25_Can(env1))

def J_MAX_25_Can(env1: Environment):
    return env1.LAI * env1.J_MAX_25_Leaf

#other subfunctions
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
    CO2_Air = env1.CO2_Air + h*t(env1)[0]
    CO2_Top = env1.CO2_Top + h*t(env1)[1]
    return round(CO2_Air,2),round(CO2_Top,2)

#rk4
def rk4(env1: Environment,h,g = dx):
    #buoc nhay
    #gia tri tai (t + h) cua co2_air
    #gia tri tai (t + h) cua co2_top
    k1 = round(g(env1)[0],2)
    l1 = round(g(env1)[1],2)
    env2 = copy.deepcopy(env1)
    env2.CO2_Air = env1.CO2_Air + 0.5 *h*k1
    env2.CO2_Top = env1.CO2_Top + 0.5 *h*l1
    k2 = round(g(env2)[0],2)
    l2 = round(g(env2)[1],2)
    env2.CO2_Air = env1.CO2_Air + 0.5 *h*k2
    env2.CO2_Top = env1.CO2_Top + 0.5 *h*l2
    k3 = round(g(env2)[0],2)
    l3 = round(g(env2)[1],2)
    env2.CO2_Air = env1.CO2_Air + h*k3
    env2.CO2_Top = env1.CO2_Top + h*l3
    k4 = round(g(env2)[0],2)
    l4 = round(g(env2)[1],2)
    co2_air = env1.CO2_Air + (h / 6.0)*(k1 + 2 * k2 + 2 * k3 + k4)
    co2_top = env1.CO2_Top + (h / 6.0)*(l1 + 2 * l2 + 2 * l3 + l4)
    return round(co2_air), round(co2_top,2)

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
    if not pandas.isnull(row["Co2ActuationRegulation"]):
        env1.U_Ext_CO2 = row["Co2ActuationRegulation"]
    if not pandas.isnull(row["OutsideTemp"]):
        env1.T_Out = row["OutsideTemp"] + 273.15
    if not pandas.isnull(row["WindSpeed"]):
        env1.v_Wind = row["WindSpeed"]

def run_Euler(env1: Environment, data_set,comparison_table,output_name):
    for i in range(5):         #for each 5 minutes
        row = data_set.iloc[i]
        #print("New environment:",row)
        update(env1,row)
        if i ==0:
            env1.CO2_Air = round(row["Co2"] / env1.eta_mg_ppm,2)
            env1.CO2_Top = round(row["Co2"] / env1.eta_mg_ppm,2)
        new_row = {"Timestamp":row["Timestamp"],"CO2_Air_real":row["Co2"],"CO2_Air_calc":env1.CO2_Air*env1.eta_mg_ppm}
        comparison_table = comparison_table.append(new_row, ignore_index = True)
        for j in range(5*60):    #for each seconds
            #print("Step:",i*5*60+j)
            #print("dx:",dx(env1))
            prediction = Euler(env1,1)
            #print("result =",prediction)
            env1.CO2_Air = prediction[0]
            env1.CO2_Top = prediction[1]
    print(comparison_table)
    comparison_table.set_index("Timestamp").plot(figsize=(10,5), grid=True, title= "Euler_plot (ppm)")
    plt.show()
    comparison_table.to_excel(output_name)

def run_rk4(env1: Environment, data_set,comparison_table,output_name):
    for i in range(5):         #for each 5 minutes
        row = data_set.iloc[i]
        #print("New environment:",row)
        update(env1,row)
        if i ==0:
            env1.CO2_Air = round(row["Co2"] / env1.eta_mg_ppm,2)
            env1.CO2_Top = round(row["Co2"] / env1.eta_mg_ppm,2)
        new_row = {"Timestamp":row["Timestamp"],"CO2_Air_real":row["Co2"],"CO2_Air_calc":env1.CO2_Air*env1.eta_mg_ppm}
        comparison_table = comparison_table.append(new_row, ignore_index = True)
        for j in range(5*60):    #for each seconds
            #print("Step:",i*5*60+j)
            #print("dx:",dx(env1))
            prediction = rk4(env1,1)
            #print("result =",prediction)
            env1.CO2_Air = prediction[0]
            env1.CO2_Top = prediction[1]
    print(comparison_table)
    comparison_table.set_index("Timestamp").plot(figsize=(10,5), grid=True, title= "rk4_plot (ppm)")
    plt.show()
    comparison_table.to_excel(output_name)


if __name__ == '__main__':
    data_set = pandas.read_csv("environment.csv",usecols=["Timestamp","Temp","Rh","Co2","VentLee","VentWind","EnergyCurtain","Co2ActuationRegulation","OutsideTemp","WindSpeed"])
    data_set = data_set.iloc[216:221,]
    data_set = data_set.reset_index(drop = True)
    print(data_set)
    comparison_table = pandas.DataFrame(columns = ["Timestamp","CO2_Air_real","CO2_Air_calc"])
    env1 = Environment()

    run_Euler(env1,data_set,comparison_table,"output_Euler.xlsx")
    run_rk4(env1,data_set,comparison_table,"output_rk4.xlsx")
