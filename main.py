# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.

import math


def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press Ctrl+F8 to toggle the breakpoint.


#equation 7
def f_th_scr(u_th_scr,k_th_scr,g,t_air,t_top,p_air,p_top):
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
