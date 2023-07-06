import os
import csv
import pandas as pd
import numpy as np
import math


## change current directory to the directory of the script
os.chdir(os.path.dirname(os.path.abspath(__file__)))

# Material properties
E = 70747.95
Sult = 530
Syield = 490
nu = 0.35
PI = 3.14159265359

# Panel Properties
t = 4
b = 200
a = 600

# Load cases
S_LC1 = 2.043
S_LC2 = 2.043
S_LC3 = 1.021

# Shell (Skin Areas) Properties
w_s_x = 200
w_s_y = 200
t_s = 4
A_x = w_s_x * t_s
A_y = w_s_y * t_s

# Stringer Specs
Sheight= 40
Sflange= 70
Sflathick = 3
Sweb = Sheight - Sflathick
Swebthick = 2
A_stringer = (Sweb*Swebthick) + (Sflange*Sflathick)
# print(A_stringer)

# Pitch Specs
# A pitch is basically the distance between two stringers
wel_x_leftpitch = 200
wel_x_rightpitch = 200
A_leftpitch = wel_x_leftpitch/2 * t_s
A_rightpitch = wel_x_rightpitch/2 * t_s

k_biax = pd.DataFrame(columns=['k_biax'])
k_shear = pd.DataFrame(columns=['k_shear'])

def find_lambda():  # ok
    
    L =  600 # Depth of the stringer
    c =1 # Hinged
    
    # Order: 
    # Skin, Flange, Web
    A_skin = 200*4
    A_flange = 70*3
    A_web = 37*2
    # Z is defined from the top of the flange and positive is downwards
    z = [-2,1.5,21.5] # 3+37/2 = 21.5
    A = [A_skin,A_flange,A_web]
    sum_A =0 
    sum_A_z=0
    for i in range(len(z)):
        sum_A_z+=A[i]*z[i]
        sum_A+=A[i]
    z_bar = sum_A_z/sum_A    
    I_y_skin = 200*4**3/12
    I_y_flange = 70*3**3/12
    I_y_web = 2*37**3/12
    z_new=[0] * 3
    for i in range(len(z)):
        z_new[i] = z_bar - z[i]
    
    I_steiner_skin = I_y_skin + A_skin*z_new[0]**2
    I_steiner_flange = I_y_flange + A_flange*z_new[1]**2
    I_steiner_web = I_y_web + A_web*z_new[2]**2
    I_total = I_steiner_skin + I_steiner_flange + I_steiner_web
    r = math.sqrt(I_total/sum_A) 
    lambda_stress = c*L/r
    
    return lambda_stress
    
def column_buckling_cr_euler(lamda_stress): # ok
    sigma_euler = (np.pi**2*E)/(lamda_stress**2)
    return sigma_euler

def choose_alpha(xi): # ok
    if (xi>= 0.4 and xi<=1.095):
        alpha = 1.4-0.628*xi
    elif (xi>1.095 and xi<=1.633):
        alpha =0.78/xi
    else:
        alpha = 0.69/pow(xi,0.75)
    return alpha

def column_buckling_cr_crip(): # ok
    r_str = 0
    b11 = Sflange/2 - (Swebthick/2) * (0.25*(Swebthick/Sflathick)-0.2*(r_str**2/(Swebthick*Sflathick))) 
    b12 = Sheight - (Sflathick/2) * (2-0.5*(Swebthick/Sflathick)-0.2*(r_str**2/(Swebthick*Sflathick))) 
    Ki = 0.41 # Supported at one side
    x11 = (b11/Sflathick) * math.sqrt(Syield/(Ki*E))
    x12 = (b12/Swebthick) * math.sqrt(Syield/(Ki*E))
    alpha11 = choose_alpha(x11)
    alpha12 = choose_alpha(x12)
    sigma_crip1 = alpha11 * Syield
    sigma_crip2 = alpha12 * Syield
    # getting the average over the three elements (two halves of the flange and one web)
    sigma = (2*(sigma_crip1*b11*Sflathick)+sigma_crip2*b12*Swebthick)/(2*b11*Sflathick+b12*Swebthick)
    # note that the equation in 6_AS_WS22_Structure_Design_Part4_filled the slide 39 is incorrect, because if smaller than Syield, the sigma_crip_i has to be inside the SIGMA symbol
    if sigma>=Syield:
        Fcrip = Syield * (2*b11*Sflathick+b12*Swebthick)
    else: 
        Fcrip = (2*sigma_crip1*b11*Sflathick+sigma_crip2*b12*Swebthick)
    sigma_crip = Fcrip / (2*b11*Sflathick+b12*Swebthick)
    return sigma_crip
    
def column_buckling_cr_crip_web(): # ok
    r_str = 0
    b12 = Sheight - (Sflathick/2) * (2-0.5*(Swebthick/Sflathick)-0.2*(r_str**2/(Swebthick*Sflathick))) 
    Ki = 0.41 # Supported at one side
    x12 = (b12/Swebthick) * math.sqrt(Syield/(Ki*E))
    alpha12 = choose_alpha(x12)
    sigma_crip2 = alpha12 * Syield
    # getting the average over the three elements (two halves of the flange and one web)
    sigma = (sigma_crip2*b12*Swebthick)/(b12*Swebthick)
    if sigma>=Syield:
        Fcrip = Syield * (b12*Swebthick)
    else: 
        Fcrip = (sigma_crip2*b12*Swebthick)
    sigma_crip = Fcrip / (b12*Swebthick)
    return sigma_crip

def column_buckling_cr_ej(): # ok
    lamda_euler_johnson = math.sqrt((2*(np.pi)**2*E)/Syield)
    lamda = min(find_lambda(),lamda_euler_johnson)
    sigma = Syield - 1/E * (Syield/(2*np.pi))**2 * lamda**2
    return sigma
    
def column_buckling_RF(combined_stringer_stress):
    lamda_stress = find_lambda()
    sigma_cr = min(column_buckling_cr_euler(lamda_stress),column_buckling_cr_ej(),column_buckling_cr_crip())
    RF = sigma_cr/combined_stringer_stress
    return RF
    
    
def optimize_K_Shear():
        alpha = a/b
        if alpha<1: 
            K = 4 + 5.34/(alpha**2)
        else: 
            K = 5.34+4/(alpha**2)
        return K

## Critical Shear in buckling
def sigma_xy_crit(): 
    sigma_xy_crit = pd.DataFrame(columns=['sigma_xy_crit'])
    for i in range(10):
        k = optimize_K_Shear()
        k_shear.loc[i] = [k]
        sigma_xy_crit.loc[i] = [(E*PI**2)/(12*(1-nu**2)) * (t/b)**2 * k ]
    return(sigma_xy_crit)
    
def optimize_K(sigma_y,sigma_x):
    n = 1
    m = 1
    alpha = a/b
    beta =sigma_y/sigma_x
    df = pd.DataFrame(columns=['K_sigma'])
    for m in range(50):
        df.loc[m] = [(m**2+n**2*alpha**2)**2/(alpha**2*(m**2+beta*n**2*alpha**2))]
    # find the minimum value of the column in df
    min_K = df['K_sigma'].min()
    # find the minimum index of that minimum K
    min_m = df['K_sigma'].idxmin()
    return min_K

def sigma_x_crit(S_YY,S_XX): 
    sigma_x_crit = pd.DataFrame(columns=['sigma_x_crit'])
    for i in range(10):
        k = optimize_K(S_YY.loc[i],S_XX.loc[i])
        k_biax.loc[i] = [k]
        sigma_x_crit.loc[i] = [(E*(np.pi)**2)/(12*(1-nu**2)) * (t/b)**2 * k]
    return sigma_x_crit


def get_combined_RF(RF_b,RF_s):
    return (1/RF_b) + (1/RF_s)**2

# a method Divide two dataframes by each other
def get_RF(df1,df2):
    df = pd.DataFrame(columns=['RF'])
    for i in range(10):
        df.loc[i] = [df1.loc[i]/df2.loc[i]]
    dfx = pd.DataFrame(columns=['RF'])
    for i in range(10):
        dfx.loc[i] = (df.loc[i])[0][0]
    return dfx

def avg_stress(df,A):
    # Calculate the sum of every three consecutive rows
    # [-2] or [4] is the last column in the txt files extracted from HW
    # those columns are the ones that has stress values in them
    
    stress = df[df.columns[4]]*A
    avg = (stress.rolling(window=3, min_periods=1).sum())/(3*A)
    # Exclude non-multiple of 3 rows from the resulting object
    avg = avg[df.index % 3 == 2]
    # Reset index starting from 1
    avg = avg.reset_index(drop=True)
    return avg

def combined_stress(Sxx,Sax,Al,Ar,As):
    combined = pd.DataFrame(columns=['Combined Stress'])
    # 𝜎_𝑎𝑣𝑔,𝑋𝑋,𝑙𝑒𝑓𝑡𝑝𝑖𝑡𝑐ℎ*𝐴_𝑙𝑒𝑓𝑡𝑝𝑖𝑡𝑐ℎ + 𝜎_𝑎𝑣𝑔,𝑋𝑋,𝑟𝑖𝑔ℎ𝑡𝑝𝑖𝑡𝑐ℎ*𝐴_𝑟𝑖𝑔ℎ𝑡𝑝𝑖𝑡𝑐ℎ + 𝜎_𝑎𝑣𝑔,𝑎𝑥𝑖𝑎𝑙,𝑠𝑡𝑟𝑖𝑛𝑔𝑒𝑟*𝐴_𝑠𝑡𝑟𝑖𝑛𝑔𝑒𝑟
    # ---------------------------------------------------------------------------------------------
    # 𝐴_𝑙𝑒𝑓𝑡𝑝𝑖𝑡𝑐ℎ+𝐴_𝑟𝑖𝑔ℎ𝑡𝑝𝑖𝑡𝑐ℎ+𝐴_𝑠𝑡𝑟𝑖𝑛𝑔𝑒𝑟
    
    for i in range(9):
        combined.loc[i] = [(Sxx.loc[i]*Al + Sxx.loc[i+1]*Ar  + Sax.loc[i]*As)/(Al+Ar+As)]
    return combined

## main function
if __name__ == "__main__":
    LC = 'LC2' 
    LC_VON = pd.read_csv(LC+'_VON.txt')
    LC_XX = pd.read_csv(LC+'_XX.txt')
    LC_YY = pd.read_csv(LC+'_YY.txt')
    LC_XY = pd.read_csv(LC+'_XY.txt')
    LC_1D = pd.read_csv(LC+'_1D.txt')
    RF_sf_skin = round(Sult/LC_VON["Contour(Element Stresses (2D & 3D))"],2)
    # changing the index of RF_sf_skin to start with 1
    RF_sf_skin.index = RF_sf_skin.index + 1
    ## Reserve factors against strength failure for the stringers
    RF_sf_stringer = round(Sult/LC_1D["Contour(Element Stresses (1D))"],2)
    # changing the index of RF_sf_stringer to start with 1 
    RF_sf_stringer.index = RF_sf_stringer.index + 1
    
    print("---------------")
    print("RF against Strength Failure for the Skin")
    print("---------------") 
    ## Reserve factors against strength failure for the skin
    print(RF_sf_skin.head(30))
    
    print("---------------")
    print("RF against Strength Failure for the Stringers")
    print("---------------") 
    print(RF_sf_stringer.tail(27))
    
    STRESS_avg_XX = avg_stress(LC_XX,A_x).head(10)
    STRESS_avg_YY = avg_stress(LC_YY,A_y).head(10)
    STRESS_avg_XY = avg_stress(LC_XY,A_x).head(10)
    
    print("---------------")
    print("Average stresses in XX")
    print("---------------")  
    print(STRESS_avg_XX)
    
    print("---------------")
    print("Average stresses in YY")
    print("---------------")  
    print(STRESS_avg_YY)
    
    print("---------------")
    print("Average stresses in XY")
    print("---------------")  
    print(STRESS_avg_XY)
    
    critical_Biaxial = sigma_x_crit(STRESS_avg_YY,STRESS_avg_XX)
    RF_Buckling_Biaxial = get_RF(critical_Biaxial,STRESS_avg_XX)
    
    print("---------------")
    print("Critical Biaxial Stress")
    print("---------------")
    print(critical_Biaxial)
    
    print("---------------")
    print("RF against Panel Buckling Biaxial")
    print("---------------")   
    # changing the index of RF_Buckling_Biaxial to start with 1
    RF_Buckling_Biaxial.index = RF_Buckling_Biaxial.index + 1
    print((RF_Buckling_Biaxial))
    
    critical_Shear = sigma_xy_crit()
    RF_Buckling_Shear = get_RF(critical_Shear,STRESS_avg_XY)
    
    print("---------------")
    print("RF against Panel Buckling in Shear")
    print("---------------")
    # changing the index of RF_Buckling_Shear to start with 1
    RF_Buckling_Shear.index = RF_Buckling_Shear.index + 1
    print(RF_Buckling_Shear)
    
    
    print("---------------")
    print("RF against Panel Buckling Combined")
    print("---------------")   
    RF_Buckling = get_combined_RF(RF_Buckling_Biaxial,RF_Buckling_Shear)
    print(RF_Buckling)    
    
    STRESS_avg_axial_stringer = avg_stress(LC_1D,A_stringer).tail(9)
    STRESS_avg_axial_stringer = STRESS_avg_axial_stringer.reset_index(drop=True)
    STRESS_avg_combined_stringer= combined_stress(STRESS_avg_XX,STRESS_avg_axial_stringer,A_leftpitch,A_rightpitch,A_stringer)
    
    print("---------------")
    print("Average stresses in Axial Stringers")
    print("---------------")  
    print(STRESS_avg_axial_stringer)
    
    print("---------------")
    print("Average stresses Combined in Stringers")
    print("---------------")  
    print(STRESS_avg_combined_stringer)
    
    print("---------------")
    print("RF against Column Panel Buckling ")
    print("---------------")
    print(column_buckling_RF(STRESS_avg_combined_stringer))
    
    print("---------------")
    print("RF against Crippling ")
    print("---------------")
    print(column_buckling_cr_crip_web()/STRESS_avg_combined_stringer)
