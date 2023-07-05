import os
import csv
import pandas as pd
import numpy as np

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
        sigma_xy_crit.loc[i] = [(E*PI**2)/(12*(1-nu**2)) * (t/b)**2 * optimize_K_Shear()]
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
    stress = df[df.columns[-2]]*A
    print('AHAHAHAHAHA')
    print(df[df.columns[4]])
    avg = (stress.rolling(window=3, min_periods=1).sum())/(3*A)
    # Exclude non-multiple of 3 rows from the resulting object
    avg = avg[df.index % 3 == 2]
    # Reset index starting from 1
    avg = avg.reset_index(drop=True)
    return avg

def avg_stress_1d(df,A):
    # Calculate the sum of every three consecutive rows
    stress = df['Contour(Element Stresses (1D))']*A
    avg = (stress.rolling(window=3, min_periods=1).sum())/(3*A)
    # Exclude non-multiple of 3 rows from the resulting object
    avg = avg[df.index % 3 == 2]
    # Reset index starting from 1
    avg = avg.reset_index(drop=True)
    return avg

## main function
if __name__ == "__main__":
    LC = 'LC2' 
    LC_VON = pd.read_csv(LC+'_VON.txt')
    LC_XX = pd.read_csv(LC+'_XX.txt')
    LC_YY = pd.read_csv(LC+'_YY.txt')
    LC_XY = pd.read_csv(LC+'_XY.txt')
    LC_1D = pd.read_csv(LC+'_1D.txt')
    print("---------------")
    print("RF against Strength Failure for the Skin")
    print("---------------") 
    ## Reserve factors against strength failure for the skin
    RF_sf_skin = round(Sult/LC_VON["Contour(Element Stresses (2D & 3D))"],2)
    # changing the index of RF_sf_skin to start with 1
    RF_sf_skin.index = RF_sf_skin.index + 1
    print(RF_sf_skin.head(30))
    
    print("---------------")
    print("RF against Strength Failure for the Stringers")
    print("---------------") 
    ## Reserve factors against strength failure for the stringers
    RF_sf_stringer = round(Sult/LC_1D["Contour(Element Stresses (1D))"],2)
    # changing the index of RF_sf_stringer to start with 1 
    RF_sf_stringer.index = RF_sf_stringer.index + 1
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
    print("---------------")
    print("HRHRHRHRHRHRHH")
    print("---------------")
    
    print(critical_Biaxial)
    RF_Buckling_Biaxial = get_RF(critical_Biaxial,STRESS_avg_XX)
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
    