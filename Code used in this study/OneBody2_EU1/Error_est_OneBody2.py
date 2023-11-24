import numpy as np
import matplotlib.pyplot as plt
import os

def Norm2(x, y):
    return (x**2 + y**2)**0.5
def get_rad(theta, y):
    if y>= 0:
        print("theta : ",theta)
        return np.arccos(theta)
    else:
        print("theta : ",theta)
        return np.pi*2 - np.arccos(theta)


FileList = list()
for i in os.listdir(os.getcwd()):
    if (i.endswith(".txt") and i.startswith("OneBody2")):
        FileList.append(i)


for file_name in FileList:
    new_file_name = f"new_{file_name}"
    print(new_file_name)
    f = open(new_file_name, "w")
    
    init = list()
    ecc =  float(file_name.split("ecc")[1].split("_")[0])
    a, b = int(), int()
    G, M = 1, 1

    
    ang_err_lis = list()
    
    with open( file_name ) as file_contents:      
        for line in file_contents:
            if ("initial value :" in line):
                var1 = line.split("= [")[1].split("]")[0].split(",")
                for i in range(4):
                    init.append(float(var1[i]))      
                
                if not(len(init)):
                    print("Something Error happened!!!")
                    exit()
                
                rp = init[0]
                ra = rp*(1+ecc)/(1-ecc)
                a = (ra+rp)/2
                b = a*(1-ecc**2)**0.5
            

            try :
                int(line.split("|")[0])
            except:
                f.write(line)
                continue
            
            var2 = line.split("|")[1].split(",")
            
            x2 = float(var2[1])
            y2 = float(var2[2])
            
            # angle = get_rad( (a*ecc + Norm2(x2,y2)*np.cos()))/a ,y2)          
            
            k_const = G*M
            l_const = Norm2(init[0], init[1])*Norm2(init[2], init[3])*np.sin(np.arccos((Norm2(init[0], init[1])**2+Norm2(init[2], init[3])**2-Norm2(init[0]-init[2], init[1]-init[3])**2) /(2*Norm2(init[0], init[1])*Norm2(init[2], init[3]))))
            A_coeff =  (1/Norm2(init[0], init[1]) - (k_const/(l_const**2)))/(init[0]/Norm2(init[0], init[1]))
            r_0 = 1/(A_coeff*np.cos(float(var2[0]))+(k_const/(l_const**2)))
            # r_0 = Norm2(a*np.cos(angle) - (a-rp), b*np.sin(angle))
            absolute_error = abs( r_0 - Norm2(x2, y2)  )
            relative_err = 100 * absolute_error / r_0
            
            f.write(f'{int(line.split("|")[0])},{float(var2[0])},{float(var2[1])},{float(var2[2])},{absolute_error},{relative_err}\n')         
            
    f.close()
    print(f"end_{file_name}")


