import numpy as np
import matplotlib.pyplot as plt
import os

def Norm2(x, y):
    return (x**2 + y**2)**0.5
def get_rad(theta, y):
    if y>= 0:
        return np.arccos(theta)
    else:
        return np.pi*2 - np.arccos(theta)

#theoretical error without direction
def theoretical_error(kappa , d):
    # r = 1/kappa
    
    # m2 = -2*r/d
    
    # dx3 = -d/(2*(1+m2**2)**0.5)
    # m3 = -(r+dx3)/(m2*dx3)
    
    # dx4 = -d/(1+m3**2)**0.5
    # m4 = -(r+dx4)/(m3*dx4)
    
    # tx2 = -d/(1+m2**2)**0.5
    # tx3 = -d/(1+m3**2)**0.5
    # tx4 = -d/(1+m4**2)**0.5
    
    # r_nume = ((r + tx2/3 + tx3/3 + tx4/6)**2 + (d/6 + m2*tx2/3 + m3*tx3/3 + m4*tx4/6)**2)**0.5
    
    r = 1/kappa
    return abs(Norm2(r, d) - r)




###main###
FileList = list()
for i in os.listdir(os.getcwd()):
    if (i.endswith(".txt") and i.startswith("new_")):
        FileList.append(i)
for file_name in FileList:
    new_file_name = f"result_{file_name}"
    f = open(new_file_name, "w")
    f.write(new_file_name)
    
    init = list()
    ecc =  float(file_name.split("ecc")[1].split("_")[0])
    a, b = int(), int()

    geod_err = list()
    newt_err = list()
    
    foo = 1
    goo = 2
    
    #read file line by line. And append data on the lists
    with open( file_name ) as file_contents:      
        for line in file_contents:
            if foo:
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
                elif("num|degree" in line):
                    foo = 0
            else:
                if goo:
                    if (int(line.split("|")[0]) == 1):
                        goo -= 1
                        continue
                        if not goo :
                            var2 = line.split("|")[1].split(",")
                            newt_err.append( [ float(var2[1]), float(var2[2]), var2[3]] )   
                            continue
                    var2 = line.split("|")[1].split(",")
                    
                    x2 = float(var2[1])
                    y2 = float(var2[2])
                    angle = get_rad( (a*ecc + Norm2(x2,y2)*np.cos(float(var2[0])))/a ,y2)         
                    x1, y1 = a*np.cos(angle) - (a-rp) , b*np.sin(angle)
                    
                    geod_err.append( [ x2, y2, var2[3], x1, y1] )
                else:
                    var2 = line.split("|")[1].split(",")
                    
                    x2 = float(var2[1])
                    y2 = float(var2[2])
                    angle = get_rad( (a*ecc + Norm2(x2,y2)*np.cos(float(var2[0])))/a ,y2)         
                    x1, y1 = a*np.cos(angle) - (a-rp) , b*np.sin(angle)
                    
                    newt_err.append( [ x2, y2, var2[3], x1, y1] )                    


    #return curvature based on the current position
    def Find_curvature(x, y):
        return a*b/(x**2 + y**2)**1.5
    
    #theoretical error of vector element
    def error_vec(x0, y0, d ):
        x = x0 + a - rp
        y = y0
        sca_err = theoretical_error(Find_curvature(x, y) , d)
        
        m_perpen =  (a**2*y)/(b**2*x)
        dx = d/(1 + m_perpen**2)**0.5
        
        if x >= 0 :
            result = [ sca_err*dx, sca_err*dx*m_perpen ]
        else:
            result = [-sca_err*dx , -sca_err*dx*m_perpen]
            
        return result
       
    
    #eliminate values that over 2pi
    geod_err.pop()
    newt_err.pop()

    #eliminate nan values
    while(1):
        if(geod_err[-1][2] == "nan"):
            geod_err.pop()
        else :
            break
    while(1):
        if( newt_err[-1][2] == "nan"):
            newt_err.pop()
        else :
            break

    
    g_rela_err = list()
    n_rela_err = list()
    
    f.write("\ngeod|\nnewt|\nnum|numerrical|theoretical|relative difference(%)\n")
    for i in range(1, len(geod_err)):
        hoo = theoretical_error(Find_curvature(geod_err[i-1][0] + a -rp, geod_err[i-1][1]) , Norm2(geod_err[i-1][3] - geod_err[i][0], geod_err[i-1][4] - geod_err[i][1]))
        ioo = float(geod_err[i][2])
        f.write(f"{i}|{ioo},{hoo},{abs(100*(ioo-hoo)/ioo)}\n")
        g_rela_err.append(abs(100*(ioo-hoo)/ioo))
        

    for i in range(1, len(newt_err)):
        hoo = theoretical_error(Find_curvature(newt_err[i-1][0] + a -rp, newt_err[i-1][1]) , Norm2(newt_err[i-1][3] - newt_err[i][0], newt_err[i-1][4] - newt_err[i][1]))
        ioo = float(newt_err[i][2])
        f.write(f"{i}|{ioo},{hoo},{abs(100*(ioo-hoo)/ioo)}\n")
        n_rela_err.append(abs(100*(ioo-hoo)/ioo))
    
    f.write(f"\n\ngeodesic relative error mean : {np.nanmean(g_rela_err)}\n")
    f.write(f"newtonian relative error mean : {np.nanmean(n_rela_err)}\n")
                
            
    f.close()
    print(f"end_{file_name}")
