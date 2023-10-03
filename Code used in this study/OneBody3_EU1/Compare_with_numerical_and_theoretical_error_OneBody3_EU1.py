import numpy as np
import matplotlib.pyplot as plt
import os

def Norm2(x, y):
    return (x**2 + y**2)**0.5

#theoretical error without direction
def theoretical_error(kappa , d):   
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
    a, b = float(), float()

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
                    
                    var2 = line.split("|")[1].split(",")
                    geod_err.append( [ float(var2[1]), float(var2[2]), var2[3]] )
                else:
                    var2 = line.split("|")[1].split(",")
                    newt_err.append( [ float(var2[1]), float(var2[2]), var2[3]] )                     


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
        
    
    #accumulate the theoretical error
    geod_theorerical_error_vecsum = [0,0]
    newt_theorerical_error_vecsum = [0,0]
    for i in range(1, len(geod_err)):
        lis = error_vec(geod_err[i][0], geod_err[i][1], Norm2(geod_err[i][0] - geod_err[i-1][0] , geod_err[i][1] - geod_err[i-1][1] ) )
        for j in range(2):
            geod_theorerical_error_vecsum[j] += lis[j] 
    for i in range(1, len(newt_err)):
        lis = error_vec(newt_err[i][0], newt_err[i][1], Norm2(newt_err[i][0] - newt_err[i-1][0] , newt_err[i][1] - newt_err[i-1][1] ) )
        for j in range(2):
            newt_theorerical_error_vecsum[j] += lis[j]
    

    #wirte a result
    var3 = float(geod_err[-1][2])
    var4 = Norm2(geod_theorerical_error_vecsum[0], geod_theorerical_error_vecsum[1])
    f.write(f"\n\n-geodesic accmulated error-\nnumerical   : {var3}\ntheoretical : {var4}\ndifference : {100*(var3-var4)/var3}%\n\n")
    var3 = float(newt_err[-1][2])
    var4 = Norm2(newt_theorerical_error_vecsum[0], newt_theorerical_error_vecsum[1])
    f.write(f"-newtonian accmulated error-\nnumerical   : {var3}\ntheoretical : {var4}\ndiffernece : {100*(var3-var4)/var3}%\n")
            
            
    f.close()
    print(f"end_{file_name}")
