import numpy as np
import matplotlib.pyplot as plt


E = -0.05
ecc = [0, 0.4, 0.8]
y_val_geod = []
y_val_newt = []

FileList = ["new_OneBody2_E-0.05_ecc0.0_EU1.txt", "new_OneBody2_E-0.05_ecc0.4_EU1.txt", "new_OneBody2_E-0.05_ecc0.8_EU1.txt"]

for file_name in FileList:

    geod_err = list()
    newt_err = list()
    
    foo = 1
    goo = 2
    
    #read file line by line. And append data on the lists
    with open( file_name ) as file_contents:      
        for line in file_contents:
            if foo:
                if("num|degree" in line):
                    foo = 0
            else:
                if goo:
                    if (int(line.split("|")[0]) == 1):
                        goo -= 1
                        continue
                    var2 = line.split("|")[1].split(",")
                    geod_err.append(float(var2[4]))
                else:
                    var2 = line.split("|")[1].split(",")
                    newt_err.append(float(var2[4]))                    
           
    #eliminate values that over 2pi
    geod_err.pop()
    newt_err.pop()

    
    print(file_name)
    print(f'plot point num : geod = {len(geod_err)} , newt = {len(newt_err)}')
    
    y_val_geod.append(np.nanmean(geod_err))
    y_val_newt.append(np.nanmean(newt_err))
    


plt.ylabel("mean of realtive error")
plt.xlabel("eccentricity")

X1 = [1, 2.5, 4]
X2 = [1.5, 3, 4.5]

plt.bar(X1, y_val_geod, color = 'r', width=0.5 , label="Geodesic")
plt.bar(X2, y_val_newt, color = 'b', width=0.5 , label="Newtonian")
plt.legend(loc = 'upper left')

X_mid  = [(X1[i] + X2[i])/2 for i in range(3)]
ticklabel = ['0.0', '0.4', '0.8']
plt.xticks(X_mid,ticklabel,fontsize=15,rotation=0 )


plt.tick_params(
    axis='x',              
    bottom=False)
plt.show()