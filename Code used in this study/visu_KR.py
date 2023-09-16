from matplotlib import pyplot as plt
from math import pi

for file_name in []:
    
    
    geod_theta_err = [[],[]]
    newt_theta_err = [[],[]]

    period = -1

    with open( file_name ) as file_contents:
        for line in file_contents:
            if line.split()[0] == "cycle":
                # if line.split()[1] == "2":
                #     break 
                # else:
                    period += 1
                    continue
            
            equ_typ = line.split("|")
            geod_line = equ_typ[0].split(",")
            newt_line = equ_typ[1].split(",")
            
            if (geod_line[0] == geod_line[1] and geod_line[2] == geod_line[3]):
                pass
            else:
                geod_theta_err[0].append(float(geod_line[2]) + period*pi*2)
                geod_theta_err[1].append(float(geod_line[3]))
            
            if (newt_line[0] == newt_line[1] and newt_line[2] == newt_line[3]):
                pass
            else:
                newt_theta_err[0].append(float(newt_line[2]) + period*pi*2)
                newt_theta_err[1].append(float(newt_line[3]))


    print(f' g_err max :{max(geod_theta_err[1])}, n_err max :{max(newt_theta_err[1])}')

    plt.plot(geod_theta_err[0][:], geod_theta_err[1][:], 'or')
    plt.plot(newt_theta_err[0][:], newt_theta_err[1][:], "+b")

    plt.savefig(file_name[:-4]+".png", dpi = 3000)

    plt.show()
    plt.pause(1)