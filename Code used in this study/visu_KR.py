from matplotlib import pyplot as plt
from math import pi

for file_name in ["OneBody1_E-0.05_ecc0.0_RK4.txt", "OneBody1_E-0.10_ecc0.0_RK4.txt", "OneBody1_E-0.15_ecc0.0_RK4.txt", "OneBody1_E-0.20_ecc0.0_RK4.txt","OneBody1_E-0.25_ecc0.0_RK4.txt" ]:
    
    lis = list()
    lis2 = list()
    lis3 = list()
    with open( file_name ) as file_contents:
        for line in file_contents:
            try :
                int(line.split("|")[0])
            except:
                continue
            
            goo = line.split("|")[1].split(",")
            lis2.append((float(goo[1])**2+float(goo[2])**2)**0.5)
            
            
            foo = line.split("|")[2].split(",")
            lis.append((float(foo[1])**2+float(foo[2])**2)**0.5)
            
            
            lis3.append(float(foo[3]) )


    # print(lis)
    plt.plot(range(len(lis)), lis)
    # plt.plot(range(len(lis3)), lis3, "g.-")
    plt.plot(range(len(lis2)), lis2, "r:")

    plt.ylim(lis[0]-1, lis[0]+1)




    plt.show()
