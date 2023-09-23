from matplotlib import pyplot as plt
from math import pi

E_lis = ["-0.05", "-0.10", "-0.15", "-0.20", "-0.25"]
pi_lis = ["0.100", "0.050", "0.010", "0.005", "0.001"]


for j in range(5):

    var = 1
    
    if(not var):
        name_lis =  [   f"OneBody4_E{ E_lis[i]}_plotinterval{pi_lis[j]}_ecc0.0_RK4.txt" for i in range(5) ]
    else:
        name_lis = [   f"OneBody4_E{E_lis[j]}_plotinterval{pi_lis[i]}_ecc0.0_RK4.txt" for i in range(5) ]

    each_lis_x = list()
    each_lis_y = list()

    for file_name in name_lis:
        
        nume_lis = list()
        anly_lis = list()
        lis3 = list()
        with open( file_name ) as file_contents:
            for line in file_contents:
                try :
                    int(line.split("|")[0])
                except:
                    continue
                
                goo = line.split("|")[1].split(",")
                anly_lis.append((float(goo[1])**2+float(goo[2])**2)**0.5)
                
                
                foo = line.split("|")[2].split(",")
                nume_lis.append((float(foo[1])**2+float(foo[2])**2)**0.5)
                
        nume = sum(nume_lis)/len(nume_lis)
        anly = sum(anly_lis)/len(anly_lis)
        if not var :
            each_lis_x.append( float(file_name.split("_E")[1][:5]) )
            each_lis_y.append(abs(nume-anly)/anly*100)
        else:
            each_lis_x.append( float(file_name.split("_plotinterval")[1][:5]) )
            each_lis_y.append(abs(nume-anly)/anly*100)

    
    if not var:
        plt.title(f"plot interval : {pi_lis[j]}")
        print(f"plot interval : {pi_lis[j]}")
    else:
        plt.title(f"Energy : {E_lis[j]}")
        print(f"Energy : {E_lis[j]}")

    for i in range(len(each_lis_x)):
        print(each_lis_x[i], each_lis_y[i])
    print()

    plt.plot(each_lis_x, each_lis_y, "r-o")
    plt.show()    
    