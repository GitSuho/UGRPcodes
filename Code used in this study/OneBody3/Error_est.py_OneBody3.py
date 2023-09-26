import numpy as np

def Norm2(x, y):
    return (x**2 + y**2)**0.5


for file_name in ["OneBody3_E-0.25_ecc0.8_RK4.txt" ]:
    new_file_name = f"{file_name[-4]}_new.txt"
    f = open(new_file_name, "w")
    
    init = list()
    ecc =  float(file_name.split("ecc")[1].split("_")[0])
    a, b = int(), int()

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
                b = a/(1-ecc**2)**0.5
                                       

            try :
                int(line.split("|")[0])
            except:
                f.write(line)
                continue
            
            var2 = line.split("|")[1].split(",")
            angle = float(var2[0])
            r_0 = Norm2(a*np.cos(angle), b*np.sin(angle))
            absolute_error = abs( r_0 - Norm2(float(var2[1]), float(var2[2]))  )
            relative_err = 100 * absolute_error / r_0
                       
            
            
            
            f.write(f'{line.split("|")}|{var2[0]},{var2[1]},{var2[2]},{absolute_error},{relative_err}')         
            




    f.close()
    print("end")
