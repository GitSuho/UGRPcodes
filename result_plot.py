
import matplotlib.pyplot as plt

x_val = []
y_val = []
temp_x = []
temp_y = []

with open("3bodysimulresult.txt") as file_data:
    for line in file_data:
        if line.startswith("Position of object"):
            in_line = line.split(":")
            if in_line[0][-1] == "1" or in_line[0][-1] == "2":
                num_line = in_line[1].split()
                temp_x.append(float(num_line[0].strip(",")))
                temp_y.append(float(num_line[1].strip(",")))
            elif in_line[0][-1] == "3":
                num_line = in_line[1].split()
                temp_x.append(float(num_line[0].strip(",")))
                temp_y.append(float(num_line[1].strip(",")))
                x_val.append(temp_x.copy())  # Use a copy of temp_x
                y_val.append(temp_y.copy())  # Use a copy of temp_y
                temp_x = []
                temp_y = []
            else:
                print("ERROR:12087dghf27hd!!!")

x_val = list(map(list, zip(*x_val)))  # Transpose x_val
y_val = list(map(list, zip(*y_val)))  # Transpose y_val

plt.plot(x_val[0], y_val[0], color='red')
plt.plot(x_val[1], y_val[1], color='green')
plt.plot(x_val[2], y_val[2], color='blue')

plt.show()
