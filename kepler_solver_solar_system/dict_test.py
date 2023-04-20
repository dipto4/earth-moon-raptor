import numpy as np

data = np.loadtxt("radii")

dict_r = {}
obj_names = ["Sun","Venus","Earth","Mars","Jupiter","Saturn","Uranus","Neptune"]

for i in range(0,len(obj_names)):
    dict_r[obj_names[i]] = data[i]

print(dict_r)
