import sys
import pandas as pd
import matplotlib.pyplot as plt
file1=sys.argv[1]
file2=sys.argv[2]
outfile=sys.argv[3]
#title = sys.argv[4]
plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"

data1 = pd.read_csv(file1)
data2 = pd.read_csv(file2)
if 'wlog' in file1:
    MIN_VAL = 0
    MAX_VAL = 1000
elif 'Tc' in file1:
    MIN_VAL = -1
    MAX_VAL = 30
elif 'lambda' in file1:
    MIN_VAL = -0.05
    MAX_VAL = 3.0
else:
    print("file not found")
data = [data1,data2]
color = ['r', 'b']
label = ['CGCNN', 'ALIGNN']
symbol = ['o', '^']
fig,ax = plt.subplots()
for i,_ in enumerate(data):
    d = data[i]
    x = d['target']
    y = d['prediction']
    plt.scatter(x, y,color=color[i],label=label[i],marker=symbol[i],s=50)
ax.plot([MIN_VAL,MAX_VAL], [MIN_VAL,MAX_VAL], "k--", lw=2)
plt.legend(loc='upper right',fontsize=20)
#plt.title(title,fontsize=20)
plt.xlabel("DFPT", fontsize=20)
plt.ylabel("ML prediction", fontsize=20)
ax.set_xlim(MIN_VAL,MAX_VAL)
ax.set_ylim(MIN_VAL,MAX_VAL)
ax.set_aspect('equal')
plt.savefig(outfile)
