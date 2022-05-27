import numpy as np
import os, sys, pickle, glob, re
from scipy.io import savemat

loc = 'sweep/'

# Angles
q12,q23 = [],[]

for f in glob.glob(loc+'dos_q12_*'):
    s = re.split('q12_|_q23_|_kcut',f)
    q12.append(float(s[1]))
    q23.append(float(s[2]))

q_dict = {'q12':np.array(q12),'q23':np.array(q23)}

savemat("q_dict.mat", q_dict)
    


