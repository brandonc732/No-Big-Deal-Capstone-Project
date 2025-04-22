#%%

import numpy as np
import matplotlib.pyplot as plt
from skimage.io import imread
from skimage.exposure import adjust_gamma
from skimage.exposure import is_low_contrast

#%%

def read_plot_frame(num, gamma=0.5, clip_quantile=0.9999):
    tmp = imread("1. Fluo_0.1s_mid_t"+num+"_ORG.tif")
    tmp[tmp > np.quantile(tmp, clip_quantile)] = np.quantile(tmp, clip_quantile)
    tmp = adjust_gamma(tmp, gamma)
    print(is_low_contrast(tmp))
    plt.imshow(tmp)
    plt.show()

# %%

read_plot_frame("050", gamma=2)

# %%
