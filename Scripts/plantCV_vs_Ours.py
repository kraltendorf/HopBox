# compare our color color correction method with plantCV's color correction method

import cv2
import numpy as np
import time
import  glob
import matplotlib.pyplot as plt

from hopbox_functions import *
from plantcv import plantcv as pcv


def load_csv_chart(f):
    '''Input CSV's shape is (24, 3), with sRGB reference values
    '''
    data = np.loadtxt(f,delimiter=",")
    assert data.shape ==(24, 3)
    return data

def color_correct_pcv(img,ref):
    '''Color corrects image using plantCV's color correction method
    '''
    _, start, space = pcv.transform.find_color_card(rgb_img=ref, background='light')
    target_mask = pcv.transform.create_color_card_mask(ref, radius=15, start_coord=start, 
                                                       spacing=space, nrows=6, ncols=4)
    source_mask = pcv.transform.create_color_card_mask(img, radius=10, start_coord=start, 
                                                       spacing=space, nrows=6, ncols=4)
    
    tm, sm, _, corrected_img = pcv.transform.correct_color(target_img=ref, 
                                                           target_mask=target_mask,
                                                           source_img=img,
                                                           source_mask=source_mask,
                                                           output_directory='pcv_results')
    
    # pcv.transform.quick_color_check(source_matrix=sm, target_matrix=tm, num_chips=24)
    
    return corrected_img
    

# load image folders with images to be corrected
img_folder = "C:/Users/collin/Desktop/HopBox/HopBox 2022/119CANON/"
img_list = glob.glob(img_folder + "*.JPG")

ref = load_csv_chart('ref_babel_avg_srgb_d50.csv') # loads reference RGB color values

img_ref = cv2.imread('Img_.png') # loads reference image

# loop through images in folder do color correction using plantCV and our method and show the images
for img_path in img_list:
    img = cv2.imread(img_path)
    
    # color correct using plantCV
    tic = time.time()
    img_pcv = color_correct_pcv(img,img_ref)
    pcv_t = (time.time()-tic)*1000
    
    # color correct using our method
    tic = time.time()
    img_ours,_ = correct_color(img,ref)
    our_t= (time.time()-tic)*1000
    
    img_number = img_list.index(img_path)+1
    
    # assign values to dataframe
    df = pd.DataFrame(columns=['img_name', 'image_number','pcv_time','our_time'])
    df.loc[1] = [img_path,img_number,pcv_t,our_t]
    
    
    # print df in table format
    print(df)
    
# show last image side by side using matplotlib
# name figure as figure 2

plt.figure(2, figsize=(10,10))
plt.subplot(1,2,1)
plt.imshow(cv2.cvtColor(img_pcv,cv2.COLOR_BGR2RGB))
plt.title("plantCV")
plt.subplot(1,2,2)
plt.imshow(img_ours)
plt.title("Ours")
plt.show()