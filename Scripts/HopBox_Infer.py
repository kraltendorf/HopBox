import os 
os.environ['KMP_DUPLICATE_LIB_OK']='True'
import sys
import numpy as np

from hopbox_functions import run_inference


Base_folder = "C:/Users/collin/Desktop/HopBox/"
inference_dir = Base_folder + "HopBox 2022/119CANON/"
# list images in directory
# print(images_dir)
filess = os.listdir(inference_dir)
print(str(len(filess)) + " image files found")
print(filess)

Save_Images = 1 # 1 for saving images, 0 for NOT saving images

model_dir = Base_folder + "Models/hoptrial1.pth"

# append output folder with date today and time now
# tn = datetime.datetime.now()
# tnn = str(tn.month) +"_"+str(tn.day) +"_"+str(tn.year)+"_T_"+str(tn.time())
output_directory = Base_folder + "Results_MAY_23/119CANON"

if Save_Images == 1:
  mask_dir = output_directory + "/Masks/"
  label_dir = output_directory + "/Mask_labels/"
  corrected_dir = output_directory + "/Corrected_imgs/"

  try:
    os.makedirs(mask_dir)
    os.makedirs(label_dir)
    os.makedirs(corrected_dir)
  except:
    print("'" + output_directory+ "' already exists")

directories = [inference_dir,
               mask_dir,
               label_dir,
               corrected_dir,
               model_dir,
               output_directory
]
# Load images, preprocess image, color correct, read barcode, predict masks, get region properties

scale = 1; # should be based image size and network input size

# Reference Data 
def load_csv_chart(f):
    '''Input CSV's shape is (24, 3), with sRGB reference values
    '''
    data = np.loadtxt(f,delimiter=",")
    assert data.shape ==(24, 3)
    return data
ref = load_csv_chart(Base_folder + 'ref_babel_avg_srgb_d50.csv') # loads reference RGB color values

# Run Inference
run_inference(
  ref=ref, 
  directories=directories, 
  filess=filess, 
  scale=scale, 
  Save_Images=Save_Images
) 