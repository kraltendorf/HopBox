
# Load more dependencies
from __future__ import print_function

import io as io
import os
import sys
import time

import cv2
import numpy as np
import pandas as pd
import pyzbar.pyzbar as pzb
import scipy as scipy
import torch
import torch.nn as nn
import torchvision
import torchvision.transforms as T
from google.colab.patches import cv2_imshow
from PIL import Image
from scipy import ndimage as ndi
from scipy.ndimage import morphology
from skimage import io, morphology
from skimage.color import label2rgb
from skimage.measure import label, regionprops
from sklearn.cross_decomposition import PLSRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import PolynomialFeatures
from torchvision.models.segmentation.fcn import FCNHead


# Materials class
class Material:
  def __init__(self, name, input_rgb_vals, output_val, confidence_threshold=0):
    self.name = name
    self.input_rgb_vals = input_rgb_vals
    self.output_val = output_val
    self.confidence_threshold = confidence_threshold


# resizes image based on scale
def resize_image(image,scale): 
  w, h = image.size
  newW, newH = int(scale * w), int(scale * h)
  image1 = image.resize((newW, newH))

  return image1, newW, newH

# Read QR code or bar code from images
def read_QR_bar(image): 

  detectObjects = pzb.decode(image)
  obj_type = ""
  obj_data = ""

  if len(detectObjects)>0:
    for obj in detectObjects:
      obj_type = obj.type
      obj_data = obj.data
  else:
    print("Could NOT find Bar/QR code in image!!!")

  return detectObjects, obj_type, obj_data

# displays image with detected barcode or QR code
def displayQR(im, decodeObjects): 
  if len(decodeObjects)>0:
      for decodeObject in decodeObjects:
          points = decodeObject.polygon
          if len(points)>4:
              hull = cv2.convexHull(np.array([point for point in points],dtype = np.float32))
              hull = list(map(tuple,np.squeeze(hull)))
          else:
              hull = points
          n = len(hull)
          hull = np.int0(hull)

          for j in range(0,n):
              cv2.line(im, hull[j], hull[(j+1) % n],(0,255,0),2)
          
      cv2_imshow(im)
      cv2.waitKey(0)

  else:
    print("Could NOT find Bar/QR code in image!!!")
    
# get RGB values at specified coordinates
def get_RGB_values(image,coords):
  RGB_values = []
  for iz in range(len(coords)-1):
    RGB_vals = image[coords[iz,0],coords[iz,1]]
    RGB_values.append(RGB_vals)

  RGB_values = tuple(RGB_values)
  RGB_values = abs(np.asarray(RGB_values))
  RGB_values[RGB_values == 0] = 1
  return RGB_values

#Extract Color chart values from image
def extract_color_chart(img):
    detector = cv2.mcc.CCheckerDetector_create()
    detector.process(img, cv2.mcc.MCC24,1)
    checkers = detector.getListColorChecker()

    for checker in checkers:
            cdraw = cv2.mcc.CCheckerDraw_create(checker)
            img_draw = img.copy()
            cdraw.draw(img_draw)
            
            chartsRGB = checker.getChartsRGB()
            width, height = chartsRGB.shape[:2]
            roi = chartsRGB[0:width,1]
            # print (roi)
            rows = int(roi.shape[:1][0])
            src = chartsRGB[:,1].copy().reshape(int(rows/3), 1, 3)
            src = src.reshape(24,3)
    return src, img_draw

# Correct image color
def correct_color(img, ref): 
  im_rgb = cv2.cvtColor(img, cv2.COLOR_BGR2RGB)
  src, img2 = extract_color_chart(img)
  model = Pipeline([('poly', PolynomialFeatures(degree = 3)),
                  ('pls', PLSRegression(n_components = 12))])

  model.fit(src, ref)
  print("Color model Fit score = ", model.score(src, ref))

  corr_img = morph = im_rgb.copy()
  corr_img = corr_img.astype(np.uint32)

  for im in corr_img:
      im[:]=model.predict(im[:])

  corr_img[corr_img >= 255] = im_rgb[corr_img >= 255]
  corr_img[corr_img <= 0] = im_rgb[corr_img <= 0]
  corr_img = corr_img.astype(np.uint8)

  return corr_img

# gets region property and color params, outputs them in a data frame
def get_morph_xtics(mat_mask, image1, cor_img, fn, qr):

    selem = morphology.disk(3)
    mat_mask = morphology.binary_opening(mat_mask, footprint=selem) # use morphological opening
    # mat_mask = morphology.binary_closing(mat_mask, selem=selem) # use morphological closing
    shX,shY,Ch = image1.shape
    min_fac = 0.0011 ######################
    mat_mask = morphology.remove_small_objects(mat_mask, int(min_fac*shX*shY)) 
    mat_mask = ndi.binary_fill_holes(mat_mask)

    img2 = label(mat_mask, background=0)
    Labelled_img = label2rgb(img2, image1, alpha=0.5, bg_label=0, bg_color = (0,0,0))
    Labelled_img = 255*Labelled_img
    Labelled_img.astype(np.uint8)

    np_mat = np.array(mat_mask)
    label_mat = label(np_mat, background = 0)

    if np.sum(label_mat) < 1:
        label_mat = np.ones_like(label_mat)
    elif np.sum(label_mat) >= 1:
        label_mat = label_mat

    props_table = regionprops(label_mat)

    df_input = {'File_name':[],
                'QR_name':[], 
                'Label':[], 
                'Location_X':[],
                'Location_Y':[], 
                'Area':[], 
                'Perimeter':[], 
                'Major_axis':[], 
                'Minor_axis':[], 
                'Orientaion':[], 
                'Eccentricity':[], 
                'Solidity':[], 
                'Median_R':[],
                'Median_G':[],
                'Median_B':[],
                'std_R':[],
                'std_G':[],
                'std_B':[],
                'Mean_R':[],
                'Mean_G':[],
                'Mean_B':[]}

    df = pd.DataFrame(df_input)

    ix = 0
    for obj in props_table:
        centX,centY = obj.centroid
        coordz = obj.coords
        colors = get_RGB_values(cor_img,coordz)

        stdRGB = np.nanstd(colors, axis = 0)
        medRGB = np.nanmedian(colors, axis = 0)
        avgRGB = np.nanmean(colors, axis = 0)

        hop_label = obj.label

        font = cv2.FONT_HERSHEY_SIMPLEX
        pos = (round(centY),round(centX))
        fontScale = 3
        color = (0, 0, 0) 
        thickness = 4
        
        Labelled_img = image = cv2.putText(Labelled_img, str(hop_label), pos, font, 
                          fontScale, color, thickness, cv2.LINE_AA)
        
        major_al = round(obj['major_axis_length'])
        minor_al = round(obj['minor_axis_length'])

        Data_list = [
            fn,
            qr,
            obj.label,
            round(centX),
            round(centY),
            obj.area,
            obj.perimeter,
            major_al,
            minor_al,
            obj.orientation,
            obj.eccentricity,
            obj.solidity,
            medRGB[0],
            medRGB[1],
            medRGB[2],
            stdRGB[0],
            stdRGB[1],
            stdRGB[2],
            avgRGB[0],
            avgRGB[1],
            avgRGB[2]
        ]
        df.loc[len(df)] = Data_list
        ix = ix+1


    return Labelled_img, mat_mask, df


def run_inference(ref, directories, filess, scale, Save_Images=True):
    """
    Run inference on the images in the given directory.
    Save the output images in the inference directory.
    Save the output data in csv file.
    
    args:
        ref: reference color chart
        directories: list of directories[inference, image, mask, model, and output directory]
        filess: list of files in the image directory
        scale: scale factor for the image
        Save_Images: boolean to save the output images
        
    """
    materials = [
             Material("background", [255,255,255], 255, 0.5),
             Material("hops", [0,0,0], 1, 0.7),
             ]
    num_materials =len(materials)

    model=torchvision.models.segmentation.fcn_resnet101(pretrained=False)
    model.classifier=FCNHead(2048, num_materials)

    kd = 0
    if torch.cuda.is_available():
        device = torch.device('cuda')
        kd = 1
    else:
        device = torch.device('cpu')
    
    outputs = []
    model.to(device)

    if kd == 1: # cuda device available
        model.load_state_dict(torch.load(directories[4]), strict=False)

    model.train()
    
    DF = pd.DataFrame()
    for file in filess:
        # 1. load and preprocess image
        image = Image.open(directories[0] + file)
        
        w, h = image.size
        print(image.size)
        image, newW, newH = resize_image(image,scale)

        # 2. Color correction
        img = cv2.cvtColor(np.array(image), cv2.COLOR_RGB2BGR)
        corrected_img = correct_color(img, ref)


        # 3. Read barcode
        decodeObjects, obj_type, plt_name = read_QR_bar(image)
        # displayQR(image,decodeObjects)

        if (len(decodeObjects)<1)|(plt_name == ""):
            plt_name = file

        # 4. Preprocess image for DL inference
        image1 = np.array(image, dtype='ubyte')
        image2 = image1

        image = np.array(image, dtype=float)
        new_im = np.zeros((3, newH, newW))
        new_im[0,:,:] = image[:,:,0]
        new_im[1,:,:] = image[:,:,1]
        new_im[2,:,:] = image[:,:,2]
        image = new_im
        
        # Either
        mean = [0.6654, 0.6710, 0.6492] # exerpt from taining code
        std = [0.2726, 0.2342, 0.2406]

        image = torch.from_numpy(image)
        image = T.Normalize(mean = mean, std = std)(image)
        
        image.unsqueeze_(0)
        image = image.to(device = device, dtype=torch.float32)

        # 5. Predict mask using DL model
        tic = time.time()
        with torch.no_grad():
            mask = model(image)['out']
            mask = nn.Sigmoid()(mask)

        toc = time.time()
        print('Inference time: '+str(toc-tic))

        k = 1
        min_conf = 0.7
        mat_mask = mask.cpu().detach().numpy()[0,k,:,:]
        mat_mask[mat_mask >= min_conf] = 1
        mat_mask[mat_mask < min_conf] = 0

        Labelled_img, mat_mask, df = get_morph_xtics(mat_mask, image1, corrected_img, file, plt_name)
        DF = pd.concat([DF,df])

        # save mask, labelled image, and segmented image
        if Save_Images == 1:
            io.imsave(directories[1] + file + "_mask.png", mat_mask)
            io.imsave(directories[2] + file + "_label_mask.png", Labelled_img.astype(np.uint8))
            io.imsave(directories[3] + file +"_corrected.png", corrected_img.astype(np.uint8))
            
    csv_name = directories[5] + "Results.csv"
    DF.to_csv(csv_name)
    print('Data saved to: '+csv_name)