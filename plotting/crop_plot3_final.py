#!/usr/bin/env python
# coding: utf-8

from PIL import Image
import os.path, sys

Folder_start = 0
number_folder = 5

# PATHS DEFINITION
path = './plot-waveform/'
target = '../data/ENZ_crop/' ## define target
dirs = os.listdir(path)
dirs = sorted(dirs)[Folder_start:Folder_start+number_folder]


# In[5]:


from PIL import Image
import os.path, sys

def crop(path,target,dirs):
    for item in dirs:
        fullpath = os.path.join(path,item) #corrected
        if os.path.isfile(fullpath):
            im = Image.open(fullpath)
            f, e = os.path.splitext(item)
            imCrop = im.crop((96, 96, 519, 712)) #corrected
            imCrop.save(target + f + '_crop.png', "png", quality=100)


# MAKE OUTPUT FOLDER
try:
    os.mkdir(target)
    print('destination created')
except:
    print('destination alredy exist')

iterasi = 1

for dir_ in dirs:
    try:
        os.mkdir(target+dir_+'crop')
        print(dir_+'crop','created', end=" ")
    except:
        print(dir_+'crop','already exist', end=" ")  
        if len(os.listdir(target+dir_+'crop/'))>5:
            continue
     
    print("{} start cropping...({}/{})".format(dir_,iterasi,len(dirs)))
    iterasi += 1
    pathdir = path+dir_+'/'
    targetdir = target+dir_+'crop/'
    imagedirs = os.listdir(pathdir)
    try:
        crop(pathdir,targetdir,imagedirs)
    except:
        print(dir_+'erorr',' exist', end=" ") 
