# Beltrami-Coefficient-Map
This pipeline, consisting of three major steps: 
1. cortical surface conformal parameterization, 
2. surface-spline-based cortical activation signal smoothing, 
3. vertex-wise Beltrami coefficient-based map description. 

By correcting most of the topological violations, the result was a “Beltrami coefficient map”(BCM) that completely characterizes the retinotopic map by quantifying the local quasiconformal mapping distortion at each visual field location. The BCM provided topological and fully reconstructable retinotopic maps.

**Human Connectome Project (HCP) data** 

  Download [Data](https://osf.io/5hvg6/files/) (data_step0) to replicate our work.
  
 ![Setup_data1](https://user-images.githubusercontent.com/97193844/174275329-509203d5-a04c-491c-ad7a-17ba6cae3bdd.jpg)
 
  Extract it into the project folder.
  
**Before starting the processing**

Before starting the procedure you need to add both dataset and toolboxes to matlab path. By running the startup.m code from the project directory you can easily be ready to run the project. Alternatively you can run the following commands in command window.

```
addpath(genpath('utilties'))
addpath('gsl_retinotopic')
addpath(genpath('data_step0'))
```

**Note:** __Since each step uses the result of previous step, all steps need to be run consecutively.__ 

## Run each step on HCP data
You can apply each step on one selected subject or all 181 subjects at once.

**Step 1: Conformal Flattening**
  
To apply step 1 which cut a patch from the cortical surface that contains V1 on a single user selected subject (by GUI) use the following command.
  
```
step1_cut_flat(0);
```  

You can run the step in batch mode (all subjects in folder 'data_step0'. 

```
step1_cut_flat(1);
```  

Folders ‘data_step1’ and ‘data_step0_subdiv’ will be made in the working directory which contains output of the step and the updated data with subdivided meshes, respectively.

**Step 2: Thin Plate Spline Smoothing**

This step, in contrast with previous step needs another input which specifies the percentage of the r-square threshold. To run this step which smooth the pRF eccentricity and polar angle data on V1 of the cortical surface patch using thin plate spline on one selected output of previous step, use command:

```
step2_smooth(0, 25);
``` 

Same as previous step you can run this step on all subjects at once. To do so use command below.

```
step2_smooth(1, 25);
``` 

Folder ‘data_step2’ will be made in the current directory by running the command above on the selected data. The folder contains three sub-folders ‘figures’, ‘smooth’ and ‘unsmooth’. You can find the smoothed and unsmoothed data in the sub-folder ‘smooth’ and ‘unsmooth’ respectively.

You can find the following figures in 'figures' folder.

![Figure1](https://user-images.githubusercontent.com/97193844/174284121-700779cd-a40b-4908-9d0a-1ba7c0feb0e8.jpg)

Also you can find the contours for both eccentricity and angle of the V1.

![Figure2](https://user-images.githubusercontent.com/97193844/174310084-57420312-b7fc-44bd-82e7-b1ad9acbd090.jpg)

**Step 3: Beltrami Coefficient Computing**

To compute the Beltrami coefficient of the smoothed V1 retinotopic map on one single selected subject, use following command.

```
step3_compute_bcm(0, 25);
``` 

Same as other two previous steps to run this step on all data at once you can run the command below.

```
step3_compute_bcm(1, 25);
``` 

Folder ‘data_step3’ will be made in the current directory by running the command above on the selected data. The folder contains a sub-folder ‘figures’. 

The V1 ROI vertices are shown in green and high confidence vertices will be marked with blue dots. To compare the number of flipped triangles before and after smoothing, the following figures will be provided in folder 'figures' which is generated after running the commands above.

![step3-fig2](https://user-images.githubusercontent.com/97193844/176073033-517b8b28-2d7a-4d93-b202-6b8ef694a943.png)

![Step3-figure1](https://user-images.githubusercontent.com/97193844/176070781-802bca4e-1a46-4994-83b3-7a683eca110f.png)

Also you can find the Beltrami coefficient map of the unsmoothed and smoothed visual field coordinates after removing the flipped triangles.

![step3-fig3](https://user-images.githubusercontent.com/97193844/176450879-a436ae80-31e0-43f9-8cf3-a56fde1d640b.png)

