# Anomaly Splitting in Gravity Method
## Project Target:
In gravity method which is a geophysical survey, Anomaly splitting operations for a slice that taken from bouguer anomaly map of UK. 
## Required Data for satisfy project target: 
for each point: Longitude, latitude, elevation (m), bouguer anomaly values.
## Desired Outputs: 
Observation of residual and regional anomally values that are splitted from bouger anomaly of slice.
## Project Phases:
  ##  1) DATA COLLECTION PHASE
  - An open source dataset which includes Longitude, latitude, elevation (m), bouguer anomaly values of UK is used. source (csv file)  https://dataunderground.org/dataset/great-britain-land-gravity-survey  

  ##  2) DATA PREPERATION PHASE
  - Dataset is imported as dataframe and data type is checked for mathematical operations.
  ##  3) DATA PROCESSING AND VISUALIZATION PHASE
- Elevation and Bouguer maps are created. <br/>
  ![output1](https://user-images.githubusercontent.com/114949587/225301124-043626e5-433d-4b01-88fe-e5de2093ca79.png)
- Slice is choosed. <br/>
  ![output2](https://user-images.githubusercontent.com/114949587/225301275-d77d17a7-03c0-4daf-b0d1-c40ecc085986.png)  
- CBA values vs coordinates are observed. <br/>
  ![output3](https://user-images.githubusercontent.com/114949587/225303231-f3bcbeef-7757-43c9-98c2-efcd3595bb6f.png)
- New coordinate plane is defined for slice.<br/>
  ![output4](https://user-images.githubusercontent.com/114949587/225303300-d88dab64-55c0-4952-b2ed-3ea51aa28021.png)
- Data of slice that is in length domain turned into wavenumber domain with the help of the Fourier transform.<br/>
  ![output5](https://user-images.githubusercontent.com/114949587/225303340-06bb5779-9aac-4452-8456-6247d02445e8.png)
- Half number of wavenumber is taken.<br/>
  ![output6](https://user-images.githubusercontent.com/114949587/225303384-937f0271-8958-4488-8372-d9056fc07fbd.png)
- Manipulated data is fitted with 2 lines.<br/>
  ![output7](https://user-images.githubusercontent.com/114949587/225303432-18d9a63a-65c9-4ba7-9fd1-cea8823dbd76.png)
- Splitted data is made and desired output is created.<br/>
  ![output8](https://user-images.githubusercontent.com/114949587/225303483-a80f7290-3695-43df-bad1-195c3a9309af.png)
