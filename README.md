This repo is to store the matlab code I use to run the SPAD512 camera as the spad_cam class. 

Extra code/functions are temporary for testing and development. 

Usage examples will be uploaded in the future. 

Current functionality is: 

1. Connecting to the camera via the tcp port (CSPAD software should be running)
2. Define acquisition parameters
3. Take an aquisition

Currently working on: 

1. Integrating data read functions into the class
2. Characterising the read/write speeds
3. Integrating correlation measurements

Long-term goals: 

1. Integrated ROI selection capability
2. Binning
3. Correlations between arbitrarily defined ROIs
4. Live-streaming of intensity into matlab
5. Live-streaming of correlations into matlab
