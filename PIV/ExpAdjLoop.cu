#include "CUDAfunctions.cuh"
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <cuda_runtime.h>
#include <cufft.h>
#include <string>
#include "Spinnaker.h"
#include "SpinGenApi/SpinnakerGenApi.h"
#include <unistd.h>
#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/core/types_c.h>

int main(int argc, char** argv){
    std::cout << argv[0] << " Starting..." << std::endl;
    
    // Parameters
    const float camExposure = 800.0;
    const float camExpRatio = 0.78;

    const int OffsetX = atoi(argv[1]);
    // const int OffsetX = 584;
    const int OffsetY = atoi(argv[2]);
    // const int OffsetY = 506;
    
    const int imgLen = 512;
    const int intrSize = imgLen/8;
    const int gridSize = imgLen/8;
    const int srchSize = imgLen/4;
    const int gridNum = (int)(imgLen/gridSize);

    const float zFront = 1000*60.0;
    const float dz = 50.0;
    const float wavLen = 0.532;
    const float dx = 3.45/0.5;

    const int blockSize = 16; 


    // Camera Init
    Spinnaker::SystemPtr system = Spinnaker::System::GetInstance();
    Spinnaker::CameraList camList = system->GetCameras();
    unsigned int numCameras = camList.GetSize();
    if (numCameras==0){
        std::cout << "No Cameras are Connected! Quitting..." << std::endl;
        exit(1);
    }
    Spinnaker::CameraPtr pCam[numCameras];
    std::cout << "Camera" << "\t" << "ModelName" << "\t\t\t" << "SerialNumber" << std::endl;
    for (int i = 0; i < numCameras; i++){
        pCam[i] = camList.GetByIndex(i);
        pCam[i]->Init();
        Spinnaker::GenICam::gcstring modelName = pCam[i]->TLDevice.DeviceModelName.GetValue();
        Spinnaker::GenICam::gcstring serialNum = pCam[i]->TLDevice.DeviceSerialNumber.GetValue();
        std::cout << i << "\t" << modelName << "\t" << serialNum << std::endl;
    }
    if (numCameras != 2){
        std::cout << "Number of Connected Cameras is not 2. Quitting..." << std::endl;
        exit(0);
    }
    std::cout << "Camera Enum OK" << std::endl;

    // Camera Setup
    cameraSetup(pCam,imgLen,OffsetX,OffsetY,camExposure,camExpRatio);

    float mean1, mean2, exp1, exp2;
    exp1 = pCam[0]->ExposureTime.GetValue();
    exp2 = pCam[1]->ExposureTime.GetValue();
    // Processing
    while(1){
        std::tie(mean1,mean2) = getCamMean(pCam,imgLen);
        if (abs(mean1-0.5)<=0.01 && abs(mean2-0.5) <= 0.01){
            break;
        }else if(abs(mean2-0.5)<=0.01){
            exp1 += (mean1-0.5)*100.0;
            pCam[0]->ExposureTime.SetValue((double)exp1);
        }else{
            exp2 += (mean2-0.5)*100.0;
            pCam[1]->ExposureTime.SetValue((double)exp2);
        }

    }
    exp1 = pCam[0]->ExposureTime.GetValue();
    exp2 = pCam[1]->ExposureTime.GetValue();
    std::cout << "Cam1 Exposure:" << exp1 << std::endl;
    std::cout << "Cam2 Exposure:" << exp2 << std::endl;
    std::cout << "Exp ratio: " << exp2/exp1 << std::endl;
    
    camList.Clear();
    system->ReleaseInstance();

    cudaDeviceReset();
    return 0;
}