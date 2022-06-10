#include "CUDAfunctions.cuh"
// #include "hostfunctions.h"
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
    // Parameters
    const float camExposure = 700.0;
    const float camExpRatio = 0.78;

    const int OffsetX = atoi(argv[1]);
    // const int OffsetX = 592;
    const int OffsetY = atoi(argv[2]);
    // const int OffsetY = 510;
    
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

    cameraSetup(pCam,imgLen,OffsetX,OffsetY,camExposure,camExpRatio);

    getImgAndBundleAdjCheck(pCam,imgLen,gridSize,intrSize,srchSize,zFront,dz,wavLen,dx,blockSize);

    
    pCam[0]->DeInit();
    pCam[1]->DeInit();
    system->ReleaseInstance();

    cudaDeviceReset();
    return 0;
}