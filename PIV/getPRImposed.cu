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
#include <signal.h>

volatile sig_atomic_t e_flag = 0;
void abrt_handler(int sig){
    e_flag = 1;
}

int main(int argc, char** argv){
    std::cout << argv[0] << " Starting..." << std::endl;

    if ( signal(SIGINT, abrt_handler) == SIG_ERR ) {
        exit(1);
    }
    
    // Parameters
    const float camExposure = 100.0;
    const float gainInit = 1.0;

    const int OffsetX = atoi(argv[1]);
    // const int OffsetX = 592;
    const int OffsetY = atoi(argv[2]);
    // const int OffsetY = 514;

    float gain1,gain2;
    std::tie(gain1,gain2) = readGain("./gain.dat");
    
    const int imgLen = 512;
    const int intrSize = imgLen/8;
    const int gridSize = imgLen/8;
    const int srchSize = imgLen/4;
    const int gridNum = (int)(imgLen/gridSize);

    const float zF = 1000*60.0;
    const float dz = 50.0;
    const float waveLen = 0.532;
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
    cameraSetup(pCam,imgLen,OffsetX,OffsetY,camExposure,gain1,gain2);


    // Processing

    // Constant Declaration
    const int datLen = imgLen*2;
    dim3 grid((int)ceil((float)datLen/(float)blockSize),(int)ceil((float)datLen/(float)blockSize)), block(blockSize,blockSize);

    // Propagation Init
    float *d_sqr;
    cufftComplex *d_transF, *d_transInt, *d_transPR, *d_transPRInv;
    CHECK(cudaMalloc((void **)&d_sqr, sizeof(float)*datLen*datLen));
    CHECK(cudaMalloc((void **)&d_transF, sizeof(cufftComplex)*datLen*datLen));
    CHECK(cudaMalloc((void **)&d_transInt, sizeof(cufftComplex)*datLen*datLen));
    CHECK(cudaMalloc((void **)&d_transF2, sizeof(cufftComplex)*datLen*datLen));
    CHECK(cudaMalloc((void **)&d_transF2, sizeof(cufftComplex)*datLen*datLen));
    CuTransSqr<<<grid,block>>>(d_sqr,datLen,waveLen,dx);
    CuTransFunc<<<grid,block>>>(d_transF,d_sqr,zF,waveLen,datLen,dx);
    CuTransFunc<<<grid,block>>>(d_transF2,d_sqr,zF*2.0,waveLen,datLen,dx);
    CuTransFunc<<<grid,block>>>(d_transInt,d_sqr,dz,waveLen,datLen,dx);
    std::cout << "Gabor Init OK" << std::endl;


    camList.Clear();
    system->ReleaseInstance();

    cudaDeviceReset();
    return 0;
}