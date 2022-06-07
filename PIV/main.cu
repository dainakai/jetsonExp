#include "CUDAfunctions.cuh"
#include "hostfunctions.h"
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <cuda_runtime.h>
#include <cufft.h>
#include "Spinnaker.h"
#include "SpinGenApi/SpinnakerGenApi.h"
// #include <opencv2/core.hpp>
// #include <opencv2/imgcodecs.hpp>
// #include <opencv2/highgui.hpp>
// #include <opencv2/core/types_c.h>
// #include <gtk/gtk.h>

int main(int argc, char** argv){
    std::cout << argv[0] << " Starting..." << std::endl;
    // printf("%s Starting...\n", argv[0]);

    // const char* pathImg1 = "./cam1.bmp";
    // const char* pathImg2 = "./cam2.bmp";

    // const int imgLen = 1024;
    // const int intrSize = imgLen/8;
    // const int gridSize = imgLen/8;
    // const int srchSize = imgLen/4;
    // const int gridNum = (int)(imgLen/gridSize);

    // const int blockSize = 16;

    // unsigned char *UIntimage1, *UIntimage2;
    // float *fimg1, *fimg2;
    // UIntimage1 = (unsigned char*)malloc(sizeof(unsigned char)*imgLen*imgLen);
    // UIntimage2 = (unsigned char*)malloc(sizeof(unsigned char)*imgLen*imgLen);
    // fimg1 = (float*)malloc(sizeof(float)*imgLen*imgLen);
    // fimg2 = (float*)malloc(sizeof(float)*imgLen*imgLen);

    // getFloatimage(fimg1,imgLen,pathImg1);
    // getFloatimage(fimg2,imgLen,pathImg2);

    // float vecArrayX[(gridNum-1)*(gridNum-1)];
    // float vecArrayY[(gridNum-1)*(gridNum-1)];
    // float *pvecArrX = (float *)vecArrayX;
    // float *pvecArrY = (float *)vecArrayY;
    // getPIVMapOnGPU(pvecArrX,pvecArrY,fimg1,fimg2,imgLen,gridSize,intrSize,srchSize,blockSize);

    // for (int i = 0; i < gridNum-1; i++)
    // {
    //     for (int j = 0; j < gridNum-1; j++)
    //     {
    //         printf("vx: %f\t vy: %f\n",vecArrayX[i*(gridNum-1) +j],vecArrayY[i*(gridNum-1) +j]);
    //     }
    //     printf("\n");
    // }

    // saveVecArray(vecArrayX,vecArrayY,gridSize,gridNum);
    // plotVecFieldOnGnuplot(imgLen);
    
    // Parameters
    const float camExposure = 400.0;
    const float camGain = 0.0;
    
    const int imgLen = 1024;
    const int intrSize = imgLen/8;
    const int gridSize = imgLen/8;
    const int srchSize = imgLen/4;
    const int gridNum = (int)(imgLen/gridSize);

    const float zFront = 1000*100.0;
    const float dz = 50.0;
    const float wavLen = 0.532;
    const float dx = 3.45/0.5;

    const int blockSize = 16; 


    // Camera Init
    Spinnaker::SystemPtr system = Spinnaker::System::GetInstance();
    Spinnaker::CameraList camList = system->GetCameras();
    unsigned int numCameras = camList.GetSize();
    if (numCameras==0){
        // printf("No Cameras are Connected! Quitting...\n");
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
        // printf("%d\t%s\t%s\n",i,modelName,serialNum);
        std::cout << i << "\t" << modelName << "\t" << serialNum << std::endl;
    }
    if (numCameras != 2){
        // printf("Number of Connected Cameras is not 2. Quitting...\n");
        std::cout << "Number of Connected Cameras is not 2. Quitting..." << std::endl;
        exit(0);
    }
    // printf("\n");
    std::cout << "Camera Enum OK" << std::endl;

    cameraSetup(pCam,1024,100,100);
    getImgAndPIV(pCam,imgLen,gridSize,intrSize,srchSize,zFront,dz,wavLen,dx,blockSize);
    
    pCam[0]->DeInit();
    pCam[1]->DeInit();
    system->ReleaseInstance();

    // GtkWidget* window;
    // gtk_init(&argc,&argv);
    // window = gtk_window_new(GTK_WINDOW_TOPLEVEL);
    // gtk_widget_set_size_request(window,100,50);
    // gtk_window_set_title(GTK_WINDOW(window),"Controller");
    // GtkWidget* button;
    // button = gtk_button_new_with_label("Stop");
    // gtk_container_add(GTK_CONTAINER(window),button);
    // int state = 1;
    // gpointer ptr = GINT_TO_POINTER(state);
    // g_signal_connect(button,"clicked",G_CALLBACK(clicked_button),ptr);
    // gtk_widget_show_all(window);
    // gtk_main();
    // // sleep(1);

    // int i = 0;
    // while (state){
    //     std::cout << i << std::endl;
    //     i += 1;
    // }

    // free(UIntimage1);
    // free(UIntimage2);
    // free(fimg1);
    // free(fimg2);
    cudaDeviceReset();
    return 0;
}