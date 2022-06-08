/**
 * @file hostfunctions.h
 * @brief Jetson上の位相回復ホログラフィによる流刑分布取得実験用関数群
 * @author Dai Nakai
 * @date May, 2022.
 */

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include "Spinnaker.h"
#include "SpinGenApi/SpinnakerGenApi.h"
// #include <gtk/gtk.h>
#include <string>

/**
 * @fn
 * @brief path の画像を読み込み正規化されたfloat配列 *img に格納する
 * @param path 画像のパス
 * @param imgLen 画像一辺の長さ
 * @param img (imgLen,imgLen)型配列のポインタ
 * @return なし
 */
void getFloatimage(float *img, const int imgLen, const char* path){
    unsigned char *UIntImage;
    UIntImage = (unsigned char *)malloc(sizeof(unsigned char)*imgLen*imgLen);

    FILE *fp;
    fp = fopen(path,"rb");
    if(fp == NULL){
        printf("%s seems not to exist! Quitting...\n",path);
        exit(1);
    }
    fseek(fp,1078,0);
    fread(UIntImage,sizeof(unsigned char),imgLen*imgLen,fp);
    fclose(fp);
    
    for (int i = 0; i < imgLen*imgLen; i++){
        img[i] = (float)UIntImage[i]/255.0;
    }

    free(UIntImage);
}

/**
 * @fn
 * @brief vecArrayをvecArray.datに保存します
 * @param vecArrayX ベクトル場配列のポインタ
 * @param imgLen 画像一辺の長さ
 * @param gridNum gridNum-1 で各次元のベクトル個数
 * @return なし
 */
void saveVecArray(float *vecArrayX, float *vecArrayY, const int gridSize, const int gridNum){
    FILE *fp;
    if ((fp = fopen("vecArray.dat", "w")) == NULL) {
	    printf("File access not available! Quitting...!\n");
	    exit(1);
    }

    int n = gridNum-1;

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            fprintf(fp,"%d %d %lf %lf\n",(j+1)*gridSize,(i+1)*gridSize,vecArrayX[i*n+j],vecArrayY[i*n+j]);
        }
    }

    fclose(fp);
}

/**
 * @fn
 * @brief vecArrayをベクトル場表示します。
 * @param imgLen 画像一辺の長さ
 * @return なし
 */
void plotVecFieldOnGnuplot(const int imgLen){
    FILE *gp;
    const char* outputPath = "vecArrayPlot.pdf";
    const char* vecArrayDataPath = "vecArray.dat";
    const char* xxlabel = "{/Times-New-Roman:Italic=20 x} [pixel]";
    const char* yylabel = "{/Times-New-Roman:Italic=20 y} [pixel]";

    const int vecLenSclr = 20;
    
    if ((gp = popen("gnuplot", "w")) == NULL) {
	    printf("gnuplot is not available! Quitting...!\n");
	    exit(1);
    }

    fprintf(gp,"set terminal pdfcairo enhanced font 'Times New Roman,15' \n");
	fprintf(gp,"set output '%s'\n",outputPath);
	fprintf(gp,"set size ratio 1\n");
    fprintf(gp,"set xrange[0:%d]\n",imgLen);
    fprintf(gp,"set yrange[0:%d]\n",imgLen);
	fprintf(gp,"set palette rgb 33,13,10\n");

  // fprintf(gp,"set yrange reverse\n");

	fprintf(gp,"set xlabel '%s'offset 0.0,0.5\n",xxlabel);
	fprintf(gp,"set ylabel '%s'offset 0.5,0.0\n",yylabel);

	fprintf(gp,"plot '%s' using 1:2:(%d*$3):(%d*$4):(sqrt($3*$3+$4*$4))  w vector lc palette ti ''\n",vecArrayDataPath,vecLenSclr,vecLenSclr);

 	fflush(gp); //Clean up Data

	fprintf(gp, "exit\n"); // Quit gnuplot
	pclose(gp);
}

/**
 * @fn
 * @brief 2カメラの設定
 * @param cam2OffSetX カメラ2のオフセット．事前にSpinViewで調整
 * @return なし
 */
void cameraSetup(Spinnaker::CameraPtr pCam[2], int imgLen, int cam2OffSetX, int cam2OffSetY, float exposure=400.0, float expratio=0.78, float gain=0.0){
    // unsigned int numCameras = camList.GetSize();
    // if (numCameras==0){
    //     printf("No Cameras are Connected! Quitting...\n");
    //     exit(1);
    // }

    // Spinnaker::CameraPtr pCam[numCameras];
    // // printf("CameraNumber\tModelName\tSerialNumber\n");
    // std::cout << "Camera" << "\t" << "ModelName" << "\t\t\t" << "SerialNumber" << std::endl;
    // for (int i = 0; i < numCameras; i++){
    //     pCam[i] = camList.GetByIndex(i);
    //     pCam[i]->Init();
    //     Spinnaker::GenICam::gcstring modelName = pCam[i]->TLDevice.DeviceModelName.GetValue();
    //     Spinnaker::GenICam::gcstring serialNum = pCam[i]->TLDevice.DeviceSerialNumber.GetValue();
    //     // printf("%d\t%s\t%s\n",i,modelName,serialNum);
    //     std::cout << i << "\t" << modelName << "\t" << serialNum << std::endl;
    // }
    // printf("\n");

    // Settings common to all cameras
    for (int i = 0; i < 2; i++){
        pCam[i]->GammaEnable.SetValue(false);
        pCam[i]->AdcBitDepth.SetValue(Spinnaker::AdcBitDepth_Bit12);
        pCam[i]->AcquisitionFrameRateEnable.SetValue(false);
        pCam[i]->PixelFormat.SetValue(Spinnaker::PixelFormat_Mono16);
        // pCam[i]->Width.SetValue(imgLen);
        // pCam[i]->Height.SetValue(imgLen);
        pCam[i]->ExposureAuto.SetValue(Spinnaker::ExposureAutoEnums::ExposureAuto_Off);
        pCam[i]->ExposureMode.SetValue(Spinnaker::ExposureModeEnums::ExposureMode_Timed);
        pCam[i]->GainAuto.SetValue(Spinnaker::GainAutoEnums::GainAuto_Off);
        pCam[i]->Gain.SetValue(gain);
        pCam[i]->AcquisitionMode.SetValue(Spinnaker::AcquisitionModeEnums::AcquisitionMode_Continuous);
    }



    // Settings for Camera 1
    pCam[0]->OffsetX.SetValue((int)((2048-imgLen)/2));
    pCam[0]->OffsetY.SetValue((int)((1536-imgLen)/2));
    pCam[0]->Width.SetValue(imgLen);
    pCam[0]->Height.SetValue(imgLen);
    pCam[0]->ExposureTime.SetValue(exposure);
    pCam[0]->ReverseX.SetValue(false);
    pCam[0]->ReverseY.SetValue(false);
    pCam[0]->TriggerMode.SetValue(Spinnaker::TriggerModeEnums::TriggerMode_On);
    pCam[0]->TriggerSource.SetValue(Spinnaker::TriggerSourceEnums::TriggerSource_Software);
    pCam[0]->TriggerSelector.SetValue(Spinnaker::TriggerSelectorEnums::TriggerSelector_FrameStart);
    pCam[0]->LineSelector.SetValue(Spinnaker::LineSelectorEnums::LineSelector_Line1);
    pCam[0]->LineMode.SetValue(Spinnaker::LineModeEnums::LineMode_Output);
    pCam[0]->LineSelector.SetValue(Spinnaker::LineSelectorEnums::LineSelector_Line2);
    pCam[0]->V3_3Enable.SetValue(true);
    pCam[0]->TriggerOverlap.SetValue(Spinnaker::TriggerOverlapEnums::TriggerOverlap_Off);

    // Settings for Camera 2
    pCam[1]->OffsetX.SetValue(cam2OffSetX);
    pCam[1]->OffsetY.SetValue(cam2OffSetY);
    pCam[1]->Width.SetValue(imgLen);
    pCam[1]->Height.SetValue(imgLen);
    pCam[1]->ExposureTime.SetValue(exposure*expratio);
    pCam[1]->ReverseX.SetValue(false);
    pCam[1]->ReverseY.SetValue(true);
    pCam[1]->TriggerMode.SetValue(Spinnaker::TriggerModeEnums::TriggerMode_On);
    pCam[1]->TriggerSource.SetValue(Spinnaker::TriggerSourceEnums::TriggerSource_Line3);
    pCam[1]->TriggerSelector.SetValue(Spinnaker::TriggerSelectorEnums::TriggerSelector_FrameStart);
    pCam[1]->TriggerOverlap.SetValue(Spinnaker::TriggerOverlapEnums::TriggerOverlap_ReadOut);

    printf("Camera Setup Completed.\n\n");
}

// /**
//  * @fn
//  * @brief ボタンが押されたときに呼び出されるコールバック関数。
//  * @return なし
//  */
// void clicked_button(GtkWidget* widget, gpointer data)
// {
//     int pn = GPOINTER_TO_INT(data);
//     pn = 0;
//     // exit(0); // プログラムが終了する
//     // gtk_main_quit(); // プログラムが終了する
// }