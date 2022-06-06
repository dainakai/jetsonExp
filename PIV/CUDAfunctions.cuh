/**
 * @file CUDAfunctions.cuh
 * @brief Jetson上の位相回復ホログラフィによる流刑分布取得実験用カーネル群
 * @author Dai Nakai
 * @date May, 2022.
 */

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <cufft.h>
#include "Spinnaker.h"
#include "SpinGenApi/SpinnakerGenApi.h"
// #include <thrust/host_vector.h>
// #include <thrust/device_vector.h>

/** @def
 * CUDA組み込み関数のチェックマクロ。cudaMalloc や cudaMemcpy に。
 */
#define CHECK(call)                                                             \
{                                                                               \
    const cudaError_t error = call;                                             \
    if(error != cudaSuccess){                                                   \
        printf("Error: %s:%d, ",__FILE__, __LINE__);                            \
        printf("code:%d, reason: %s\n", error, cudaGetErrorString(error));      \
        exit(1);                                                                \
    }                                                                           \
}

__global__ void CuGetCrossCor(float *corArray, float *img1, float *img2, int gridNum, int srchSize, int intrSize, int gridSize, int imgLen){
    int x = blockIdx.x*blockDim.x +threadIdx.x;
    int y = blockIdx.y*blockDim.y +threadIdx.y;

    if (x < (srchSize-intrSize+1)*(gridNum-1) && y < (srchSize-intrSize+1)*(gridNum-1)){
        int gridIdxx = x/(srchSize-intrSize+1);
        int gridIdxy = y/(srchSize-intrSize+1);
        int idxx = x - gridIdxx*(srchSize-intrSize+1);
        int idxy = y - gridIdxy*(srchSize-intrSize+1);

        int a1 = gridIdxy*gridSize + (int)(intrSize/2);
        int a2 = gridIdxx*gridSize + (int)(intrSize/2);
        int b1 = gridIdxy*gridSize + idxy;
        int b2 = gridIdxx*gridSize + idxx;

        float meanA = 0.0;
        float meanB = 0.0;
        float num = 0.0;
        float denomA = 0.0;
        float denomB = 0.0;

        for (int i = 0; i < intrSize; i++){
            for (int j = 0; j < intrSize; j++){
                meanA += img1[(a1+i)*imgLen + a2+j];
                meanB += img2[(b1+i)*imgLen + b2+j];
            }
        }

        meanA /= (float)(intrSize*intrSize);
        meanB /= (float)(intrSize*intrSize);

        for (int i = 0; i < intrSize; i++){
            for (int j = 0; j < intrSize; j++){
                num += (img1[(a1+i)*imgLen + a2+j]-meanA)*(img2[(b1+i)*imgLen + b2+j]-meanB);
                denomA += (img1[(a1+i)*imgLen + a2+j]-meanA)*(img1[(a1+i)*imgLen + a2+j]-meanA);
                denomB += (img2[(b1+i)*imgLen + b2+j]-meanB)*(img2[(b1+i)*imgLen + b2+j]-meanB);
            }
        }

        corArray[y*(srchSize-intrSize+1)*(gridNum-1) + x] = num/(sqrt(denomA)*sqrt(denomB));
    }
}

__global__ void CuGetVector(float *vecArrayX, float *vecArrayY, float *corArray, int gridNum, int corArrSize, int intrSize){
    int gridIdxx = blockIdx.x*blockDim.x +threadIdx.x;
    int gridIdxy = blockIdx.y*blockDim.y +threadIdx.y;

    if (gridIdxx < gridNum-1 && gridIdxy < gridNum -1){
        // printf("%d,%d\n",gridIdxx,gridIdxy);
        int x0 = 0;
        int y0 = 0;
        float tmp = 0.0;

        for (int i = 0; i < corArrSize; i++){
            for (int j = 0; j < corArrSize; j++){
                if (corArray[corArrSize*(gridNum-1)*(corArrSize*gridIdxy+i)+corArrSize*gridIdxx+j]>tmp){
                    x0 = corArrSize*gridIdxx + j;
                    y0 = corArrSize*gridIdxy + i;
                    tmp = corArray[corArrSize*(gridNum-1)*(corArrSize*gridIdxy+i)+corArrSize*gridIdxx+j];
                }
            }
        }

        if (x0==0 || y0==0 || x0==corArrSize*(gridNum-1)-1 || y0==corArrSize*(gridNum-1)-1){
            vecArrayX[(gridNum-1)*gridIdxy+gridIdxx] = (float)x0 - (float)(intrSize)/2.0 - gridIdxx*corArrSize;
            vecArrayY[(gridNum-1)*gridIdxy+gridIdxx] = (float)y0 - (float)(intrSize)/2.0 - gridIdxy*corArrSize;
        }else{
            float valy1x0 = corArray[corArrSize*(gridNum-1)*(y0+1) + x0];
            float valy0x0 = corArray[corArrSize*(gridNum-1)*y0 + x0];
            float valyInv1x0 = corArray[corArrSize*(gridNum-1)*(y0-1) + x0];
            float valy0x1 = corArray[corArrSize*(gridNum-1)*y0 + x0+1];
            float valy0xInv1 = corArray[corArrSize*(gridNum-1)*y0 + x0-1];

            if (valy1x0-2.0*valy0x0+valyInv1x0==0.0 || valy0x1-2.0*valy0x0+valy0xInv1==0.0){
                valy0x0 += 0.00001;
            }

            vecArrayX[(gridNum-1)*gridIdxy+gridIdxx] = (float)x0 - (valy0x1-valy0xInv1)/(valy0x1-2.0*valy0x0+valy0xInv1)/2.0 - (float)(intrSize)/2.0 - gridIdxx*corArrSize;
            vecArrayY[(gridNum-1)*gridIdxy+gridIdxx] = (float)y0 - (valy1x0-valyInv1x0)/(valy1x0-2.0*valy0x0+valyInv1x0)/2.0 - (float)(intrSize)/2.0 - gridIdxy*corArrSize;
        }
    }
}

__global__ void CuErrorCorrect(float *corArrayIn, float *corArrayOut, int corArrSize, int gridNum){
    int x = blockIdx.x*blockDim.x +threadIdx.x;
    int y = blockIdx.y*blockDim.y +threadIdx.y;

    if (x<corArrSize*(gridNum-1) && y<corArrSize*(gridNum-1)){
        int gridIdxx = x/corArrSize;
        int gridIdxy = y/corArrSize;

        // Four Corners
        if ((gridIdxx==0 || gridIdxx==(gridNum-2))&&(gridIdxy==0 || gridIdxy==(gridNum-2))){
            corArrayOut[y*corArrSize*(gridNum-1)+x] = corArrayIn[y*corArrSize*(gridNum-1)+x];
        }else if (gridIdxx==0 || gridIdxx==(gridNum-2)){
            corArrayOut[y*corArrSize*(gridNum-1)+x] = 0.0;
            corArrayOut[y*corArrSize*(gridNum-1)+x] += corArrayIn[y*corArrSize*(gridNum-1)+x-corArrSize]/3.0;
            corArrayOut[y*corArrSize*(gridNum-1)+x] += corArrayIn[y*corArrSize*(gridNum-1)+x+corArrSize]/3.0;
            corArrayOut[y*corArrSize*(gridNum-1)+x] += corArrayIn[y*corArrSize*(gridNum-1)+x]/3.0;
        }else if (gridIdxy==0 || gridIdxy==(gridNum-2)){
            corArrayOut[y*corArrSize*(gridNum-1)+x] = 0.0;
            corArrayOut[y*corArrSize*(gridNum-1)+x] += corArrayIn[y*corArrSize*(gridNum-1)+x]/3.0;
            corArrayOut[y*corArrSize*(gridNum-1)+x] += corArrayIn[y*corArrSize*(gridNum-1+1)+x]/3.0;
            corArrayOut[y*corArrSize*(gridNum-1)+x] += corArrayIn[y*corArrSize*(gridNum-1-1)+x]/3.0;
        }else{
            corArrayOut[y*corArrSize*(gridNum-1)+x] = 0.0;
            corArrayOut[y*corArrSize*(gridNum-1)+x] += corArrayIn[y*corArrSize*(gridNum-1)+x]/9.0;
            corArrayOut[y*corArrSize*(gridNum-1)+x] += corArrayIn[y*corArrSize*(gridNum-1-1)+x]/9.0;
            corArrayOut[y*corArrSize*(gridNum-1)+x] += corArrayIn[y*corArrSize*(gridNum-1+1)+x]/9.0;
            corArrayOut[y*corArrSize*(gridNum-1)+x] += corArrayIn[y*corArrSize*(gridNum-1)+x-corArrSize]/9.0;
            corArrayOut[y*corArrSize*(gridNum-1)+x] += corArrayIn[y*corArrSize*(gridNum-1)+x+corArrSize]/9.0;
            corArrayOut[y*corArrSize*(gridNum-1)+x] += corArrayIn[y*corArrSize*(gridNum-1-1)+x-corArrSize]/9.0;
            corArrayOut[y*corArrSize*(gridNum-1)+x] += corArrayIn[y*corArrSize*(gridNum-1-1)+x+corArrSize]/9.0;
            corArrayOut[y*corArrSize*(gridNum-1)+x] += corArrayIn[y*corArrSize*(gridNum-1+1)+x-corArrSize]/9.0;
            corArrayOut[y*corArrSize*(gridNum-1)+x] += corArrayIn[y*corArrSize*(gridNum-1+1)+x+corArrSize]/9.0;
        }
    }
}

/**
 * @fn
 * @brief ホストメモリ上のfloat画像配列img1,img2に対してPIVを行う。エラーコレクトはパラボラ。Jetson(Maxwell)では blockSize = 32 はダメそう？
 * @param vecArray 出力ベクトル配列のポインタ。((gridNum-1),(gridNum-1),2)型で(:,:,0)はx成分、(:,:,1)はy成分。gridNum = floor(imgLen/gridSize).
 * @param intrSize 参照窓の大きさ
 * @param srchSize 探査窓の大きさ
 * @param gridSize 参照窓の充填グリッド
 * @param blockSize CUDAブロックサイズ
 * @return なし
 */
void getPIVMapOnGPU(float *vecArrayX, float *vecArrayY, float *img1, float *img2, int imgLen, int gridSize, int intrSize, int srchSize, int blockSize){
    const int gridNum = (int)(imgLen/gridSize);
    
    float *dev_corArray, *dev_corArray2, *dev_vecArrayX, *dev_vecArrayY;
    CHECK(cudaMalloc((void **)&dev_corArray, sizeof(float)*(srchSize-intrSize+1)*(gridNum-1)*(srchSize-intrSize+1)*(gridNum-1)));
    CHECK(cudaMalloc((void **)&dev_corArray2, sizeof(float)*(srchSize-intrSize+1)*(gridNum-1)*(srchSize-intrSize+1)*(gridNum-1)));
    CHECK(cudaMalloc((void **)&dev_vecArrayX, sizeof(float)*(gridNum-1)*(gridNum-1)));
    CHECK(cudaMalloc((void **)&dev_vecArrayY, sizeof(float)*(gridNum-1)*(gridNum-1)));

    dim3 grid((int)ceil((float)(srchSize-intrSize+1)*(gridNum-1)/(float)blockSize), (int)ceil((float)(srchSize-intrSize+1)*(gridNum-1)/(float)blockSize)), block(blockSize,blockSize);
    dim3 grid2((int)ceil((float)(gridNum-1)/(float)blockSize), (int)ceil((float)(gridNum-1)/(float)blockSize));
    // printf("%d\n",(int)ceil((float)(gridNum-1)/(float)blockSize));

    float *dev_img1, *dev_img2;
    CHECK(cudaMalloc((void **)&dev_img1, sizeof(float)*imgLen*imgLen));
    CHECK(cudaMalloc((void **)&dev_img2, sizeof(float)*imgLen*imgLen));

    CHECK(cudaMemcpy(dev_img1, img1, sizeof(float)*imgLen*imgLen, cudaMemcpyHostToDevice));
    CHECK(cudaMemcpy(dev_img2, img2, sizeof(float)*imgLen*imgLen, cudaMemcpyHostToDevice));

    CuGetCrossCor<<<grid,block>>>(dev_corArray,dev_img1,dev_img2,gridNum,srchSize,intrSize,gridSize,imgLen);
    CuErrorCorrect<<<grid,block>>>(dev_corArray,dev_corArray2,srchSize-intrSize+1,gridNum);
    CuGetVector<<<grid2,block>>>(dev_vecArrayX,dev_vecArrayY,dev_corArray,gridNum,(srchSize-intrSize+1),intrSize);
    
    CHECK(cudaMemcpy(vecArrayX,dev_vecArrayX,sizeof(float)*(gridNum-1)*(gridNum-1),cudaMemcpyDeviceToHost));
    CHECK(cudaMemcpy(vecArrayY,dev_vecArrayY,sizeof(float)*(gridNum-1)*(gridNum-1),cudaMemcpyDeviceToHost));

    CHECK(cudaFree(dev_corArray));
    CHECK(cudaFree(dev_vecArrayX));
    CHECK(cudaFree(dev_vecArrayY));
    CHECK(cudaFree(dev_img1));
    CHECK(cudaFree(dev_img2));
}

__global__ void CuTransFunc(cufftComplex *output, float *sqr, float trans_z, float waveLen, int datLen, float dx){
    // dim3 grid((width+31)/32, (height+31)/32), block(32,32)
    // CHECH deviceQuery and make sure threads per block are 1024!!!!

	int x = blockIdx.x*blockDim.x + threadIdx.x;
    int y = blockIdx.y*blockDim.y + threadIdx.y;
    float tmp;
    float tmpx, tmpy;
    float uband = 1.0/waveLen/sqrt(2*trans_z/(float)datLen/dx + 1);

    if( (x < datLen) && (y < datLen) ){
        tmp = 2.0*3.14159265358979*trans_z/waveLen*sqrt(sqr[x + datLen*y]);
        output[x + datLen*y].x = cos(tmp);
        output[x + datLen*y].y = sin(tmp);
        tmpx = abs(((float)x - (float)datLen/2.0)*waveLen/(float)datLen/dx);
        tmpy = abs(((float)y - (float)datLen/2.0)*waveLen/(float)datLen/dx);
        if (tmpx > uband || tmpy > uband){
            output[x + datLen*y].x = 0.0;
            output[x + datLen*y].y = 0.0;
        }
    }
}

__global__ void CuTransSqr(float *d_sqr, int datLen, float waveLen, float dx){
    // dim3 grid((width+31)/32, (height+31)/32), block(32,32)
    // CHECH deviceQuery and make sure threads per block are 1024!!!!

	int x = blockIdx.x*blockDim.x + threadIdx.x;
    int y = blockIdx.y*blockDim.y + threadIdx.y;
    float x_f = (float)x;
    float y_f = (float)y;
    float w_f = (float)datLen;

    if( (x < datLen) && (y < datLen) ){
        d_sqr[x + datLen*y] = 1.0 - ((x_f - w_f/2.0)*waveLen/w_f/dx)*((x_f - w_f/2.0)*waveLen/w_f/dx) - ((y_f - w_f/2.0)*waveLen/w_f/dx)*((y_f - w_f/2.0)*waveLen/w_f/dx);
    }
}

void getImgAndPIV(Spinnaker::CameraPtr pCam[2],const int imgLen, const int gridSize, const int intrSize, const int srchSize, const float zF, const float dz, const float waveLen, const float dx, const int blockSize){
    // Constant Declaretion
    const int datLen = imgLen*2;
    dim3 grid((int)ceil((float)datLen/(float)blockSize),(int)ceil((float)datLen/(float)blockSize)), block(blockSize,blockSize);

    // Gabor Init
    float *d_sqr;
    cufftComplex *d_transF, *d_transInt;
    CHECK(cudaMalloc((void **)&d_sqr, sizeof(float)*datLen*datLen));
    CHECK(cudaMalloc((void **)&d_transF, sizeof(cufftComplex)*datLen*datLen));
    CHECK(cudaMalloc((void **)&d_transInt, sizeof(cufftComplex)*datLen*datLen));
    CuTransSqr<<<grid,block>>>(d_sqr,datLen,waveLen,dx);
    CuTransFunc<<<grid,block>>>(d_transF,d_sqr,zF,waveLen,datLen,dx);
    CuTransFunc<<<grid,block>>>(d_transInt,d_sqr,dz,waveLen,datLen,dx);
    std::cout << "Gabor Init OK" << std::endl;

    // Camera Init
    Spinnaker::CameraPtr cam1 = pCam[0];
    // Spinnaker::CameraPtr cam1 = camList.GetByIndex(0);
    Spinnaker::CameraPtr cam2 = pCam[1];
    // Spinnaker::CameraPtr cam2 = camList.GetByIndex(1);
    cam1->BeginAcquisition();
    cam2->BeginAcquisition();
    cam1->TriggerSoftware.Execute();
    Spinnaker::ImagePtr pimg1 = cam1->GetNextImage();
    Spinnaker::ImagePtr pimg2 = cam2->GetNextImage();
    cam1->EndAcquisition();
    cam2->EndAcquisition();
    pimg1->Convert(Spinnaker::PixelFormat_Mono8);
    pimg2->Convert(Spinnaker::PixelFormat_Mono8);
    pimg1->Save("./outimg1.png");
    pimg1->Save("./outimg2.png");
    cam1->DeInit();
    cam2->DeInit();

    // Finalize
    CHECK(cudaFree(d_sqr));
    CHECK(cudaFree(d_transF));
    CHECK(cudaFree(d_transInt));

    std::cout << "OK" << std::endl;
}

// void getGaborImposed(float *out, char *in, cufftComplex *transF, cufftComplex *transInt, int imgLen, int blockSize=16){
//     int datLen = imgLen*2;
//     dim3 grid((int)ceil((float)datLen/(float)blockSize), (int)ceil((float)datLen/(float)blockSize)), block(blockSize,blockSize);
    
//     thrust::device_ptr<char> 
// }