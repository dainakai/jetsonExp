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
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/device_ptr.h>
#include <thrust/reduce.h>
#include "hostfunctions.h"

/** @def
 * CUDA組み込み関数のチェックマクロ。cudaMalloc や cudaMemcpy に。
 */
#define CHECK(call)                                                             \
{                                                                               \
    const cudaError_t error = call;                                             \
    if(error != cudaSuccess){                                                   \
        std::cout << "Error: " << __FILE__ << ":" << __LINE__ << ", ";          \
        std::cout << "code:" << error << ", reason: " << cudaGetErrorString(error) << std::endl;\
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

        float tmp[3][3];
        for (int i = 0; i < 3; i++){
            for (int j = 0; j < 3; j++){
                tmp[i][j] = 1.0;
            }
        }
        
        if (gridIdxx==0){
            tmp[0][0]=0.0;
            tmp[1][0]=0.0;
            tmp[2][0]=0.0;
        }else if (gridIdxx==(gridNum-2)){
            tmp[0][2]=0.0;
            tmp[1][2]=0.0;
            tmp[2][2]=0.0;
        }else if (gridIdxy==0){
            tmp[0][0]=0.0;
            tmp[0][1]=0.0;
            tmp[0][2]=0.0;
        }else if (gridIdxy==(gridNum-2)){
            tmp[2][0]=0.0;
            tmp[2][1]=0.0;
            tmp[2][2]=0.0;
        }

        corArrayOut[y*corArrSize*(gridNum-1)+x] = 0.0;

        for (int i = 0; i < 3; i++){
            for (int j = 0; j < 3; j++){
                if (tmp[i][j]==1.0){
                    corArrayOut[y*corArrSize*(gridNum-1)+x]+=corArrayIn[y*corArrSize*(gridNum-1+(i-1))+x+(j-1)*corArrSize];
                }
            }
        }

        // Four Corners
        // if ((gridIdxx==0 || gridIdxx==(gridNum-2))&&(gridIdxy==0 || gridIdxy==(gridNum-2))){
        //     corArrayOut[y*corArrSize*(gridNum-1)+x] = corArrayIn[y*corArrSize*(gridNum-1)+x];
        // }else if (gridIdxx==0 || gridIdxx==(gridNum-2)){
        //     corArrayOut[y*corArrSize*(gridNum-1)+x] = 0.0;
        //     corArrayOut[y*corArrSize*(gridNum-1)+x] += corArrayIn[y*corArrSize*(gridNum-1)+x-corArrSize]/3.0;
        //     corArrayOut[y*corArrSize*(gridNum-1)+x] += corArrayIn[y*corArrSize*(gridNum-1)+x+corArrSize]/3.0;
        //     corArrayOut[y*corArrSize*(gridNum-1)+x] += corArrayIn[y*corArrSize*(gridNum-1)+x]/3.0;
        // }else if (gridIdxy==0 || gridIdxy==(gridNum-2)){
        //     corArrayOut[y*corArrSize*(gridNum-1)+x] = 0.0;
        //     corArrayOut[y*corArrSize*(gridNum-1)+x] += corArrayIn[y*corArrSize*(gridNum-1)+x]/3.0;
        //     corArrayOut[y*corArrSize*(gridNum-1)+x] += corArrayIn[y*corArrSize*(gridNum-1+1)+x]/3.0;
        //     corArrayOut[y*corArrSize*(gridNum-1)+x] += corArrayIn[y*corArrSize*(gridNum-1-1)+x]/3.0;
        // }else{
        //     corArrayOut[y*corArrSize*(gridNum-1)+x] = 0.0;
        //     corArrayOut[y*corArrSize*(gridNum-1)+x] += corArrayIn[y*corArrSize*(gridNum-1)+x]/9.0;
        //     corArrayOut[y*corArrSize*(gridNum-1)+x] += corArrayIn[y*corArrSize*(gridNum-1-1)+x]/9.0;
        //     corArrayOut[y*corArrSize*(gridNum-1)+x] += corArrayIn[y*corArrSize*(gridNum-1+1)+x]/9.0;
        //     corArrayOut[y*corArrSize*(gridNum-1)+x] += corArrayIn[y*corArrSize*(gridNum-1)+x-corArrSize]/9.0;
        //     corArrayOut[y*corArrSize*(gridNum-1)+x] += corArrayIn[y*corArrSize*(gridNum-1)+x+corArrSize]/9.0;
        //     corArrayOut[y*corArrSize*(gridNum-1)+x] += corArrayIn[y*corArrSize*(gridNum-1-1)+x-corArrSize]/9.0;
        //     corArrayOut[y*corArrSize*(gridNum-1)+x] += corArrayIn[y*corArrSize*(gridNum-1-1)+x+corArrSize]/9.0;
        //     corArrayOut[y*corArrSize*(gridNum-1)+x] += corArrayIn[y*corArrSize*(gridNum-1+1)+x-corArrSize]/9.0;
        //     corArrayOut[y*corArrSize*(gridNum-1)+x] += corArrayIn[y*corArrSize*(gridNum-1+1)+x+corArrSize]/9.0;
        // }
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

__global__ void CuCharToNormFloatArr(float *out, char16_t *in, int datLen, float Norm){
    // dim3 grid((width+31)/32, (height+31)/32), block(32,32)
    // CHECH deviceQuery and make sure threads per block are 1024!!!!

	int x = blockIdx.x*blockDim.x + threadIdx.x;
    int y = blockIdx.y*blockDim.y + threadIdx.y;

    if( (x < datLen) && (y < datLen) ){
        out[y*datLen + x] = (float)((int)in[y*datLen + x])/Norm;
        // printf("%lf ",out[y*datLen+x]);
    }
}

__global__ void CuNormFloatArrToChar(unsigned char *out, float *in, int datLen, float Norm){
    // dim3 grid((width+31)/32, (height+31)/32), block(32,32)
    // CHECH deviceQuery and make sure threads per block are 1024!!!!

	int x = blockIdx.x*blockDim.x + threadIdx.x;
    int y = blockIdx.y*blockDim.y + threadIdx.y;

    if( (x < datLen) && (y < datLen) ){
        out[y*datLen + x] = (unsigned char)(in[y*datLen + x]*Norm);
        // printf("%lf ",out[y*datLen+x]);
    }
}

__global__ void CuSetArrayCenterHalf(cufftComplex *out, float *img, int imgLen){
    // dim3 grid((width+31)/32, (height+31)/32), block(32,32)
    // CHECH deviceQuery and make sure threads per block are 1024!!!!

	int x = blockIdx.x*blockDim.x + threadIdx.x;
    int y = blockIdx.y*blockDim.y + threadIdx.y;

    if( (x < imgLen) && (y < imgLen) ){
        out[(y+imgLen/2)*imgLen*2 + (x+imgLen/2)].x = img[y*imgLen+x]; 
        out[(y+imgLen/2)*imgLen*2 + (x+imgLen/2)].y = 0.0; 
        // printf("%lf ",out[(y+imgLen/2)*imgLen*2 + (x+imgLen/2)].x);
    }
}

__global__ void CuFFTshift(cufftComplex *data, int datLen){
    // dim3 grid((width+31)/32, (height+31)/32), block(32,32)
    // CHECH deviceQuery and make sure threads per block are 1024!!!!

	int x = blockIdx.x*blockDim.x + threadIdx.x;
    int y = blockIdx.y*blockDim.y + threadIdx.y;
    cufftComplex temp1,temp2;
    
    if((x < datLen/2) && (y < datLen/2)){
        temp1 = data[x + datLen*y];
        data[x + datLen*y] = data[x + datLen/2 + datLen*(y + datLen/2)];
        data[x + datLen/2 + datLen*(y + datLen/2)] = temp1;
    }
    if((x < datLen/2) && (y >= datLen/2)){
        temp2 = data[x + datLen*y];
        data[x + datLen*y] = data[x + datLen/2 + datLen*(y - datLen/2)];
        data[x + datLen/2 + datLen*(y - datLen/2)] = temp2;
    }
}

__global__ void CuComplexMul(cufftComplex *out, cufftComplex *inA, cufftComplex *inB, int datLen){
    // dim3 grid((width+31)/32, (height+31)/32), block(32,32)
    // CHECH deviceQuery and make sure threads per block are 1024!!!!

	int x = blockIdx.x*blockDim.x + threadIdx.x;
    int y = blockIdx.y*blockDim.y + threadIdx.y;
    cufftComplex tmp1, tmp2;

    if( (x < datLen) && (y < datLen) ){
        tmp1 = inA[y*datLen + x];
        tmp2 = inB[y*datLen + x];
        out[y*datLen + x].x = tmp1.x * tmp2.x - tmp1.y * tmp2.y;
        out[y*datLen + x].y = tmp1.x * tmp2.y + tmp1.y * tmp2.x;
    }
}

__global__ void CuGetAbsFromComp(float *out, cufftComplex *in, int datLen){
	int x = blockIdx.x*blockDim.x + threadIdx.x;
    int y = blockIdx.y*blockDim.y + threadIdx.y;

    if( (x < datLen) && (y < datLen) ){
        cufftComplex tmp = in[y*datLen + x];
        out[y*datLen + x] = sqrt(tmp.x * tmp.x + tmp.y * tmp.y); // Need Sqrt() ?
    }
}

__global__ void CuUpdateImposed(float *imp, float *tmp, int datLen){
	int x = blockIdx.x*blockDim.x + threadIdx.x;
    int y = blockIdx.y*blockDim.y + threadIdx.y;

    if( (x < datLen) && (y < datLen) ){
        if (tmp[y*datLen+x] < imp[y*datLen+x]){
            imp[y*datLen+x] = tmp[y*datLen+x];
        }
    }
}

__global__ void CuGetCenterHalf(float *out, float *in, int imgLen){
	int x = blockIdx.x*blockDim.x + threadIdx.x;
    int y = blockIdx.y*blockDim.y + threadIdx.y;

    if( (x < imgLen) && (y < imgLen) ){
        out[y*imgLen + x] = in[(y+imgLen/2)*imgLen*2 + x+imgLen/2];
    }
}

__global__ void CuFillArrayComp(cufftComplex* array, float value, int datLen){
	int x = blockIdx.x*blockDim.x + threadIdx.x;
    int y = blockIdx.y*blockDim.y + threadIdx.y;

    if( (x < datLen) && (y < datLen) ){
        array[y*datLen + x].x = value;
        array[y*datLen + x].y = 0.0;
        // printf("%lf ",array[y*datLen+x].x);
    }
}

__global__ void CuFillArrayFloat(float* array, float value, int datLen){
	int x = blockIdx.x*blockDim.x + threadIdx.x;
    int y = blockIdx.y*blockDim.y + threadIdx.y;

    if( (x < datLen) && (y < datLen) ){
        array[y*datLen + x] = value;
    }
}

__global__ void CuInvFFTDiv(cufftComplex* array, float div, int datLen){
	int x = blockIdx.x*blockDim.x + threadIdx.x;
    int y = blockIdx.y*blockDim.y + threadIdx.y;

    if( (x < datLen) && (y < datLen) ){
        array[y*datLen + x].x /= div;
        array[y*datLen + x].y /= div;
    }
}

void getGaborImposed(float *floatout, unsigned char *charout, char16_t *in, cufftComplex *transF, cufftComplex *transInt, int imgLen, int loopCount, int blockSize=16){
    int datLen = imgLen*2;
    dim3 gridImgLen((int)ceil((float)imgLen/(float)blockSize), (int)ceil((float)imgLen/(float)blockSize)), block(blockSize,blockSize);
    dim3 gridDatLen((int)ceil((float)datLen/(float)blockSize), (int)ceil((float)datLen/(float)blockSize));
    
    char16_t *dev_in;
    CHECK(cudaMalloc((void**)&dev_in,sizeof(char16_t)*imgLen*imgLen));
    float *dev_img;
    CHECK(cudaMalloc((void**)&dev_img,sizeof(float)*imgLen*imgLen));
    CHECK(cudaMemcpy(dev_in, in, sizeof(char16_t)*imgLen*imgLen, cudaMemcpyHostToDevice));
    CuCharToNormFloatArr<<<gridImgLen,block>>>(dev_img,dev_in,imgLen,65535.0);
    thrust::device_ptr<float> thimg(dev_img);
    float meanImg = thrust::reduce(thimg,thimg+imgLen*imgLen, (float)0.0, thrust::plus<float>());
    meanImg /= (float)(imgLen*imgLen);
    std::cout << "Image mean: " << meanImg << std::endl;

    cufftComplex *dev_holo;
    CHECK(cudaMalloc((void**)&dev_holo,sizeof(cufftComplex)*datLen*datLen));
    CuFillArrayComp<<<gridDatLen,block>>>(dev_holo,meanImg,datLen);
    CuSetArrayCenterHalf<<<gridImgLen,block>>>(dev_holo,dev_img,imgLen);

    float *dev_imp;
    CHECK(cudaMalloc((void**)&dev_imp,sizeof(float)*datLen*datLen));
    CuFillArrayFloat<<<gridDatLen,block>>>(dev_imp,1.0,datLen);

    cufftHandle plan;
    cufftPlan2d(&plan, datLen, datLen, CUFFT_C2C);

    cufftExecC2C(plan, dev_holo, dev_holo, CUFFT_FORWARD);
    CuFFTshift<<<gridDatLen,block>>>(dev_holo, datLen);
    CuComplexMul<<<gridDatLen,block>>>(dev_holo, dev_holo, transF, datLen);

    cufftComplex *tmp_holo;
    float *tmp_imp;
    CHECK(cudaMalloc((void**)&tmp_holo,sizeof(cufftComplex)*datLen*datLen));
    CHECK(cudaMalloc((void**)&tmp_imp,sizeof(float)*datLen*datLen));
    
    for (int itr = 0; itr < loopCount; itr++){
        CuComplexMul<<<gridDatLen,block>>>(dev_holo,dev_holo,transInt,datLen);
        CHECK(cudaMemcpy(tmp_holo,dev_holo,sizeof(cufftComplex)*datLen*datLen,cudaMemcpyDeviceToDevice));
        CuFFTshift<<<gridDatLen,block>>>(tmp_holo,datLen);
        cufftExecC2C(plan, tmp_holo, tmp_holo, CUFFT_INVERSE);
        CuInvFFTDiv<<<gridDatLen,block>>>(tmp_holo,(float)(datLen*datLen),datLen);
        CuGetAbsFromComp<<<gridDatLen,block>>>(tmp_imp,tmp_holo,datLen);
        CuUpdateImposed<<<gridDatLen,block>>>(dev_imp,tmp_imp,datLen);
    }

    float *dev_outImp;
    CHECK(cudaMalloc((void**)&dev_outImp,sizeof(float)*imgLen*imgLen));
    CuGetCenterHalf<<<gridImgLen,block>>>(dev_outImp,dev_imp,imgLen);

    CHECK(cudaMemcpy(floatout, dev_outImp, sizeof(float)*imgLen*imgLen, cudaMemcpyDeviceToHost));

    unsigned char *saveImp;
    CHECK(cudaMalloc((void**)&saveImp,sizeof(unsigned char)*imgLen*imgLen));
    CuNormFloatArrToChar<<<gridImgLen,block>>>(saveImp,dev_outImp,imgLen,255.0);

    CHECK(cudaMemcpy(charout, saveImp, sizeof(unsigned char)*imgLen*imgLen, cudaMemcpyDeviceToHost));

    std::cout << (float)floatout[0] << std::endl;
    std::cout << (float)floatout[10] << std::endl;
    std::cout << (float)floatout[100] << std::endl;
    std::cout << (float)floatout[1000] << std::endl;

    cufftDestroy(plan);
    CHECK(cudaFree(dev_in));
    CHECK(cudaFree(dev_img));
    CHECK(cudaFree(dev_holo));
    CHECK(cudaFree(dev_imp));
    CHECK(cudaFree(tmp_holo));
    CHECK(cudaFree(tmp_imp));
    CHECK(cudaFree(dev_outImp));
    CHECK(cudaFree(saveImp));
}

void getImgAndPIV(Spinnaker::CameraPtr pCam[2],const int imgLen, const int gridSize, const int intrSize, const int srchSize, const float zF, const float dz, const float waveLen, const float dx, const int blockSize){
    // Constant Declaretion
    const int datLen = imgLen*2;
    dim3 grid((int)ceil((float)datLen/(float)blockSize),(int)ceil((float)datLen/(float)blockSize)), block(blockSize,blockSize);

    // Gabor Init
    float *d_sqr;
    cufftComplex *d_transF, *d_transF2, *d_transInt;
    CHECK(cudaMalloc((void **)&d_sqr, sizeof(float)*datLen*datLen));
    CHECK(cudaMalloc((void **)&d_transF, sizeof(cufftComplex)*datLen*datLen));
    CHECK(cudaMalloc((void **)&d_transF2, sizeof(cufftComplex)*datLen*datLen));
    CHECK(cudaMalloc((void **)&d_transInt, sizeof(cufftComplex)*datLen*datLen));
    CuTransSqr<<<grid,block>>>(d_sqr,datLen,waveLen,dx);
    CuTransFunc<<<grid,block>>>(d_transF,d_sqr,zF,waveLen,datLen,dx);
    CuTransFunc<<<grid,block>>>(d_transF2,d_sqr,zF*2.0,waveLen,datLen,dx);
    CuTransFunc<<<grid,block>>>(d_transInt,d_sqr,dz,waveLen,datLen,dx);
    std::cout << "Gabor Init OK" << std::endl;

    // Camera Init
    Spinnaker::CameraPtr cam1 = pCam[0];
    Spinnaker::CameraPtr cam2 = pCam[1];
    cam1->BeginAcquisition();
    cam2->BeginAcquisition();
    cam1->TriggerSoftware.Execute();
    Spinnaker::ImagePtr pimg1 = cam1->GetNextImage();
    Spinnaker::ImagePtr pimg2 = cam2->GetNextImage();
    cam1->EndAcquisition();
    cam2->EndAcquisition();
    // unsigned char *charimg1 = (unsigned char *)pimg1->GetData();
    char16_t *charimg1 = (char16_t *)pimg1->GetData();
    // unsigned char *charimg2 = (unsigned char *)pimg2->GetData();
    char16_t *charimg2 = (char16_t *)pimg2->GetData();
    std::cout << (int)charimg1[0] << std::endl;

    float *floatimp1, *floatimp2;
    floatimp1 = (float *)malloc(sizeof(float)*imgLen*imgLen);
    floatimp2 = (float *)malloc(sizeof(float)*imgLen*imgLen);

    unsigned char *charimp1, *charimp2;
    charimp1 = (unsigned char *)malloc(sizeof(unsigned char)*imgLen*imgLen);
    charimp2 = (unsigned char *)malloc(sizeof(unsigned char)*imgLen*imgLen);



    //Save Imposed Image
    Spinnaker::ImagePtr saveImg1 = Spinnaker::Image::Create(imgLen,imgLen,0,0,Spinnaker::PixelFormatEnums::PixelFormat_Mono8,charimp1);
    getGaborImposed(floatimp1,charimp1,charimg1,d_transF2,d_transInt,imgLen,10);
    saveImg1->Convert(Spinnaker::PixelFormat_Mono8);
    
    Spinnaker::ImagePtr saveImg2 = Spinnaker::Image::Create(imgLen,imgLen,0,0,Spinnaker::PixelFormatEnums::PixelFormat_Mono8,charimp2);
    getGaborImposed(floatimp2,charimp2,charimg2,d_transF,d_transInt,imgLen,10);
    saveImg2->Convert(Spinnaker::PixelFormat_Mono8);

    saveImg1->Save("./imposed1.jpg");
    saveImg2->Save("./imposed2.jpg");

    std::cout << "getGaborImposed OK" << std::endl;

    // Original image
    pimg1->Convert(Spinnaker::PixelFormat_Mono8);
    pimg2->Convert(Spinnaker::PixelFormat_Mono8);
    pimg1->Save("./outimg1.jpg");
    pimg2->Save("./outimg2.jpg");

    // PIV
    int gridNum = imgLen/gridSize;
    float vecArrayX[(gridNum-1)*(gridNum-1)];
    float vecArrayY[(gridNum-1)*(gridNum-1)];
    float *pvecArrX = (float *)vecArrayX;
    float *pvecArrY = (float *)vecArrayY;
    getPIVMapOnGPU(pvecArrX,pvecArrY,floatimp1,floatimp2,imgLen,gridSize,intrSize,srchSize,blockSize);
    saveVecArray(pvecArrX,pvecArrY,gridSize,gridNum);
    plotVecFieldOnGnuplot(imgLen);

    // Finalize
    free(floatimp1);
    free(floatimp2);
    free(charimp1);
    free(charimp2);
    CHECK(cudaFree(d_sqr));
    CHECK(cudaFree(d_transF));
    CHECK(cudaFree(d_transF2));
    CHECK(cudaFree(d_transInt));

    // std::cout << "OK" << std::endl;
}

void getNewImage(char16_t *out, char16_t *img, float a[12], int imgLen){
    int bkg = 0;
    for (int j = 0; j < imgLen*imgLen; j++){
        bkg += (int)img[j];
    }
    bkg = (int)(round((float)bkg/(float)(imgLen*imgLen)));


    for (int i = 0; i < imgLen; i++){
        for (int j = 0; j < imgLen; j++){
            int tmpX = (int)(round(a[0]+a[1]*j+a[2]*i+a[3]*j*j+a[4]*i*j+a[5]*i*i));
            int tmpY = (int)(round(a[6]+a[7]*j+a[8]*i+a[9]*j*j+a[10]*i*j+a[11]*i*i));
            if (tmpX>=0 && tmpX<imgLen && tmpY>=0 && tmpY <imgLen){
                out[i*imgLen+j] = img[tmpY*imgLen+tmpX];
            }else{
                out[i*imgLen+j] = bkg;
            }
        }
    }
}

void getImgAndBundleAdjCheck(Spinnaker::CameraPtr pCam[2],const int imgLen, const int gridSize, const int intrSize, const int srchSize, const float zF, const float dz, const float waveLen, const float dx, const int blockSize){
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
    Spinnaker::CameraPtr cam2 = pCam[1];
    cam1->BeginAcquisition();
    cam2->BeginAcquisition();
    cam1->TriggerSoftware.Execute();
    Spinnaker::ImagePtr pimg1 = cam1->GetNextImage();
    Spinnaker::ImagePtr pimg2 = cam2->GetNextImage();
    cam1->EndAcquisition();
    cam2->EndAcquisition();
    // unsigned char *charimg1 = (unsigned char *)pimg1->GetData();
    char16_t *charimg1 = (char16_t *)pimg1->GetData();
    // unsigned char *charimg2 = (unsigned char *)pimg2->GetData();
    char16_t *charimg2 = (char16_t *)pimg2->GetData();

    // Bundle Adj Check
    float coefa[12];
    char *coefPath = "./coefa.dat";
    readCoef(coefPath,coefa);
    char16_t *charimg3;
    charimg3 = (char16_t *)malloc(sizeof(char16_t)*imgLen*imgLen);
    getNewImage(charimg3,charimg2,coefa,imgLen);

    float *floatimp1, *floatimp2;
    floatimp1 = (float *)malloc(sizeof(float)*imgLen*imgLen);
    floatimp2 = (float *)malloc(sizeof(float)*imgLen*imgLen);

    unsigned char *charimp1, *charimp2;
    charimp1 = (unsigned char *)malloc(sizeof(unsigned char)*imgLen*imgLen);
    charimp2 = (unsigned char *)malloc(sizeof(unsigned char)*imgLen*imgLen);

    //Save Imposed Image
    Spinnaker::ImagePtr saveImg1 = Spinnaker::Image::Create(imgLen,imgLen,0,0,Spinnaker::PixelFormatEnums::PixelFormat_Mono8,charimp1);
    getGaborImposed(floatimp1,charimp1,charimg1,d_transF,d_transInt,imgLen,500);
    saveImg1->Convert(Spinnaker::PixelFormat_Mono8);
    
    Spinnaker::ImagePtr saveImg2 = Spinnaker::Image::Create(imgLen,imgLen,0,0,Spinnaker::PixelFormatEnums::PixelFormat_Mono8,charimp2);
    getGaborImposed(floatimp2,charimp2,charimg3,d_transF,d_transInt,imgLen,500);
    saveImg2->Convert(Spinnaker::PixelFormat_Mono8);

    saveImg1->Save("./imposed1.jpg");
    saveImg2->Save("./imposed2.jpg");

    std::cout << "getGaborImposed OK" << std::endl;

    // Original image
    pimg1->Convert(Spinnaker::PixelFormat_Mono8);
    pimg2->Convert(Spinnaker::PixelFormat_Mono8);
    pimg1->Save("./outimg1.jpg");
    pimg1->Save("./outimg2.jpg");

    // PIV
    int gridNum = imgLen/gridSize;
    float vecArrayX[(gridNum-1)*(gridNum-1)];
    float vecArrayY[(gridNum-1)*(gridNum-1)];
    float *pvecArrX = (float *)vecArrayX;
    float *pvecArrY = (float *)vecArrayY;
    getPIVMapOnGPU(pvecArrX,pvecArrY,floatimp1,floatimp2,imgLen,gridSize,intrSize,srchSize,blockSize);
    saveVecArray(pvecArrX,pvecArrY,gridSize,gridNum);
    plotVecFieldOnGnuplot(imgLen);

    // Finalize
    free(floatimp1);
    free(floatimp2);
    free(charimp1);
    free(charimp2);
    free(charimg3);
    CHECK(cudaFree(d_sqr));
    CHECK(cudaFree(d_transF));
    CHECK(cudaFree(d_transInt));
}