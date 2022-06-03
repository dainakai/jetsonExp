/**
 * @file CUDAfunctions.cuh
 * @brief Jetson上の位相回復ホログラフィによる流刑分布取得実験用カーネル群
 * @author Dai Nakai
 * @date May, 2022.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

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

__global__ void errorCorrect(float *corArrayIn, float *corArrayOut, int corArrSize, int gridNum){
    int x = blockIdx.x*blockDim.x +threadIdx.x;
    int y = blockIdx.y*blockDim.y +threadIdx.y;

    if (x<corArrSize*(gridNum-1) && y<corArrSize*(gridNum-1)){
        
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
    
    float *dev_corArray, *dev_vecArrayX, *dev_vecArrayY;
    CHECK(cudaMalloc((void **)&dev_corArray, sizeof(float)*(srchSize-intrSize+1)*(gridNum-1)*(srchSize-intrSize+1)*(gridNum-1)));
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
    CuGetVector<<<grid2,block>>>(dev_vecArrayX,dev_vecArrayY,dev_corArray,gridNum,(srchSize-intrSize+1),intrSize);
    
    CHECK(cudaMemcpy(vecArrayX,dev_vecArrayX,sizeof(float)*(gridNum-1)*(gridNum-1),cudaMemcpyDeviceToHost));
    CHECK(cudaMemcpy(vecArrayY,dev_vecArrayY,sizeof(float)*(gridNum-1)*(gridNum-1),cudaMemcpyDeviceToHost));

    CHECK(cudaFree(dev_corArray));
    CHECK(cudaFree(dev_vecArrayX));
    CHECK(cudaFree(dev_vecArrayY));
    CHECK(cudaFree(dev_img1));
    CHECK(cudaFree(dev_img2));
}