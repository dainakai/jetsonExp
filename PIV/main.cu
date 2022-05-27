#include "CUDAfunctions.cuh"
#include "hostfunctions.h"
#include <stdlib.h>
#include <stdio.h>
#include <cufft.h>

int main(int argc, char** argv){
    printf("%s Starting...\n", argv[0]);

    const char* pathImg1 = "./cam1.bmp";
    const char* pathImg2 = "./cam2.bmp";

    const int imgLen = 512;
    const int intrSize = 128;
    const int gridSize = 128;
    const int srchSize = 256;
    const int gridNum = (int)(imgLen/gridSize);

    const int blockSize = 32;

    unsigned char *UIntimage1, *UIntimage2;
    float *fimg1, *fimg2;
    UIntimage1 = (unsigned char*)malloc(sizeof(unsigned char)*imgLen*imgLen);
    UIntimage2 = (unsigned char*)malloc(sizeof(unsigned char)*imgLen*imgLen);
    fimg1 = (float*)malloc(sizeof(float)*imgLen*imgLen);
    fimg2 = (float*)malloc(sizeof(float)*imgLen*imgLen);

    getFloatimage(fimg1,imgLen,pathImg1);
    getFloatimage(fimg2,imgLen,pathImg2);
    printf("%f\n",fimg1[imgLen*(imgLen-1)]);

    float vecArrayX[(gridNum-1)*(gridNum-1)];
    float vecArrayY[(gridNum-1)*(gridNum-1)];
    float *pvecArrX = (float *)vecArrayX;
    float *pvecArrY = (float *)vecArrayY;
    getPIVMapOnGPU(pvecArrX,pvecArrY,fimg1,fimg2,imgLen,gridSize,intrSize,srchSize,blockSize);

    for (int i = 0; i < gridNum-1; i++)
    {
        for (int j = 0; j < gridNum-1; j++)
        {
            printf("vx: %f\t vy: %f\n",vecArrayX[i*(gridNum-1) +j],vecArrayY[i*(gridNum-1) +j]);
        }
        printf("\n");
    }

    free(UIntimage1);
    free(UIntimage2);
    free(fimg1);
    free(fimg2);
    return 0;
}