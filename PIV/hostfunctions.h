/**
 * @file hostfunctions.h
 * @brief Jetson上の位相回復ホログラフィによる流刑分布取得実験用関数群
 * @author Dai Nakai
 * @date May, 2022.
 */

#include <stdlib.h>
#include <stdio.h>

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

// /**
//  * @fn
//  * @brief vecArrayをベクトル場表示します。
//  * @param vecArray ベクトル場配列のポインタ
//  * @param imgLen 画像一辺の長さ
//  * @param gridNum gridNum-1 で各次元のベクトル個数
//  * @return なし
//  */
// void plotVecFieldonGnuplot(float *vecArray, const int imgLen, const int gridNum){
//     if ((gp = popen("gnuplot", "w")) == NULL) {
// 	    printf("gnuplot is not available! Quitting...!\n");
// 	    exit(1);
//     }

//     fprintf(gp,"set terminal pngcairo enhanced font 'Times New Roman,15' \n");
// 	fprintf(gp,"set multiplot\n");

// }