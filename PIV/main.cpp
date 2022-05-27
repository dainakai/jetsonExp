/**********************************************************************************************
BUNDLE AJUSTMENT
**********************************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#ifdef __linux
 #include <sys/stat.h>
#endif
#ifdef _WIN32
 #include <direct.h>
#endif

/********** variable **********/

const int width = 1024;
const int height = 1024;
const int grid_size_x = 64; //格子サイズ
const int grid_size_y = 64; //格子サイズ
const int interrogation_window_size_x = 128; //参照窓サイズ
const int interrogation_window_size_y = 128; //参照窓サイズ
const int search_window_size_x = 256; //探査窓サイズ
const int search_window_size_y = 256; //探査窓サイズ

const char* input_folder01 = "./";
const char* input_folder02 = "./";
const char* input_data01 = "cam1.bmp";
const char* input_data02 = "cam2.bmp";

const char* output_folder01 = "vector";
const char* output_data_file1 = "vorticity_and_vector.dat";
const char* output_folder02 = "coefficient_a";
const char* output_data_file2 = "coefficient_a.dat";

/********** array **********/

int i,j,k,l,m,n,x,y;
int num;
int x_keep,y_keep,S_sum_x_back,S_sum_x_front,S_sum_y_back,S_sum_y_front;
double I_ave,I_sum,a,b,c,C_keep;

int cam1[height][width];
int cam2[height][width];
const int number_x = width/grid_size_x; //=64 //格子数
const int number_y = height/grid_size_y; //=64 //格子数
int search_window[search_window_size_y][search_window_size_x]; //探査窓配列
int interrogation_window[interrogation_window_size_y][interrogation_window_size_x]; //参照窓配列
int judge_window[interrogation_window_size_y][interrogation_window_size_x]; //審判窓配列
int S_sum[interrogation_window_size_y][interrogation_window_size_x];
double S_ave[interrogation_window_size_y][interrogation_window_size_x];
double C_xy[interrogation_window_size_y][interrogation_window_size_x];

const int N_x = number_x-3;
const int N_y = number_y-3;
const int parameter=6;

double vector_x[N_y][N_x];
double vector_y[N_y][N_x];

double target_x[N_y*N_x],target_y[N_y*N_x]; ///目標の座標
double before_x[N_y*N_x],before_y[N_y*N_x]; ///変換前の座標
double approx_x[N_y*N_x],approx_y[N_y*N_x];
int after_x[height*width],after_y[height*width]; ///変換後の座標
double vec_g[parameter*2],vec_a[parameter*2],delta_a[parameter*2];
double jacobian[2*N_y*N_x][parameter*2];
double mat_H[parameter*2][parameter*2];
double vec_d[parameter*2],vec_id[parameter*2],mat_L[parameter*2][parameter*2],mat_tL[parameter*2][parameter*2];
double vec_e[2*N_y*N_x];
double lld,ld;
double bly,yux;
double error;
double error_x,error_y,keep_error_x,keep_error_y;
int keep_x[height][width],keep_y[height][width];

unsigned char header_buf[1078];
unsigned char image_in01[height][width];
unsigned char image_in02[height][width];
char read_file01[100];
char read_file02[100];
char write_file01[100];
char write_file02[100];

/********** file **********/

FILE *infile_1;
FILE *infile_2;
FILE *outfile_1;
FILE *outfile_2;

/********** MAIN **********/
 int main() {
///linux
    mode_t mode = S_IRUSR | S_IRGRP | S_IXUSR | S_IXGRP | S_IWUSR | S_IWGRP ;
    mkdir(output_folder01,mode);
    mkdir(output_folder02,mode);

///1.画像読み込み***********************************************************************************************
    ///Open PIV images from "img" directry
    sprintf(read_file01,"%s/%s",input_folder01,input_data01);
    infile_1=fopen(read_file01,"rb");
    if(infile_1==NULL){printf("nothing \n");return(0);}

    fread ( header_buf , sizeof ( unsigned char ) ,1078 , infile_1 );
    fread(image_in01 ,sizeof(image_in01),1,infile_1);
    fclose(infile_1);

    sprintf(read_file02,"%s/%s",input_folder02,input_data02);
    infile_2=fopen(read_file02,"rb");
    if(infile_2==NULL){printf("nothing \n");return(0);}

    fread ( header_buf , sizeof ( unsigned char ) ,1078 , infile_2 );
    fread(image_in02 ,sizeof(image_in02),1,infile_2);
    fclose(infile_2);

///2.PIV********************************************************************************************************
    ///2-1.検査窓を動かす*************************************************************************
    for(i=0;i<N_y;i++){
        for(j=0;j<N_x;j++){
            I_ave=0.0;I_sum=0.0;C_keep=0.0;x_keep=0;y_keep=0;

    ///2-2.検査窓の配列を作る*********************************************************************
            for(k=0;k<search_window_size_y;k++){
                for(l=0;l<search_window_size_x;l++){
                    search_window[k][l]=0;
                    search_window[k][l]=image_in02[i*grid_size_y+k][j*grid_size_x+l];
                    }
                }

    ///2-3.参照窓の配列を作る*********************************************************************
            for(k=0;k<interrogation_window_size_y;k++){
                for(l=0;l<interrogation_window_size_x;l++){
                    interrogation_window[k][l]=0;
                    interrogation_window[k][l]=image_in01[(i+1)*grid_size_y+k][(j+1)*grid_size_x+l];
                    }
                }

    ///2-4.参照窓内の平均値を計算する*************************************************************
            I_sum=0.0;I_ave=0.0;
            for(k=0;k<interrogation_window_size_y;k++){
                for(l=0;l<interrogation_window_size_x;l++){
                    I_sum=I_sum+double(interrogation_window[k][l]);
                    }
                }
            I_ave=I_sum/double((search_window_size_x*search_window_size_y));

    ///2-5.審判窓の平均値を計算する***************************************************************
        ///2-5-1.x=0,y=0の合計**********************************************************
            S_sum[0][0]=0;
            for(k=0;k<interrogation_window_size_y;k++){
                for(l=0;l<interrogation_window_size_x;l++){
                    S_sum[0][0]=S_sum[0][0]+search_window[k][l];
                    }
                }

        ///2-5-2.x=0,y=yの合計**********************************************************
            for(y=1;y<interrogation_window_size_y;y++){
                S_sum_y_back=0;
                S_sum_y_front=0;
                for(l=0;l<interrogation_window_size_x;l++){
                    S_sum_y_back=S_sum_y_back+search_window[y-1][l+0];
                    S_sum_y_front=S_sum_y_front+search_window[y+interrogation_window_size_y-1][l+0];
                    }
                S_sum[y][0]=S_sum[y-1][0]-S_sum_y_back+S_sum_y_front;
                }

        ///2-5-3.x=x,y=yの合計**********************************************************
            for(y=0;y<interrogation_window_size_y;y++){
                for(x=1;x<interrogation_window_size_x;x++){
                    S_sum_x_back=0;
                    S_sum_x_front=0;
                    for(k=0;k<interrogation_window_size_y;k++){
                        S_sum_x_back=S_sum_x_back+search_window[y+k][x-1];
                        S_sum_x_front=S_sum_x_front+search_window[y+k][x+interrogation_window_size_x-1];
                        }
                    S_sum[y][x]=S_sum[y][x-1]-S_sum_x_back+S_sum_x_front;
                    }
                }

    ///2-6.審判窓の配列を作る*********************************************************************
            for(y=0;y<interrogation_window_size_y;y++){
                for(x=0;x<interrogation_window_size_x;x++){
                    for(k=0;k<interrogation_window_size_y;k++){
                        for(l=0;l<interrogation_window_size_x;l++){
                            judge_window[k][l]=0;
                            judge_window[k][l]=search_window[y+k][x+l];
                            }
                        }

    ///2-7.相関係数の計算*************************************************************************
        ///2-7-1.審判窓の平均値の計算***************************************************
                    S_ave[y][x]=double(S_sum[y][x])/double(interrogation_window_size_y*interrogation_window_size_x);

        ///2-7-2.相関係数の計算*********************************************************
                    a=0.0;b=0.0;c=0.0;
                    for(k=0;k<interrogation_window_size_y;k++){
                        for(l=0;l<interrogation_window_size_x;l++){
                            a=a+((double(judge_window[k][l])-S_ave[y][x])*(double(interrogation_window[k][l])-I_ave));
                            b=b+((double(judge_window[k][l])-S_ave[y][x])*(double(judge_window[k][l])-S_ave[y][x]));
                            c=c+((double(interrogation_window[k][l])-I_ave)*(double(interrogation_window[k][l])-I_ave));
                            }
                        }
                    if(b*c>0){
                        C_xy[y][x]=a/sqrt(b*c);
                        }
                    else{
                        C_xy[y][x]=0.0;
                        }

    ///2-8.最良値の取り出し**********************************************************************
                    if(C_keep<C_xy[y][x]){
                        C_keep=C_xy[y][x]; //最良相関値
                        x_keep=x; //最良x座標
                        y_keep=y; //最良y座標
                        }
                    }
                }

    ///2-9.3点ピークフィットでベクトル成分を決定*************************************************
            if(C_xy[y_keep+1][x_keep]-2.0*C_xy[y_keep][x_keep]+C_xy[y_keep-1][x_keep]==0.0){
                C_xy[y_keep][x_keep]=C_xy[y_keep][x_keep]+0.00001;
                }
            if(C_xy[y_keep][x_keep+1]-2.0*C_xy[y_keep][x_keep]+C_xy[y_keep][x_keep-1]==0.0){
                C_xy[y_keep][x_keep]=C_xy[y_keep][x_keep]+0.00001;
                }
            vector_y[i][j]=double(y_keep)-double(C_xy[y_keep+1][x_keep]-C_xy[y_keep-1][x_keep])/double(C_xy[y_keep+1][x_keep]-2.0*C_xy[y_keep][x_keep]+C_xy[y_keep-1][x_keep])/2.0-double(grid_size_y);
            vector_x[i][j]=double(x_keep)-double(C_xy[y_keep][x_keep+1]-C_xy[y_keep][x_keep-1])/double(C_xy[y_keep][x_keep+1]-2.0*C_xy[y_keep][x_keep]+C_xy[y_keep][x_keep-1])/2.0-double(grid_size_x);
            printf("%d \t %d \n",i,j);
            }
        }

///3.bandle ajustment*******************************************************************************************
    for(i=0;i<N_y;i++){
        for(j=0;j<N_x;j++){
            target_x[i*N_y+j]=(j+2)*grid_size_x;
            target_y[i*N_y+j]=(i+2)*grid_size_y;
            //printf("%.3lf \n",before_y[i*N+j]);
            }
        }

    for(i=0;i<N_y;i++){
        for(j=0;j<N_x;j++){
            before_x[i*N_y+j]=target_x[i*N_y+j]-vector_x[i][j];
            before_y[i*N_y+j]=target_y[i*N_y+j]-vector_y[i][j];
            //printf("%.3lf \n",before_y[i*N+j]);
            }
        }
/*
    for(i=0;i<N_y;i++){
        for(j=0;j<N_x;j++){
            before_x[i*N_y+j]=target_x[i*N_y+j]+vector_x[i][j];
            before_y[i*N_y+j]=target_y[i*N_y+j]+vector_y[i][j];
            //printf("%.3lf \n",before_y[i*N+j]);
            }
        }
*/
    ///ヤコビ行列
    for(i=0;i<N_y*N_x;i=i+1){
        jacobian[i][0]=-1.0;jacobian[i][1]=-before_x[i];jacobian[i][2]=-before_y[i];jacobian[i][3]=-before_x[i]*before_x[i];jacobian[i][4]=-before_x[i]*before_y[i];jacobian[i][5]=-before_y[i]*before_y[i];
        jacobian[i][6]=0.0;jacobian[i][7]=0.0;jacobian[i][8]=0.;jacobian[i][9]=0.0;jacobian[i][10]=0.0;jacobian[i][11]=0.0;
        jacobian[N_y*N_x+i][0]=0.0;jacobian[N_y*N_x+i][1]=0.0;jacobian[N_y*N_x+i][2]=0.0;jacobian[N_y*N_x+i][3]=0.0;jacobian[N_y*N_x+i][4]=0.0;jacobian[N_y*N_x+i][5]=0.0;
        jacobian[N_y*N_x+i][6]=-1.0;jacobian[N_y*N_x+i][7]=-before_x[i];jacobian[N_y*N_x+i][8]=-before_y[i];jacobian[N_y*N_x+i][9]=-before_x[i]*before_x[i];jacobian[N_y*N_x+i][10]=-before_x[i]*before_y[i];jacobian[N_y*N_x+i][11]=-before_y[i]*before_y[i];
        }
    for(i=0;i<2*N_y*N_x;i++){
        //printf("%d %.3lf %.3lf %.3lf %.3lf %.3lf %.3lf \n",i,jacobian[i][0],jacobian[i][1],jacobian[i][2],jacobian[i][3],jacobian[i][4],jacobian[i][5]);
        }

    ///H行列
    for(i=0;i<parameter*2;i=i+1){
        for(j=0;j<parameter*2;j=j+1){
            mat_H[i][j]=0.0;
            for(k=0;k<2*N_y*N_x;k=k+1){
                mat_H[i][j]+=jacobian[k][i]*jacobian[k][j];
                }
            }
        }
    for(i=0;i<parameter*2;i=i+1){
        for(j=0;j<parameter*2;j=j+1){
            //printf("%d %d %.3lf \n",i,j,mat_H[i][j]);
        }
    }
    for(i=0;i<parameter*2;i=i+1){
        for(j=0;j<parameter*2;j=j+1){
            if(mat_H[i][j]!=mat_H[i][j]){
                printf("%d %d %.3lf %.3lf \n",i,j,mat_H[i][j],mat_H[j][i]);///
            }
        }
    }

    ///修正コレスキー分解
    vec_d[0]=mat_H[0][0];
    mat_L[0][0]=1.0;
    for(i=1;i<parameter*2;i=i+1){
        for(j=0;j<i;j=j+1){
            lld=mat_H[i][j];
            for(k=0;k<j;k=k+1){
                lld-=mat_L[i][k]*mat_L[j][k]*vec_d[k];
                }
            mat_L[i][j]=(1.0/vec_d[j])*lld;
            }
        ld=mat_H[i][i];
        for(k=0;k<i;k=k+1){
            ld-=mat_L[i][k]*mat_L[i][k]*vec_d[k];
            }
        vec_d[i]=ld;
        mat_L[i][i]=1.0;
        }
    for(i=0;i<parameter*2;i=i+1){
        vec_id[i]=1.0/vec_d[i];
        //printf("%.3lf \n",vec_id[i]);
        for(j=0;j<parameter*2;j=j+1){
            //printf("%.3lf \n",mat_L[i][j]);
            }
        }


    ///vec_aの初期値
    for(i=0;i<parameter*2;i=i+1){
        vec_a[i]=1.0;
        }

    num=0;
    do{
    ///近似
    for(i=0;i<N_y*N_x;i=i+1){
        approx_x[i]=vec_a[0]+vec_a[1]*before_x[i]+vec_a[2]*before_y[i]+vec_a[3]*before_x[i]*before_x[i]+vec_a[4]*before_x[i]*before_y[i]+vec_a[5]*before_y[i]*before_y[i];
        approx_y[i]=vec_a[6]+vec_a[7]*before_x[i]+vec_a[8]*before_y[i]+vec_a[9]*before_x[i]*before_x[i]+vec_a[10]*before_x[i]*before_y[i]+vec_a[11]*before_y[i]*before_y[i];
        //printf("%.3lf %.3lf \n",approx_x[i],approx_y[i]);
        }
    ///vec_eの作成
        for(i=0;i<N_y*N_x;i=i+1){
            vec_e[i]=target_x[i]-approx_x[i];
            vec_e[i+N_y*N_x]=target_y[i]-approx_y[i];
            }
        for(i=0;i<2*N_y*N_y;i=i+1){
            //printf("%d %.3lf \n",i,vec_e[i]);
        }

    ///絶対値の計算
        error=0.0;
        for(i=0;i<2*N_y*N_x;i=i+1){
            error+=vec_e[i]*vec_e[i];
            }
    //printf("%d %.3lf \n",num,error);

    ///vec_gの作成
        for(j=0;j<parameter*2;j=j+1){
            vec_g[j]=0.0;
            for(i=0;i<2*N_y*N_x;i=i+1){
                vec_g[j]+=(jacobian[i][j]*vec_e[i]);
                }
            }
        for(i=0;i<parameter*2;i=i+1){
            //printf("%d %d %.3lf \n",num,i,vec_g[i]);
            }

    ///連立1次方程式
        for(i=0;i<parameter*2;i=i+1){
            bly=-vec_g[i];
            for(j=0;j<i;j=j+1){
                bly-=mat_L[i][j]*delta_a[j];
                }
            delta_a[i]=bly/mat_L[i][i];
            }
        for(i=0;i<parameter*2;i=i+1){
            //printf("%.3lf %.3lf \n",vec_id[i],delta_a[i]);
            delta_a[i]=vec_id[i]*delta_a[i];
            }
    for(i=0;i<parameter*2;i=i+1){
        printf("%d \t %.3lf \n",i,delta_a[i]);
        }

    ///vec_aの更新
    for(i=0;i<parameter*2;i=i+1){
        vec_a[i]+=delta_a[i];
        }

    num++;
    }while(num<10);

///2-11.datファイルに出力
    sprintf(write_file01,"%s//%s",output_folder01,output_data_file1);///
    outfile_1=fopen(write_file01,"w");
    for(i=0;i<number_y-3;i=i+1){
        for(j=0;j<number_x-3;j=j+1){
            fprintf(outfile_1,"%d \t %d \t %.3lf \t %.3lf \n",(j+2)*grid_size_x,(i+2)*grid_size_y,vector_x[i][j],vector_y[i][j]);
            }
        fprintf(outfile_1,"\n");
        }
    fclose(outfile_1);

    sprintf(write_file02,"%s//%s",output_folder02,output_data_file2);///
    outfile_2=fopen(write_file02,"w");
    for(i=0;i<2*parameter;i=i+1){
        fprintf(outfile_2,"%.10lf \n",vec_a[i]);
        }
    fclose(outfile_2);

    return(0);
}
