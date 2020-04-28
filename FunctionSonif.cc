//------------------------------------------------------
// module  : FunctionSonif.cc
// auteur  : Mignotte Max, Toffa Ohini
// date    :
// version : 1.0
// langage : C++
// labo    : DIRO
// note    :
//------------------------------------------------------
// 

//------------------------------------------------
// FICHIERS INCLUS -------------------------------
//------------------------------------------------
#include "FunctionSonif.h"

//>Tabular[Saturation]-(C1-C8)
const int TabSat[8][7]={ {2,3,4,5,6,7,8},
                         {1,3,0,4,5,6,7},
                         {2,4,1,5,0,6,7},
                         {5,3,6,2,7,1,8},
                         {4,6,3,7,2,8,1},
                         {5,7,4,3,2,1,8},
                         {6,5,4,3,2,1,8},
                         {7,6,5,4,3,2,1}, };

const int DO=262;
const int RE=294;
const int MI=330;
const int FA=349;
const int SOL=392;
const int LA=440;
const int SI=494;
const int TabFreqNote[7]={ DO,RE,MI,FA,SOL,LA,SI };
const float MAX_AMP=1.0;
const float MIN_AMP=-1.0;

const int F_ECH=16384;
const int NB_ECH=F_ECH;

//>Model Parameters
const float MAGDIS=5.0; //0.0
const float PHADIS=5.0;  //0.0
const float IMPGRDLOC=0.90;

//--------------------------//
//-- Matrice de Flottant ---//
//--------------------------//
//---------------------------------------------------------
//  alloue de la memoire pour une matrice 1d de float
//----------------------------------------------------------
float* fmatrix_allocate_1d(int hsize)
{
    float* matrix;
    matrix=new float[hsize]; return matrix; }

//----------------------------------------------------------
//  alloue de la memoire pour une matrice 2d de float
//----------------------------------------------------------
float** fmatrix_allocate_2d(int vsize,int hsize)
{
    float** matrix;
    float *imptr;

    matrix=new float*[vsize];
    imptr=new  float[(hsize)*(vsize)];
    for(int i=0;i<vsize;i++,imptr+=hsize) matrix[i]=imptr;
    return matrix;
}

//----------------------------------------------------------
// alloue de la memoire pour une matrice 3d de int
//----------------------------------------------------------
float*** fmatrix_allocate_3d(int dsize,int vsize,int hsize)
{
    float*** matrix;
    matrix=new float**[dsize];

    for(int i=0;i<dsize;i++)
        matrix[i]=fmatrix_allocate_2d(vsize,hsize);
    return matrix;
}

//----------------------------------------------------------
// libere la memoire de la matrice 1d de float
//----------------------------------------------------------
void free_fmatrix_1d(float* pmat)
{ delete[] pmat; }

//----------------------------------------------------------
// libere la memoire de la matrice 2d de float
//----------------------------------------------------------
void free_fmatrix_2d(float** pmat)
{ delete[] (pmat[0]);
    delete[] pmat;}

//----------------------------------------------------------
// libere la memoire de la matrice 3d de float
//----------------------------------------------------------
void free_fmatrix_3d(float*** pmat,int dsize)
{ for(int i=0;i<dsize;i++)
    {
        delete[] (pmat[i][0]);
        delete[] (pmat[i]);
    }
    delete[] (pmat); }

//--------------------//
//-- Matrix Gestion --//
//--------------------//
//----------------------------------------------------------
// copie une matrice dans une autre
//----------------------------------------------------------
void CopyMat(float*** pmatl1,float*** pmatl2,int lgth,int wdth)
{
    for(int k=0;k<3;k++) for(int i=0;i<lgth;i++) for(int j=0;j<wdth;j++)
    { pmatl2[k][i][j]=pmatl1[k][i][j]; }
}

//----------------------------------------------------------
// copie une matrice dans une autre
//----------------------------------------------------------
void CopyMat(float** pmatl1,float** pmatl2,int lgth,int wdth)
{
    for(int i=0;i<lgth;i++) for(int j=0;j<wdth;j++)
    { pmatl2[i][j]=pmatl1[i][j]; }
}

//----------------------------------------------------------
// Recal                                                    
//----------------------------------------------------------
void Recal(float** mat,int lgth,int wdth)
{
    int i,j;
    float max,min;

    /*Initialisation*/
    min=mat[0][0];

    /*Recherche du min*/
    for(i=0;i<lgth;i++) for(j=0;j<wdth;j++)
        if (mat[i][j]<min) min=mat[i][j];

    /*plus min*/
    for(i=0;i<lgth;i++) for(j=0;j<wdth;j++)
        mat[i][j]-=min;

    max=mat[0][0];
    /*Recherche du max*/
    for(i=0;i<lgth;i++) for(j=0;j<wdth;j++)
        if (mat[i][j]>max) max=mat[i][j];

    /*Recalibre la matrice*/
    for(i=0;i<lgth;i++) for(j=0;j<wdth;j++)
        mat[i][j]*=(GREY_LEVEL/max);
}

//----------------------------------------------------------
// Recal                                                    
//----------------------------------------------------------
void Recal(float*** mat,int lgth,int wdth)
{
    Recal(mat[0],lgth,wdth);
    Recal(mat[1],lgth,wdth);
    Recal(mat[2],lgth,wdth);
}

//----------------------------------------------------------
// ImgLog                                                    
//----------------------------------------------------------
void LogImg(float** mat,int lgth,int wdth)
{
    int i,j;
    for(i=0;i<lgth;i++) for(j=0;j<wdth;j++) mat[i][j]=log(1+mat[i][j]);
}

//----------//
//-- FILE --//
//----------//
//----------------------------------------------------------
// Get Length and Width
//----------------------------------------------------------
void GetLengthWidth(char* path,int* length,int* width)
{
    unsigned char var;
    int   temp;
    char* tempc;
    char stringTmp1[100];
    char stringTmp2[100];
    int ta1,ta2;
    FILE *fic;

    //ouverture du fichier
    fic=fopen(path,"r");
    if (fic==NULL)
    { printf("\n- Grave erreur a l'ouverture de %s -\n",path);
        exit(-1); }

    //recuperation de l'entete
    tempc=fgets(stringTmp1,100,fic);
    for(;;) { temp=fread(&var,1,1,fic); if (var==35) tempc=fgets(stringTmp2,100,fic);
        else break; }
    fseek(fic,-1,SEEK_CUR);
    temp=fscanf(fic,"%d %d",&ta1,&ta2);
    if (0) printf("[%d][%s]",temp,tempc);

    //enregistrement
    (*length)=ta2;
    (*width)=ta1;

    //fermeture du fichier
    fclose(fic);
}

//----------------------------------------------------------
// load pgm image 
//----------------------------------------------------------
void LoadImagePgm(char* path,float** data,int length,int width)
{
    int i,j;
    int   temp;
    char* tempc;
    unsigned char var;
    char header[100];
    char* ptr;
    int ta1,ta2,ta3;
    FILE *fic;

    //Open file
    fic=fopen(path,"r");
    if (fic==NULL)
    { printf("\n -> Grave erreur a l'ouverture de %s !\n",path);
        exit(-1); }

    tempc=fgets(header,100,fic);
    if ( (header[0]!=80) ||    /* 'P' */
         (header[1]!=53) ) {   /* '5' */
        fprintf(stderr,"Image %s is not PGM.\n",path);
        exit(1); }

    tempc=fgets(header,100,fic);
    while(header[0]=='#') tempc=fgets(header,100,fic);

    ta1=strtol(header,&ptr,0);
    ta2=atoi(ptr);
    tempc=fgets(header,100,fic);
    ta3=strtol(header,&ptr,0);
    
    //Load
    for(i=0;i<length;i++) for(j=0;j<width;j++)
    { temp=fread(&var,1,1,fic);
        data[i][j]=var; }

    //Close file
    if (0) printf("[%d][%s][%d:%d:%d]",temp,tempc,ta1,ta2,ta3);
    fclose(fic);
}


//----------------------------------------------------------
// load ppm image 
//----------------------------------------------------------
void LoadImagePpm(char* path,float*** data,int length,int width)
{
    int i,j;
    int   temp;
    char* tempc;
    unsigned char varr,varb,varv;
    char header[100];
    char* ptr;
    int ta1,ta2,ta3;
    FILE *fic;

    //Open file
    fic=fopen(path,"r");
    if (fic==NULL)
    { printf("\n -> Grave erreur a l'ouverture de %s !\n",path);
        exit(-1); }

    tempc=fgets(header,100,fic);
    if ( (header[0]!=80) ||    /* 'P' */
         (header[1]!=54) ) {   /* '6' */
        fprintf(stderr,"Image %s is not PPM.\n",path);
        exit(1); }

    tempc=fgets(header,100,fic);
    while(header[0]=='#') tempc=fgets(header,100,fic);

    ta1=strtol(header,&ptr,0);
    ta2=atoi(ptr);
    tempc=fgets(header,100,fic);
    ta3=strtol(header,&ptr,0);
    
    //Load
    for(i=0;i<length;i++) for(j=0;j<(width*3);j+=3)
    {
        temp=fread(&varr,1,1,fic);
        temp=fread(&varv,1,1,fic);
        temp=fread(&varb,1,1,fic);
        data[0][i][j/3]=varr;
        data[1][i][j/3]=varv;
        data[2][i][j/3]=varb;
    }

    //Close file
    if (0) printf("[%d][%s][%d:%d:%d]",temp,tempc,ta1,ta2,ta3);
    fclose(fic);
}

//----------------------------------------------------------
//  LoadPrintImagePpm  
//----------------------------------------------------------
void LoadPrintImagePpm(char* path,int lgth,int wdth)
{
    int i,j;
    int temp;
    char* tempc;
    unsigned char varb,varv,varr;
    char stringTmp1[200],stringTmp2[200];
    int ta1,ta2,ta3;
    FILE *fic;

    //Allocation
    float*** img=fmatrix_allocate_3d(3,lgth,wdth);

    //enregistrement du chemin;
    if (strstr(path,".ppm")==NULL) strcat(path,".ppm");

    //ouverture du fichier
    fic=fopen(path,"r");
    if (fic==NULL)
    { printf("\n ->> Grave erreur a l'ouverture de %s -\n",path);
        exit(-1); }

    //recuperation de l'entete
    tempc=fgets(stringTmp1,100,fic); //P6
    tempc=fgets(stringTmp2,100,fic); //Comment

    for(;;)
    { if (((int)stringTmp2[0]==10)||((int)stringTmp2[0]==35))
            tempc=fgets(stringTmp2,100,fic);
        else break; }

    temp=sscanf(tempc," %d %d",&ta1,&ta2);
    temp=fscanf(fic," %d\n",&ta3);

    //printf("\n > Entete: %s  %d %d \n %d\n",stringTmp1,ta1,ta2,ta3);
    if (0) printf("[%d][%c]",temp,*tempc);

    //chargement dans la matrice
    for(i=0;i<lgth;i++) for(j=0;j<(wdth*3);j+=3)
    { temp=fread(&varr,1,1,fic);
        temp=fread(&varv,1,1,fic);
        temp=fread(&varb,1,1,fic);
        img[0][i][j/3]=varr;
        img[1][i][j/3]=varv;
        img[2][i][j/3]=varb; }

    //Print
    for(i=0;i<lgth;i++)
    { varr=img[0][i][wdth/2];
        varv=img[1][i][wdth/2];
        varb=img[2][i][wdth/2];
        printf("\n {%d,%d,%d},",varr,varv,varb); }

    //liberationMemoire
    if (img) free_fmatrix_3d(img,3);

    //fermeture du fichier
    fclose(fic);
}

//----------------------------------------------------------
// Sauvegarde de l'image de nom <name> au format pgm                        
//----------------------------------------------------------                
void SaveImagePgm(char* name,float** mat,int lgth,int wdth)
{
    int i,j;
    char buff[300];
    FILE* fic;

    //--extension--
    strcpy(buff,name);
    strcat(buff,".pgm");

    //--ouverture fichier--
    fic=fopen(buff,"wb");
    if (fic==NULL)
    { printf("Probleme dans la sauvegarde de %s",buff);
        exit(-1); }
    printf("\n Sauvegarde de %s au format pgm\n",buff);

    //--sauvegarde de l'entete--
    fprintf(fic,"P5");
    fprintf(fic,"\n# IMG Module");
    fprintf(fic,"\n%d %d",wdth,lgth);
    fprintf(fic,"\n255\n");

    //--enregistrement--
    for(i=0;i<lgth;i++) for(j=0;j<wdth;j++)
        fprintf(fic,"%c",(char)mat[i][j]);

    //--fermeture fichier--
    fclose(fic);
}

//----------------------------------------------------------
// Sauvegarde de l'image de nom <name> au format ppm        
//----------------------------------------------------------
void SaveImagePpm(char* Name,float*** matrvb,int lgth,int wdth)
{
    int i,j;
    char buff[200];
    FILE* fuser;

    //extension
    strcpy(buff,Name);
    if (strstr(buff,".ppm")==NULL) strcat(buff,".ppm");

    //ouverture fichier
    fuser=fopen(buff,"w");
    if (fuser==NULL)
    { printf(" probleme dans la sauvegarde de %s",buff);
        exit(-1); }

    //affichage
    printf("\n sauvegarde dans -> %s au format %s [%d][%d]",buff,".ppm",wdth,lgth);
    fflush(stdout);

    //sauvegarde de l'entete
    fprintf(fuser,"P6");
    fprintf(fuser,"\n# IMG Module");
    fprintf(fuser,"\n%d %d",wdth,lgth);
    fprintf(fuser,"\n255\n");

    //enregistrement
    for(i=0;i<lgth;i++) for(j=0;j<wdth;j++)
    {
        fprintf(fuser,"%c",(char)matrvb[0][i][j]);
        fprintf(fuser,"%c",(char)matrvb[1][i][j]);
        fprintf(fuser,"%c",(char)matrvb[2][i][j]);
    }

    //fermeture fichier
    fclose(fuser);
}


/*--------------*/
/* FOURIER -----*/
/*--------------*/
/*------------------------------------------------*/
/*  FOURN ----------------------------------------*/
/*------------------------------------------------*/
void fourn(float data[], unsigned long nn[], int ndim, int isign)
{
    int idim;
    unsigned long i1,i2,i3,i2rev,i3rev,ip1,ip2,ip3,ifp1,ifp2;
    unsigned long ibit,k1,k2,n,nprev,nrem,ntot;
    float tempi,tempr;
    double theta,wi,wpi,wpr,wr,wtemp;

    for (ntot=1,idim=1;idim<=ndim;idim++)
        ntot *= nn[idim];
    nprev=1;
    for (idim=ndim;idim>=1;idim--) {
        n=nn[idim];
        nrem=ntot/(n*nprev);
        ip1=nprev << 1;
        ip2=ip1*n;
        ip3=ip2*nrem;
        i2rev=1;
        for (i2=1;i2<=ip2;i2+=ip1) {
            if (i2 < i2rev) {
                for (i1=i2;i1<=i2+ip1-2;i1+=2) {
                    for (i3=i1;i3<=ip3;i3+=ip2) {
                        i3rev=i2rev+i3-i2;
                        SWAP(data[i3],data[i3rev]);
                        SWAP(data[i3+1],data[i3rev+1]);
                    }
                }
            }
            ibit=ip2 >> 1;
            while (ibit >= ip1 && i2rev > ibit) {
                i2rev -= ibit;
                ibit >>= 1;
            }
            i2rev += ibit;
        }
        ifp1=ip1;
        while (ifp1 < ip2) {
            ifp2=ifp1 << 1;
            theta=isign*6.28318530717959/(ifp2/ip1);
            wtemp=sin(0.5*theta);
            wpr = -2.0*wtemp*wtemp;
            wpi=sin(theta);
            wr=1.0;
            wi=0.0;
            for (i3=1;i3<=ifp1;i3+=ip1) {
                for (i1=i3;i1<=i3+ip1-2;i1+=2) {
                    for (i2=i1;i2<=ip3;i2+=ifp2) {
                        k1=i2;
                        k2=k1+ifp1;
                        tempr=(float)wr*data[k2]-(float)wi*data[k2+1];
                        tempi=(float)wr*data[k2+1]+(float)wi*data[k2];
                        data[k2]=data[k1]-tempr;
                        data[k2+1]=data[k1+1]-tempi;
                        data[k1] += tempr;
                        data[k1+1] += tempi;
                    }
                }
                wr=(wtemp=wr)*wpr-wi*wpi+wr;
                wi=wi*wpr+wtemp*wpi+wi;
            }
            ifp1=ifp2;
        }
        nprev *= n;
    }
}


/*----------------------------------------------------------*/
/* FFTDD                                                    */
/*----------------------------------------------------------*/
void FFTDD(float** mtxR,float** mtxI,int lgth, int wdth)
{
    int i,j;
    int posx,posy;

    float* data;
    float* ImgFreqR;
    float* ImgFreqI;
    unsigned long* nn;

    /*allocation memoire*/
    data=(float*)malloc(sizeof(float)*(2*wdth*lgth)+1);
    ImgFreqR=(float*)malloc(sizeof(float)*(wdth*lgth));
    ImgFreqI=(float*)malloc(sizeof(float)*(wdth*lgth));
    nn=(unsigned long*)malloc(sizeof(unsigned long)*(FFT2D+1));

    /*Remplissage de nn*/
    nn[1]=lgth; nn[2]=wdth;

    /*Remplissage de data*/
    for(i=0;i<lgth;i++) for(j=0;j<wdth;j++)
    { data[2*(i*lgth+j)+1]=mtxR[i][j];
        data[2*(i*lgth+j)+2]=mtxI[i][j]; }

    /*FFTDD*/
    fourn(data,nn,FFT2D,FFT);

    /*Remplissage*/
    for(i=0;i<(wdth*lgth);i++)
    { ImgFreqR[i]=data[(2*i)+1];
        ImgFreqI[i]=data[(2*i)+2];  }

    /*Conversion en Matrice*/
    for(i=0;i<(wdth*lgth);i++)
    { posy=(int)(i/wdth);
        posx=(int)(i%wdth);

        mtxR[posy][posx]=ImgFreqR[i];
        mtxI[posy][posx]=ImgFreqI[i];}

    /*Liberation memoire*/
    free(data);
    free(ImgFreqR);
    free(ImgFreqI);
    free(nn);
}


/*----------------------------------------------------------*/
/* IFFTDD                                                   */
/*----------------------------------------------------------*/
void IFFTDD(float** mtxR,float**  mtxI,int lgth,int wdth)
{
    int i,j;
    int posx,posy;

    float* data;
    float* ImgFreqR;
    float* ImgFreqI;
    unsigned long* nn;

    /*allocation memoire*/
    data=(float*)malloc(sizeof(float)*(2*wdth*lgth)+1);
    ImgFreqR=(float*)malloc(sizeof(float)*(wdth*lgth));
    ImgFreqI=(float*)malloc(sizeof(float)*(wdth*lgth));
    nn=(unsigned long*)malloc(sizeof(unsigned long)*(FFT2D+1));

    /*Recadrege*/

    /*Remplissage de nn*/
    nn[1]=lgth; nn[2]=wdth;

    /*Remplissage de data*/
    for(i=0;i<lgth;i++) for(j=0;j<wdth;j++)
    { data[2*(i*lgth+j)+1]=mtxR[i][j];
        data[2*(i*lgth+j)+2]=mtxI[i][j]; }

    /*FFTDD*/
    fourn(data,nn,FFT2D,IFFT);

    /*Remplissage*/
    for(i=0;i<(wdth*lgth);i++)
    { ImgFreqR[i]=data[(2*i)+1];
        ImgFreqI[i]=data[(2*i)+2]; }

    /*Conversion en Matrice*/
    for(i=0;i<(wdth*lgth);i++)
    { posy=(int)(i/wdth);
        posx=(int)(i%wdth);

        mtxR[posy][posx]=ImgFreqR[i]/(wdth*lgth);
        mtxI[posy][posx]=ImgFreqI[i]/(wdth*lgth); }

    /*Liberation memoire*/
    free(data);
    free(ImgFreqR);
    free(ImgFreqI);
    free(nn);
}

/*----------------------------------------------------------*/
/* FFT1D                                                    */
/*----------------------------------------------------------*/
void FFT1D(float* vctR,float* vctI,int lgth)
{
    int i;
    float* data;
    unsigned long* nn;

    /*allocation memoire*/
    data=(float*)malloc(sizeof(float)*((2*lgth)+1));
    nn=(unsigned long*)malloc(sizeof(unsigned long)*2);
    nn[1]=lgth;

    /*Remplissage de data*/
    for(i=0;i<lgth;i++)
    { data[(2*i)+1]=vctR[i];
        data[(2*i)+2]=vctI[i];
    }

    /*FFTDD*/
    fourn(data,nn,1,FFT);

    /*Resultat*/
    for(i=0;i<lgth;i++)
    { vctR[i]=data[(2*i)+1];
        vctI[i]=data[(2*i)+2];  }

    /*Liberation memoire*/
    free(data);
    free(nn);
}

/*----------------------------------------------------------*/
/* IFFT1D                                                   */
/*----------------------------------------------------------*/
void IFFT1D(float* vctR,float* vctI,int lgth)
{
    int i;
    float* data;
    unsigned long* nn;

    /*allocation memoire*/
    data=(float*)malloc(sizeof(float)*((2*lgth)+1));
    nn=(unsigned long*)malloc(sizeof(unsigned long)*2);
    nn[1]=lgth;

    /*Remplissage de data*/
    for(i=0;i<lgth;i++)
    { data[(2*i)+1]=vctR[i];
        data[(2*i)+2]=vctI[i]; }

    /*FFTDD*/
    fourn(data,nn,1,IFFT);

    /*Resultat*/
    for(i=0;i<lgth;i++)
    { vctR[i]=data[(2*i)+1];
        vctI[i]=data[(2*i)+2];  }

    /*Liberation memoire*/
    free(data);
    free(nn);
}

/*----------------------------------------------------------*/
/* Mod2                                                     */
/*----------------------------------------------------------*/
void Mod(float** matM,float** matR,float** matI,int lgth,int wdth)
{
    int i,j;

    /*Calcul du module*/
    for(i=0;i<lgth;i++) for(j=0;j<wdth;j++)
        matM[i][j]=sqrt((matR[i][j]*matR[i][j])+(matI[i][j]*matI[i][j]));
}

/*----------------------------------------------------------*/
/* Mod2                                                     */
/*----------------------------------------------------------*/
void Mod(float* VctM,float* VctR,float* VctI,int lgth)
{
    int i;

    /*Calcul du module*/
    for(i=0;i<lgth;i++) VctM[i]=sqrt((VctR[i]*VctR[i])+(VctI[i]*VctI[i]));
}

/*----------------------------------------------------------*/
/* Mult                                                     */
/*----------------------------------------------------------*/
void Mult(float** mat,float coef,int lgth,int wdth)
{
    int i,j;

    /*multiplication*/
    for(i=0;i<lgth;i++) for(j=0;j<wdth;j++)
    { mat[i][j]*=coef;
        if (mat[i][j]>GREY_LEVEL) mat[i][j]=GREY_LEVEL; }
}

/*----------------------------------------------------------*/
/* Mult 2 matrices complexes                                */
/*----------------------------------------------------------*/
void MultMatrix(float** matRout,float** matIout,float** mat1Rin,float** mat1Iin,
                float** mat2Rin,float** mat2Iin,int lgth,int wdth)
{
    int i,j;

    for(i=0;i<lgth;i++) for(j=0;j<wdth;j++)
    { matRout[i][j]=mat1Rin[i][j]*mat2Rin[i][j]-mat1Iin[i][j]*mat2Iin[i][j];
        matIout[i][j]=mat1Rin[i][j]*mat2Iin[i][j]+mat2Rin[i][j]*mat1Iin[i][j]; }
}

/*----------------------------------------------------------*/
/* Matrice complexe au carre                                */
/*----------------------------------------------------------*/
void SquareMatrix(float** matRout,float** matIout,float** matRin,float** matIin,int lgth,int wdth)
{
    int i,j;

    for(i=0;i<lgth;i++) for(j=0;j<wdth;j++)
    { matRout[i][j]=SQUARE(matRin[i][j])-SQUARE(matIin[i][j]);
        matIout[i][j]=2*matRin[i][j]*matIin[i][j]; }
}

/*----------------------------------------------------------*/
/*  CenterImg                                               */
/*   > Utilisant l'algo des  4 quadrants remis              */
/*----------------------------------------------------------*/
void CenterImg(float** mat,int lgth,int wdth)
{
    int i,j;
    int ci,cj;
    float** mattmp;

    /*Initialisation*/
    ci=(int)(lgth/2);
    cj=(int)(wdth/2);

    /*Allocation memoire*/
    mattmp=fmatrix_allocate_2d(lgth,wdth);

    /*Recadrage*/
    for(i=0;i<ci;i++) for(j=0;j<cj;j++)
        mattmp[ci+i][cj+j]=mat[i][j];

    for(i=ci;i<lgth;i++) for(j=cj;j<wdth;j++)
        mattmp[i-ci][j-cj]=mat[i][j];

    for(i=0;i<ci;i++) for(j=cj;j<wdth;j++)
        mattmp[ci+i][j-cj]=mat[i][j];

    for(i=ci;i<lgth;i++) for(j=0;j<cj;j++)
        mattmp[i-ci][cj+j]=mat[i][j];

    /*Transfert*/
    for(i=0;i<lgth;i++) for(j=0;j<wdth;j++)
        mat[i][j]=mattmp[i][j];

    /*desallocation memoire*/
    free_fmatrix_2d(mattmp);
}

/*----------------------------------------------------------*/
/*  CenterVct                                               */
/*----------------------------------------------------------*/
void CenterVct(float* vct,int lgth)
{
    int i;
    int ci;
    float* vcttmp;

    /*Initialisation*/
    ci=(int)(lgth/2);

    /*Allocation memoire*/
    vcttmp=fmatrix_allocate_1d(lgth);

    /*Recadrage*/
    for(i=0;i<ci;i++) vcttmp[ci+i]=vct[i];
    for(i=ci;i<lgth;i++) vcttmp[i-ci]=vct[i];

    /*Transfert*/
    for(i=0;i<lgth;i++) vct[i]=vcttmp[i];

    /*desallocation memoire*/
    free_fmatrix_1d(vcttmp);
}

//-----------------------------//
// -LECTURE/SAUVEGARDE SIGNAL -//
//-----------------------------//
/*----------------------------------------------------------*/
/* Chargement du signal de nom <name>.dat                   */
/*----------------------------------------------------------*/
float* LoadSignalDat(char* name,int* length)
{
    int i;
    float ech1,ech2;
    char buff[NBCHAR];
    float Tech;
    FILE *fic;

    //Allocation
    float** Vct2D=fmatrix_allocate_2d(2,10000000);

    //nom du fichier dat
    strcpy(buff,name);
    strcat(buff,".dat");
    printf("\n  >> Ouverture de %s",buff);

    //ouverture du fichier
    fic=fopen(buff,"r");
    if (fic==NULL)
    { printf("\n- Grave erreur a l'ouverture de %s  -\n",buff);
        exit(-1); }

    //Lecture Donnée & Longueur & Periode Ech
    for(i=0;;i++)
    { fscanf(fic,"%f %f",&ech1,&ech2);
        if (feof(fic)) break;
        //printf("\n[%f::%f]",ech1,ech2);
        Vct2D[0][i]=ech1;
        Vct2D[1][i]=ech2; }

    (*length)=i;
    Tech=Vct2D[0][1];
    Tech=1.0/Tech;
    Tech=(int)Tech;
    printf(" (%d echantillons)",(*length));
    printf("\n  >> Techantillonnage:: %.0f echantillons/seconde",Tech);

    //Chargement
    float* VctFinal=fmatrix_allocate_1d((*length));
    for(i=0;i<(*length);i++) VctFinal[i]=Vct2D[1][i];

    //End
    fclose(fic);
    (*length)=i;
    free_fmatrix_2d(Vct2D);
    return VctFinal;
}

/*----------------------------------------------------------*/
/* Sauvegarde du signal de nom <name>.dat                   */
/*----------------------------------------------------------*/
void SaveSignalDat(char* name,float* vct,int length)
{
    int i;
    char buff[NBCHAR];
    FILE* fic;

    /*--extension--*/
    strcpy(buff,name);
    strcat(buff,".dat");

    /*--ouverture fichier--*/
    fic=fopen(buff,"w");
    if (fic==NULL)
    { printf(" Probleme dans la sauvegarde de %s",buff);
        exit(-1); }
    printf("\n Sauvegarde de %s au format dat\n",name);

    /*--enregistrement--*/
    for(i=0;i<length;i++) fprintf(fic,"%f %f\n",(float)i,vct[i]);

    /*--fermeture fichier--*/
    fclose(fic);
}

/*----------------------------------------------------------*/
/* Sauvegarde du signal de nom <name>.dat                   */
/*----------------------------------------------------------*/
void SaveSignalDatWav(char* name,float* vct,int length,int SampRate)
{
    int i;
    char buffName[NBCHAR];
    char buffSyst[NBCHAR];
    FILE* fic;
    const int DEBUG=0;

    /*--extension--*/
    strcpy(buffName,name);

    /*--ouverture fichier--*/
    fic=fopen("tmp.dat","w");
    if (fic==NULL)
    { printf(" Probleme dans la sauvegarde de %s",buffName);
        exit(-1); }
    if (DEBUG) printf("\n Sauvegarde de %s au format wav",name);

    /*--enregistrement--*/
    fprintf(fic,"; Sample Rate %d",SampRate);
    fprintf(fic,"; Channels 1\n");
    for(i=0;i<length;i++) fprintf(fic,"%f %f\n",(float)i,vct[i]);

    /*--fermeture fichier--*/
    fclose(fic);

    /*--Conversion Dat -> Wav--*/
    //sprintf(buffSyst,"sox tmp.dat %s.wav&",buffName);
    //sprintf(buffSyst,"sox tmp.dat tmp.wav && mv tmp.wav %s.wav",buffName);
    sprintf(buffSyst,"sox tmp.dat tmp.wav && mv tmp.wav %s.wav&",buffName);
    //sprintf(buffSyst,"sox tmp.dat %s.wav",buffName);
    system(buffSyst);
    if (DEBUG) printf("\n Appel Systeme::sox tmp.dat %s.wav\n\n",buffName);
}

//--------------------
//-- SPACE COLOR -----
//--------------------
//----------------------------------------------------------           
//   RGB ->  [R,G,B] = [0,255]
//   HSV ->  H=[0::360]   S=[0::1]  V=[0::255]  
//   HSL ->  H=[0::3.142] S=[0::1]  L=[0::1]    
//   Recalage sur [0-255] avec :: (C-min)*255/(max-min)
//----------------------------------------------------------    
void ImgRGBToImgHsl(float*** ImgRGB,float*** ImgHSL,int Lgth,int Wdth)
{
    int i,j;
    float col1,col2,col3;

    //>HSV
    if (0)
        for(i=0;i<Lgth;i++) for(j=0;j<Wdth;j++)
        { RGBToHSV(ImgRGB[0][i][j],ImgRGB[1][i][j],ImgRGB[2][i][j],&col1,&col2,&col3);
            ImgHSL[0][i][j]=(col1/360.0)*255.0;
            ImgHSL[1][i][j]=col2*255.0;
            ImgHSL[2][i][j]=col3; }

    //>HSL
    if (0)
        for(i=0;i<Lgth;i++) for(j=0;j<Wdth;j++)
        { RGBToHSL(ImgRGB[0][i][j],ImgRGB[1][i][j],ImgRGB[2][i][j],&col1,&col2,&col3);
            float h = 255*col1/3.142;
            ImgHSL[0][i][j]=h;
            ImgHSL[1][i][j]=col2*254;
            ImgHSL[2][i][j]=col3*255; }

    //>HSL2
    if (1)
        for(i=0;i<Lgth;i++) for(j=0;j<Wdth;j++)
        { RGBToHSL2(ImgRGB[0][i][j],ImgRGB[1][i][j],ImgRGB[2][i][j],&col1,&col2,&col3);
            float h = 255*(col1/360.0);
            ImgHSL[0][i][j]=h;
            ImgHSL[1][i][j]=col2*254;
            ImgHSL[2][i][j]=col3*255; }
}

//----------------------------------------------------------*/
// RGBToHSV 
// --------
// H=[0::360]   S=[0::1]  V=[0::255]
//                        
//----------------------------------------------------------*/
void RGBToHSV(float r,float g,float b,float* h,float* s,float* v)
{
    float maxRGB=MAX(MAX(r,g),MAX(g,b));
    float minRGB=MIN(MIN(r,g),MIN(g,b));

    (*v)=maxRGB;

    if (maxRGB==0)   (*s)=0.0;
    else             (*s)= (maxRGB-minRGB)/maxRGB;


    if ((*s)==0)     (*h)=0.0;
    else
    { if (r==maxRGB)       (*h)=(g-b)/(maxRGB-minRGB);
        else if (g==maxRGB)  (*h)=2+(b-r)/(maxRGB-minRGB);
        else                 (*h)=4+(r-g)/(maxRGB-minRGB);

        (*h)*=60;
        if ((*h)<0.0)        (*h)+=360;  }
}

//----------------------------------------------------------*/
// RGBToHSL
//
//  H=ACOS[ 0.5(r-g)+(r-b)) / sqrt[ (r-g)^2=(r-b)*(g-b) ] ]
// 
//  L=(1.0/2.0)* [MAX(r,g,b)+MIN(r,g,bl)] 
//
// if (l<=0.5) S=MAX(r,g,b)-MIN(r,g,b)/(MAX(r,g,b)+MIN(r,g,b))
// if (l>0.5)  S=MAX(r,g,b)-MIN(r,g,b)/(2-(MAX(r,g,b)+MIN(r,g,b)))
//                   
//----------------------------------------------------------*/
void RGBToHSL(float r,float g,float bl,float* h,float* s,float* l)
{
    float tmp1,tmp2;
    // r=0; g=255; bl=153;
    //RGB en [0-1]
    r/=255.0;
    g/=255.0;
    bl/=255.0;

    //min/maxRGB
    float maxRGB=MAX(MAX(r,g),MAX(g,bl));
    float minRGB=MIN(MIN(r,g),MIN(g,bl));

    //Conversion RGB -> HSV
    tmp1=(0.5)*((r-g)+(r-bl));
    tmp2=sqrt((CARRE(r-g))+(r-bl)*(g-bl));
    if (tmp2) (*h)=acos(tmp1/tmp2);
    else      (*h)=0.0;

    (*l)=(1.0/2.0)*(maxRGB+minRGB);

    (*s)=0.0;
    if ((*l)<=0.5) if(maxRGB+minRGB) (*s)=(maxRGB-minRGB)/(maxRGB+minRGB);
    if ((*l)>0.5)  (*s)=(maxRGB-minRGB)/(2.0-(maxRGB+minRGB));

    if (isnan(*(h))) *h=0.0;
    if (isnan(*(s))) *s=0.0;
    if (isnan(*(l))) *l=0.0;
}

//----------------------------------------------------------*/
// RGBToHSL
//
//  H=[0::360]
//
//  L=(1.0/2.0)* [MAX(r,g,b)+MIN(r,g,bl)]
//
// if (l<=0.5) S=MAX(r,g,b)-MIN(r,g,b)/(MAX(r,g,b)+MIN(r,g,b))
// if (l>0.5)  S=MAX(r,g,b)-MIN(r,g,b)/(2-(MAX(r,g,b)+MIN(r,g,b)))
//
//----------------------------------------------------------*/
void RGBToHSL2(float r,float g,float b,float* h,float* s,float* l)
{
    float tmp1,tmp2;
    //r=0; g=255; b=153; //(104)68 255 0  (156)0 255 153
    //RGB en [0-1]
    r/=255.0;
    g/=255.0;
    b/=255.0;

    //min/maxRGB
    float maxRGB=MAX(MAX(r,g),MAX(g,b));
    float minRGB=MIN(MIN(r,g),MIN(g,b));

    //Conversion RGB -> HSV
    tmp1=(0.5)*((r-g)+(r-b));
    tmp2=sqrt((CARRE(r-g))+(r-b)*(g-b));

    // Luminance
    (*l)=(1.0/2.0)*(maxRGB+minRGB);

    // Saturation
    (*s)=0.0;
    if ((*l)<=0.5) if(maxRGB+minRGB) (*s)=(maxRGB-minRGB)/(maxRGB+minRGB);
    if ((*l)>0.5)  (*s)=(maxRGB-minRGB)/(2.0-(maxRGB+minRGB));

    // Teinte
    if ((*s)==0)     (*h)=0.0;
    else
    { if (r==maxRGB)       (*h)=(g-b)/(maxRGB-minRGB);
        else if (g==maxRGB)  (*h)=2+(b-r)/(maxRGB-minRGB);
        else                 (*h)=4+(r-g)/(maxRGB-minRGB);

        (*h)*=60;
        if ((*h)<0.0)        (*h)+=360;
    }

    if (isnan(*(h))) *h=0.0;
    if (isnan(*(s))) *s=0.0;
    if (isnan(*(l))) *l=0.0;
}

//---------------------------
//- Aleat
//---------------------------
//----------------------------------------------------------
// etourne un nombre aleatoire entre zero et un 
//----------------------------------------------------------
float randomize(void)
{ 
    return ((float)rand()/RAND_MAX);
}

//------------------//
//-- SONIFICATION --//
//------------------//
//---------------------------------------------------
// RecalSignal                                  
//---------------------------------------------------
void RecalSignal(float* sig,int lgth,float minVal,float maxVal)
{
    int k;
    float min,max;
    float scale=maxVal-minVal;
    const int DEBUG=0;

    //>DEBUG
    min=100000;
    max=-100000;
    for(k=0;k<lgth;k++)
    { if (sig[k]<min) min=sig[k];
        if (sig[k]>max) max=sig[k]; }
    if (DEBUG) printf("[%.1f::%.1f]",min,max);
    //fflush(stdout);

    //>Recal [0-scale]>[min-max]
    min=100000;
    max=-100000;
    for(k=0;k<lgth;k++) if (sig[k]<min) min=sig[k];
    for(k=0;k<lgth;k++) sig[k]+=-min;
    for(k=0;k<lgth;k++) if (sig[k]>max) max=sig[k];
    for(k=0;k<lgth;k++) sig[k]*=scale/max;
    for(k=0;k<lgth;k++) sig[k]+=minVal;

    //>DEBUG
    min=100000;
    max=-100000;
    for(k=0;k<lgth;k++)
    { if (sig[k]<min) min=sig[k];
        if (sig[k]>max) max=sig[k]; }
    if (DEBUG) printf(">[%.1f::%.1f]",min,max);
    //fflush(stdout);
}

//---------------------------------------------------
// RecalSignal                                  
//---------------------------------------------------
void RecalSignal_(float* sig,int lgth,float minVal,float maxVal)
{
    int k;
    float moy;
    float prop;

    //>Recal_[Null_Mean]
    moy=0.0;
    for(k=0;k<lgth;k++) moy+=sig[k];
    moy/=lgth;
    for(k=0;k<lgth;k++) sig[k]-=moy;

    //>Recal_[minVal-maxVal]
    for(;;)
    { prop=0;
        for(k=0;k<lgth;k++) if ((sig[k]>maxVal)||(sig[k]<minVal)) prop++;
        prop/=lgth;

        if (prop<0.01) break;
        else { for(k=0;k<lgth;k++) sig[k]*=0.9; }

        //if (1) printf(">");
    }
}

//----------------------------------------------------------
// ValHist1D
//----------------------------------------------------------
void ValHist1D(float** Img,int Lgth,int Wdth,float** ImgSeg,int label,float* VctValHist,int NBBINS)
{
    int i,j,k;
    int nb;
    int ValH;
    const float DIVBIN=256/NBBINS;

    //>Init
    for(k=0;k<NBBINS;k++) VctValHist[k]=0.0;

    //>MainLoop
    nb=0;
    for(i=0;i<Lgth;i++) for(j=0;j<Wdth;j++) if (ImgSeg[i][j]==label)
    { nb++;
        ValH=(int)(Img[i][j]/DIVBIN);
        if (ValH>(NBBINS-1)) ValH=(NBBINS-1);
        VctValHist[ValH]++; }

    for(k=0;k<NBBINS;k++) if (nb) VctValHist[k]/=(nb);
}

//----------------------------------------------------------
// GradHist1D
//----------------------------------------------------------
void GradHist1D(float*** Img,int Lgth,int Wdth,float** ImgSeg,int label,float* VctGradHist,int NBBINS)
{
    int i,j,k;
    int nb;
    int ValH;
    const float DIVBIN=256/NBBINS;

    //>Init
    float** MtGrd=fmatrix_allocate_2d(Lgth,Wdth);
    for(i=0;i<Lgth;i++) for(j=0;j<Wdth;j++)  MtGrd[i][j]=0.0;
    for(i=1;i<(Lgth-1);i++) for(j=1;j<(Wdth-1);j++)
    { for(k=0;k<TROIS;k++)
            MtGrd[i][j]+=fabs(Img[k][i][j]-Img[k][i-1][j])+fabs(Img[k][i][j]-Img[k][i][j-1])+fabs(Img[k][i][j]-Img[k][i-1][j-1]);
        MtGrd[i][j]/=3.0; };

    //>MainLoop
    nb=0;
    for(i=0;i<Lgth;i++) for(j=0;j<Wdth;j++) if (ImgSeg[i][j]==label)
    { nb++;
        ValH=(int)(MtGrd[i][j]/DIVBIN);
        if (ValH>(NBBINS-1)) ValH=(NBBINS-1);
        VctGradHist[ValH]++; }

    for(k=0;k<NBBINS;k++) if (nb) VctGradHist[k]/=(nb);

    //FreeMemory
    if (MtGrd) free_fmatrix_2d(MtGrd);
}

//----------------------------------------------------------
// ComputeHueMixture [float** Img::Hue]
//----------------------------------------------------------
void ComputeHueMixture(float** Img,float** ImgSeg,int Lgth,int Wdth,float** HueMixture,int NBREGMAX,int NBBINS)
{
    int k,l;
    //MainLoop
    for(k=0;k<NBREGMAX;k++) for(l=0;l<NBBINS;l++) HueMixture[k][l]=0.0;
    for(k=0;k<NBREGMAX;k++) ValHist1D(Img,Lgth,Wdth,ImgSeg,k,HueMixture[k],NBBINS);
}

//----------------------------------------------------------
// ComputeGradMixture [float*** Img::RGB Img]
//----------------------------------------------------------
void ComputeGrdMixture(float*** Img,float** ImgSeg,int Lgth,int Wdth,float** GradMixture,int NBREGMAX,int NBBINS)
{
    int k,l;
    //MainLoop
    for(k=0;k<NBREGMAX;k++) for(l=0;l<NBBINS;l++) GradMixture[k][l]=0.0;
    for(k=0;k<NBREGMAX;k++) GradHist1D(Img,Lgth,Wdth,ImgSeg,k,GradMixture[k],NBBINS);
}

//----------------------------------------------------------
// ComputeSatLuMixture [float*** Img::HSL Img]
//----------------------------------------------------------
void ComputeSatLuMixture(float*** Img,float** ImgSeg,int Lgth,int Wdth,float** SatLuMixture,int NBREGMAX)
{
    int i,j,k,l;
    int nb;
    float ValSat,ValLum;

    //MainLoop
    for(k=0;k<NBREGMAX;k++) for(l=0;l<3;l++) SatLuMixture[k][l]=0.0;

    for(k=0;k<NBREGMAX;k++)
    {
        nb=0;
        ValSat=ValLum=0.0;
        for(i=0;i<Lgth;i++) for(j=0;j<Wdth;j++) if (ImgSeg[i][j]==k)
        { nb++;
            ValSat+=Img[1][i][j];
            ValLum+=Img[2][i][j]; }
        if (nb)
        { ValSat/=nb;
            ValLum/=nb; }
        SatLuMixture[k][0]=ValSat;
        SatLuMixture[k][1]=ValLum;
        SatLuMixture[k][2]=nb;
    }
}

//----------------------------------------------------------
//  ComputeGrdMoy[float*** Img::RGB Img]
//----------------------------------------------------------
void ComputeGrdMoy(float*** Img,float** ImgSeg,int Lgth,int Wdth,float* VctGrdMoy,int NBREGMAX)
{
    int i,j,k;
    int nb;
    float GrdMoy;

    //>Init
    float** MtGrd=fmatrix_allocate_2d(Lgth,Wdth);
    for(i=0;i<Lgth;i++) for(j=0;j<Wdth;j++)  MtGrd[i][j]=0.0;
    for(i=1;i<(Lgth-1);i++) for(j=1;j<(Wdth-1);j++)
    { for(k=0;k<TROIS;k++)
            MtGrd[i][j]+=fabs(Img[k][i][j]-Img[k][i-1][j])+fabs(Img[k][i][j]-Img[k][i][j-1])+fabs(Img[k][i][j]-Img[k][i-1][j-1]);
        MtGrd[i][j]/=9.0; }

    //>Estimate
    for(k=0;k<NBREGMAX;k++)
    { nb=0;
        GrdMoy=0.0;
        for(i=0;i<Lgth;i++) for(j=0;j<Wdth;j++) if (ImgSeg[i][j]==k) { nb++; GrdMoy+=MtGrd[i][j]; }
        if (nb) GrdMoy/=nb;
        VctGrdMoy[k]=GrdMoy; }

    //>FreeMemory
    if (MtGrd) free_fmatrix_2d(MtGrd);
}

//----------------------------------------------------------
//  Compute OrientGradVar[float*** Img::RGB Img]
//----------------------------------------------------------
void ComputeOrientGradVar(float*** Img,float** ImgSeg,int Lgth,int Wdth,float* VctOrientGrdVar,int NBREGMAX)
{
    int i,j,k,l;
    int nb;
    float moy,var;
    int orient = 4;

    //>Init
    float*** MtGrd=fmatrix_allocate_3d(orient,Lgth,Wdth);
    float**  VctHist=fmatrix_allocate_2d(NBREGMAX,orient);
    for(k=0;k<orient;k++) for(i=0;i<Lgth;i++) for(j=0;j<Wdth;j++) MtGrd[k][i][j]=0.0;
    for(k=0;k<NBREGMAX;k++) for(l=0;l<orient;l++) VctHist[k][l]=0.0;

    for(i=1;i<(Lgth-1);i++) for(j=1;j<(Wdth-1);j++) for(k=0;k<TROIS;k++) MtGrd[0][i][j]+=fabs(Img[k][i][j]-Img[k][i-1][j]);
    for(i=1;i<(Lgth-1);i++) for(j=1;j<(Wdth-1);j++) for(k=0;k<TROIS;k++) MtGrd[1][i][j]+=fabs(Img[k][i][j]-Img[k][i][j-1]);
    for(i=1;i<(Lgth-1);i++) for(j=1;j<(Wdth-1);j++) for(k=0;k<TROIS;k++) MtGrd[2][i][j]+=fabs(Img[k][i][j]-Img[k][i-1][j-1]);
    for(i=1;i<(Lgth-1);i++) for(j=1;j<(Wdth-1);j++) for(k=0;k<TROIS;k++) MtGrd[3][i][j]+=fabs(Img[k][i][j]-Img[k][i+1][j-1]);

    for(k=0;k<NBREGMAX;k++)
    {
        for(l=0;l<orient;l++)
        { nb=0;
            for(i=0;i<Lgth;i++) for(j=0;j<Wdth;j++) if (ImgSeg[i][j]==k)
            { nb++;
                VctHist[k][l]+=MtGrd[l][i][j]; }
            if (nb) VctHist[k][l]/=nb;
        }
    }

    for(k=0;k<NBREGMAX;k++)
    { moy=0.0;
        for(l=0;l<orient;l++) moy+=VctHist[k][l];
        moy/=(float)orient;

        var=0.0;
        for(l=0;l<orient;l++) var+=CARRE(VctHist[k][l]-moy);
        var/=(float)orient;

        VctOrientGrdVar[k]=var; }

    //>FreeMemory
    if (MtGrd) free_fmatrix_3d(MtGrd,orient);
    if (VctHist) free_fmatrix_2d(VctHist);
}

//----------------------------------------------------
// BilinearInterpolation                              
//----------------------------------------------------
float BilinearInterp(float* sig,int wdth,float col)  
{
    int col_inf,col_sup;
    float t,val;

    col_inf=(int)(col);
    col_sup=(int)(col+1);
    t=(col-col_inf)/(col_sup-col_inf);

    val=0.0;
    if ((col_inf>=0)&&(col_sup<wdth))
    { val=(1-t)*sig[col_inf];
        val+=t*sig[col_sup]; }

    return val;
}

//----------------------------------------------------
// MakeEvenSpectrum                                  
//----------------------------------------------------
void MakeEvenSpectrum(float* Spectrum,int Sz)
{
    int k;
    for(k=1;k<(Sz/2);k++) Spectrum[Sz-k]=Spectrum[k];
}

//----------------------------------------------------
// MakeOddSpectrum                                  
//----------------------------------------------------
void MakeOddSpectrum(float* Spectrum,int Sz)
{
    int k;
    for(k=1;k<(Sz/2);k++) Spectrum[Sz-k]=-Spectrum[k];
}

//----------------------------------------------------
// MakeStat
//----------------------------------------------------
void MkStat(float* Vct,int Sz)
{
    int k;
    float min,max;

    max=-10000000;
    min=10000000;
    for(k=0;k<Sz;k++)
    { if (Vct[k]<min) min=Vct[k];
        if (Vct[k]>max) max=Vct[k]; }

    printf(" > StatSound[%.2f::%.2f]         ",min,max);
    fflush(stdout);
}

//----------------------------------------------------
// MakeSpectrogram
//----------------------------------------------------
void MakeSpectrogram(float** ImgSeg,int Lgth,int Wdth,int row,int col, \
                     float** HueMix,float** GrdMix,float** SatLuNbMix,float* GrdMoy,\
                     float* OrientGrdVar,int HueNbBins,int GrdNbBins)
{
    int i,j,k,l,m;
    int LgthSig;
    int CurRegionLabel;
    float tmp;
    float* AudioSig=NULL;
    const int NB_ECH=16384;

    //>Include-ThermicScale
#include "ThermicScale.h"

    //>AllocateMemory
    float** ImgSpectrumModRow=fmatrix_allocate_2d(NB_ECH/2,Wdth);
    float** ImgSpectrumModCol=fmatrix_allocate_2d(Lgth,NB_ECH/2);
    float* AudioSigIm=fmatrix_allocate_1d(NB_ECH);
    float* AudioSigModSpectr=fmatrix_allocate_1d(NB_ECH);
    float** ImgSpectrumRow=fmatrix_allocate_2d(100,Wdth);
    float** ImgSpectrumCol=fmatrix_allocate_2d(Lgth,100);
    float*** ImgCSpectrumRow=fmatrix_allocate_3d(3,100,Wdth);
    float*** ImgCSpectrumCol=fmatrix_allocate_3d(3,Lgth,100);

    //>Column--Spectrogram
    for(k=0;k<Wdth;k++)
    { CurRegionLabel=ImgSeg[row][k];
        AudioSig=BuildSonif4(CurRegionLabel,HueMix,GrdMix,SatLuNbMix,GrdMoy,OrientGrdVar,HueNbBins,GrdNbBins,&LgthSig);
        for(l=0;l<NB_ECH;l++) AudioSigIm[l]=0.0;
        FFT1D(AudioSig,AudioSigIm,NB_ECH);
        Mod(AudioSigModSpectr,AudioSig,AudioSigIm,NB_ECH);
        for(l=0;l<(NB_ECH/2);l++) ImgSpectrumModRow[l][k]=AudioSigModSpectr[l]; }
    //Recal(ImgSpectrumModRow,NB_ECH/2,Wdth);
    //LogImg(ImgSpectrumModRow,NB_ECH/2,Wdth);
    //Recal(ImgSpectrumModRow,NB_ECH/2,Wdth);
    //SaveImagePgm((char*)"ImgSpectrCol",ImgSpectrumModRow,NB_ECH/2,Wdth);
    //exit(-1);

    //-DeltaNu=20Hz--entre:[0-2000Hz]
    float DN=20;
    for(k=0;k<Wdth;k++)
    { for(l=0;l<2000;l+=DN)
        { tmp=0.0;
            for(m=l;m<l+DN;m++) tmp+=ImgSpectrumModRow[m][k];
            tmp/=DN;
            if ((l/DN)<100) ImgSpectrumRow[int(l/DN)][k]=tmp; } }
    Recal(ImgSpectrumRow,100,Wdth);
    LogImg(ImgSpectrumRow,100,Wdth);
    Recal(ImgSpectrumRow,100,Wdth);
    //SaveImagePgm((char*)"ImgSpectrRow",ImgSpectrumRow,100,Wdth);
    for(i=0;i<100;i++) for(j=0;j<Wdth;j++)
    { ImgCSpectrumRow[0][i][j]=ThermicScale[(int)(ImgSpectrumRow[i][j])][0];
        ImgCSpectrumRow[1][i][j]=ThermicScale[(int)(ImgSpectrumRow[i][j])][1];
        ImgCSpectrumRow[2][i][j]=ThermicScale[(int)(ImgSpectrumRow[i][j])][2]; }
    SaveImagePpm((char*)"ImgCSpectrRow",ImgCSpectrumRow,100,Wdth);

    //>Row--Spectrogram
    for(k=0;k<Lgth;k++)
    { CurRegionLabel=ImgSeg[k][col];
        AudioSig=BuildSonif4(CurRegionLabel,HueMix,GrdMix,SatLuNbMix,GrdMoy,OrientGrdVar,HueNbBins,GrdNbBins,&LgthSig);
        for(l=0;l<NB_ECH;l++) AudioSigIm[l]=0.0;
        FFT1D(AudioSig,AudioSigIm,NB_ECH);
        Mod(AudioSigModSpectr,AudioSig,AudioSigIm,NB_ECH);
        for(l=0;l<(NB_ECH/2);l++) ImgSpectrumModCol[k][l]=AudioSigModSpectr[l]; }

    //-DeltaNu=20Hz--entre:[0-2000Hz]
    for(k=0;k<Lgth;k++)
    { for(l=0;l<2000;l+=DN)
        { tmp=0.0;
            for(m=l;m<l+DN;m++) tmp+=ImgSpectrumModCol[k][m];
            tmp/=DN;
            if ((l/DN)<100) ImgSpectrumCol[k][int(l/DN)]=tmp; } }
    Recal(ImgSpectrumCol,Lgth,100);
    LogImg(ImgSpectrumCol,Lgth,100);
    Recal(ImgSpectrumCol,Lgth,100);
    //SaveImagePgm((char*)"ImgSpectrCol",ImgSpectrumCol,Lgth,100);
    for(i=0;i<Lgth;i++) for(j=0;j<100;j++)
    { ImgCSpectrumCol[0][i][j]=ThermicScale[(int)(ImgSpectrumCol[i][j])][0];
        ImgCSpectrumCol[1][i][j]=ThermicScale[(int)(ImgSpectrumCol[i][j])][1];
        ImgCSpectrumCol[2][i][j]=ThermicScale[(int)(ImgSpectrumCol[i][j])][2]; }
    SaveImagePpm((char*)"ImgCSpectrCol",ImgCSpectrumCol,Lgth,100);

    //>FreeMemory
    if (ImgSpectrumRow) free_fmatrix_2d(ImgSpectrumRow);
    if (ImgSpectrumCol) free_fmatrix_2d(ImgSpectrumCol);
    if (ImgSpectrumModRow) free_fmatrix_2d(ImgSpectrumModRow);
    if (ImgSpectrumCol) free_fmatrix_2d(ImgSpectrumModCol);
    if (AudioSig)  free_fmatrix_1d(AudioSig);
    if (AudioSigIm)  free_fmatrix_1d(AudioSigIm);
    if (AudioSigModSpectr) free_fmatrix_1d(AudioSigModSpectr);
    if (ImgCSpectrumRow) free_fmatrix_3d(ImgCSpectrumRow,3);
    if (ImgCSpectrumCol) free_fmatrix_3d(ImgCSpectrumCol,3);
}


//------------------------------------------------------------------------------------------
// BuildSonif3
//------------------------------------------------------------------------------------------
float* BuildSonif3(int lab,float** HueMix,float** GrdMix,float** SatLuNbMix,float* GradMoy, \
                   int HueNbBins,int GrdNbBins,int* LgthSig)
{
    int k,l;
    float ValSat;
    float ValLum;
    float PosVol;

    const int DO=262;
    const int RE=294;
    const int MI=330;
    const int FA=349;
    const int SOL=392;
    const int LA=440;
    const int SI=494;
    const int OCTAVE=262;

    const int F_ECH=16384;
    const int NB_ECH=F_ECH;

    const int SVALTOSAT=16;
    const float MAX_AMP=1.0;
    const float MIN_AMP=-1.0;
    const int LVALTOLUM=16.0;
    const int DIV_GRD_NBBINS=NB_ECH/10;
    const float IMPGRDLOC=0.90;

    //>AllocateMemory_&_Init
    float* AudioSig=fmatrix_allocate_1d(F_ECH);
    float* VctAudioFreqR=fmatrix_allocate_1d(F_ECH);
    float* VctAudioFreqI=fmatrix_allocate_1d(F_ECH);
    float* VctTmp=fmatrix_allocate_1d(HueNbBins);
    float* VctGrdTmp=fmatrix_allocate_1d(2*GrdNbBins);
    for(k=0;k<F_ECH;k++)
    { AudioSig[k]=0.0;
        VctAudioFreqR[k]=0.0;
        VctAudioFreqI[k]=0.0; }
    for(k=0;k<HueNbBins;k++) VctTmp[k]=0.0;

    //:::::::::::::::::::
    //::::BUILD_AUDIO::::
    //:::::::::::::::::::
    //>-[1]-_Interpret_<Hue>_
    VctAudioFreqR[DO]=VctTmp[0]=HueMix[lab][0];
    VctAudioFreqR[RE]=VctTmp[1]=HueMix[lab][1];
    VctAudioFreqR[MI]=VctTmp[2]=HueMix[lab][2];
    VctAudioFreqR[FA]=VctTmp[3]=HueMix[lab][3];
    VctAudioFreqR[SOL]=VctTmp[4]=HueMix[lab][4];
    VctAudioFreqR[LA]=VctTmp[5]=HueMix[lab][5];
    VctAudioFreqR[SI]=VctTmp[6]=HueMix[lab][6];

    //>-[2]-_Interpret_<Saturation>_
    ValSat=(int)(SatLuNbMix[lab][0]/SVALTOSAT);
    for(l=0;l<ValSat;l++) for(k=0;k<HueNbBins;k++) VctAudioFreqR[DO+(OCTAVE*l)+k]=VctTmp[k];

    //>-[3]-_Interpret_<GradMoy>_(distorsion)
    for(k=0;k<(3*OCTAVE);k++) VctAudioFreqR[k]+=(GradMoy[lab]/100000.0);
    //printf("<%.3f>",GradMoy[lab]);

    //>Freq.Space_To_Temp.Space::>AudioSig
    for(k=0;k<F_ECH;k++) VctAudioFreqI[k]=0.0;
    MakeEvenSpectrum(VctAudioFreqR,F_ECH);
    IFFT1D(VctAudioFreqR,VctAudioFreqI,F_ECH);
    for(k=0;k<F_ECH;k++) AudioSig[k]=VctAudioFreqR[k];
    RecalSignal(AudioSig,F_ECH,MIN_AMP,MAX_AMP);

    //>-[4]-_Interpret_<Gradient>:>:Local_Volume
    for(k=0;k<GrdNbBins;k++) VctGrdTmp[k]=GrdMix[lab][k];
    for(k=0;k<GrdNbBins;k++) VctGrdTmp[GrdNbBins+k]=GrdMix[lab][GrdNbBins-1-k];
    for(k=0;k<NB_ECH;k++)
    { PosVol=VctGrdTmp[(int)((2*k)/DIV_GRD_NBBINS)];
        //PosVol=1.0-VctGrdTmp[(int)((2*k)/DIV_GRD_NBBINS)];
        AudioSig[k]=AudioSig[k]*((1.0-IMPGRDLOC)+IMPGRDLOC*PosVol); }

    //>-[5]-_Interpret_<Luminance>:>:Global_Volume
    ValLum=SatLuNbMix[lab][1]/LVALTOLUM;
    for(k=0;k<NB_ECH;k++) AudioSig[k]*=ValLum;

    //>Stat
    MkStat(AudioSig,NB_ECH);

    //>Free_Memory
    if (VctAudioFreqR) free_fmatrix_1d(VctAudioFreqR);
    if (VctAudioFreqI) free_fmatrix_1d(VctAudioFreqI);
    if (VctTmp) free_fmatrix_1d(VctTmp);
    if (VctGrdTmp) free_fmatrix_1d(VctGrdTmp);

    //>Return
    (*LgthSig)=F_ECH;
    return AudioSig;
}

//------------------------------------------------------------------------------------------
// BuildSonif
//
// Do(262 Hz) Re(294 Hz) Mi(330 Hz) Fa(349 Hz) Sol(392 Hz) La(440 Hz) Si(494 Hz)
// Si-Do (494-262=232 Hz)
// Do of the next octave 524 (262x2=524)
// Do of the next next octave 1048 (262x2^2)
//
// T=1 seconde => Delta_Nu=1Hz
// F_ECH=16384 Hz => [0--fmax=fe/2=8192 Hz]
//
// HueMix[NbRegMx][HueNbBins=7]
// GrdMix[NbRegMx][GRD_NBBINS=10]  //> 0.1 seconde by bin
// SatLuNbMix[NbRegMx][0]=Sat
//                    [1]=Lum
//                   ([2]=Nb)
//
//------------------------------------------------------------------------------------------
float* BuildSonif4(int lab,float** HueMix,float** GrdMix,float** SatLuNbMix,float* GradMoy, \
                   float* OrientGradVar,int HueNbBins,int GrdNbBins,int* LgthSig,bool octave,bool purity,bool distorsion,bool rythm)
{
    int k,l;
    float ValSat,ValLum;
    float PosVol;
    float DivLum;
    float tmp;
    int DIV_GRD_NBBINS=NB_ECH/GrdNbBins;

    //>AllocateMemory_&_Init
    float* AudioSig=fmatrix_allocate_1d(F_ECH);
    float* VctAudioFreqR=fmatrix_allocate_1d(F_ECH);
    float* VctAudioFreqI=fmatrix_allocate_1d(F_ECH);
    float* VctAudioFreqM=fmatrix_allocate_1d(F_ECH);
    float* VctAudioFreqP=fmatrix_allocate_1d(F_ECH);
    float* VctTmp=fmatrix_allocate_1d(HueNbBins);
    float* VctGrdTmp=fmatrix_allocate_1d(2*GrdNbBins);
    for(k=0;k<F_ECH;k++)
    { AudioSig[k]=0.0;
        VctAudioFreqR[k]=0.0;
        VctAudioFreqI[k]=0.0;
        VctAudioFreqM[k]=0.0;
        VctAudioFreqP[k]=0.0; }
    for(k=0;k<HueNbBins;k++)  VctTmp[k]=0.0;
    for(k=0;k<(2*GrdNbBins);k++) VctGrdTmp[k]=0.0;

    //:::::::::::::::::::
    //::::BUILD_AUDIO::::::::::::::::::::
    //:::::::::::::::::::
    const int DEBUG=0;

    //>-[1]-_Interpret_<Hue>_<Luminance>_
    //-----------------------------------
    if(octave)
    {
        // Interpret the luminance as level of octave
        ValLum=ceil(SatLuNbMix[lab][1]/32);
    }
    else
    {
        ValLum = 4;
    }

    DivLum=powf(2,(ValLum-4));
    if (DEBUG) printf("[Lab=%d][Lum=%.1f][ValLum=%.1f][DivLum=%.3f] ",lab,SatLuNbMix[lab][1],ValLum,DivLum); /**/
    for(k=0;k<HueNbBins;k++)
    {
        float val = HueMix[lab][k];
        VctAudioFreqM[ (int)(TabFreqNote[k]*DivLum) ]=val;
    }
    //if (1) for(k=0;k<NB_ECH;k++) if (VctAudioFreqM[k]) printf(" <%d:%.2f>",k,VctAudioFreqM[k]);

    //>-[2]-_Interpret_<Saturation>_
    //------------------------------
    if(purity)
    {
        // Interpret the saturation as purity
        ValSat=7-floor(SatLuNbMix[lab][0]/32);
    }
    else
    {
        ValSat = 0;
    }

    if (DEBUG) printf("[Sat=%.2f][ValLum=%.1f][ValSat=%.1f]> ",SatLuNbMix[lab][0],ValLum,ValSat); /**/
    //if (1) for(l=0;l<ValSat;l++) printf("(%d)",TabSat[(int)ValLum-1][l]);
    for(l=0;l<ValSat;l++) for(k=0;k<HueNbBins;k++)
    { DivLum=powf(2,(TabSat[(int)ValLum-1][l])-4);
        VctAudioFreqM[ (int)(TabFreqNote[k]*DivLum) ]=(1.0/(float(ValSat)))*HueMix[lab][k]; }
    //if (1) for(k=0;k<NB_ECH;k++) if (VctAudioFreqM[k]) printf(" <%d:%.2f>",k,VctAudioFreqM[k]);

    //>-[3a]-_Interpret_<GradMoy>_(Magnitude_Distorsion)_
    //---------------------------------------------------
    if (DEBUG) printf("<Gm:%.0f> ",GradMoy[lab]); /**/
    if(distorsion)
    {
        // magnitude distorsion
        for(k=0;k<(NB_ECH/2);k++) if (!VctAudioFreqM[k]) VctAudioFreqM[k]=((randomize()*GradMoy[lab]*MAGDIS)/NB_ECH);
    }

    //if (1) for(k=0;k<(NB_ECH/2);k++) printf(" <%d:%.5f>",k,VctAudioFreqM[k]);
    //if (1) for(k=0;k<300;k++) printf(" <%d:%.5f>",k,VctAudioFreqM[k]);

    //>-[3b]-_Interpret_<OrientGradVar>_(Phase_Distorsion)_
    //------------------------------------------------------
    if (DEBUG) printf("<Ogv:%.0f> ",OrientGradVar[lab]); /**/
    if(distorsion)
    {
        // phase distorsion
        for(k=0;k<(NB_ECH/2);k++) VctAudioFreqP[k]=(((randomize()-0.5)*OrientGradVar[lab]*PHADIS)/NB_ECH);
    }

    //if (1) for(k=0;k<(NB_ECH/2);k++) printf(" <%d:%.5f>",k,VctAudioFreqP[k]);
    //if (1) for(k=0;k<300;k++) printf(" <%d:%.5f>",k,VctAudioFreqP[k]);

    //>-[4]-_Freq.Space_To_Temp.Space::>AudioSig_
    //-------------------------------------------
    MakeEvenSpectrum(VctAudioFreqM,F_ECH);
    MakeOddSpectrum(VctAudioFreqP,F_ECH);
    for(k=0;k<NB_ECH;k++)
    { VctAudioFreqR[k]=VctAudioFreqM[k]*cos(VctAudioFreqP[k]);
        VctAudioFreqI[k]=VctAudioFreqM[k]*sin(VctAudioFreqP[k]); }
    IFFT1D(VctAudioFreqR,VctAudioFreqI,F_ECH);
    for(k=0;k<F_ECH;k++) AudioSig[k]=VctAudioFreqR[k];
    RecalSignal_(AudioSig,F_ECH,MIN_AMP,MAX_AMP);
    //if (1) SaveSignalDat((char*)"SigTest",AudioSig,F_ECH);

    //>-[5]-_Interpret_<Gradient>:>:Local_Volume_
    //-------------------------------------------
    if(rythm)
    {
        for(k=0;k<GrdNbBins;k++) VctGrdTmp[k]=GrdMix[lab][k];
        for(k=0;k<GrdNbBins;k++) VctGrdTmp[GrdNbBins+k]=GrdMix[lab][GrdNbBins-1-k];
        for(k=0;k<NB_ECH;k++)
        { PosVol=VctGrdTmp[(int)((2*k)/DIV_GRD_NBBINS)];
            tmp=((1.0-IMPGRDLOC)+IMPGRDLOC*PosVol);
            AudioSig[k]=tmp*AudioSig[k];
            //if (1) if (!(k%100)) printf("\n [%d::%.2f>%.2f]",k,PosVol,tmp);
        }
        //if (1) for(k=0;k<GrdNbBins;k++) printf("[%d:%.3f]",k,GrdMix[lab][k]);
        //if (1) for(k=0;k<(GrdNbBins*2);k++) printf("(%d:%.3f)",k,VctGrdTmp[k]);
        //if (1) SaveSignalDat((char*)"SigTest2",AudioSig,F_ECH);
    }

    if (DEBUG) printf("                              "); /**/

    //::::END_BUILD_AUDIO::::::::::::::::::::

    //>Stat
    //if (1) MkStat(AudioSig,NB_ECH);

    //>Free_Memory
    if (VctAudioFreqR) free_fmatrix_1d(VctAudioFreqR);
    if (VctAudioFreqI) free_fmatrix_1d(VctAudioFreqI);
    if (VctAudioFreqM) free_fmatrix_1d(VctAudioFreqM);
    if (VctAudioFreqP) free_fmatrix_1d(VctAudioFreqP);
    if (VctTmp) free_fmatrix_1d(VctTmp);
    if (VctGrdTmp) free_fmatrix_1d(VctGrdTmp);

    //>Return
    (*LgthSig)=F_ECH;
    return AudioSig;
}
