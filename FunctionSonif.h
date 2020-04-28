//------------------------------------------------------
// module  : FunctionSonif.h
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
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <new>

//------------------------------------------------
// CONSTANTES & DEFINITIONS ----------------------
//------------------------------------------------
#define NCHAR 200
#define NBCHAR 200
#define PI  3.14159
#define GREY_LEVEL 255
#define TROIS 3

#define FFT   1
#define IFFT -1
#define FFT2D 2

#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr
#define SQUARE(X) ((X)*(X))
#define CARRE(X)  ((X)*(X))
#define MAX(a,b)   (((a)>(b))?(a):(b))
#define MIN(a,b)   (((a)>(b))?(b):(a))

//------------------------------------------------ */
// PROTOTYPE ------------------------------------- */
//------------------------------------------------ */
//>Matrix_Allocation
float*    fmatrix_allocate_1d(int);
float**   fmatrix_allocate_2d(int,int);
float***  fmatrix_allocate_3d(int,int,int);
void  free_fmatrix_1d(float*);
void  free_fmatrix_2d(float**);
void  free_fmatrix_3d(float***,int);

//>Matrix_Gestion 
void  CopyMat(float***,float***,int,int);
void  CopyMat(float**,float**,int,int);
void Recal(float**,int,int);
void Recal(float***,int,int);
void LogImg(float**,int,int);

//>Load/Save File
void GetLengthWidth(char*,int*,int*);
void LoadImagePgm(char*,float**,int,int);
void LoadImagePpm(char*,float***,int,int);
void LoadPrintImagePpm(char*,int,int);
void SaveImagePgm(char*,float**,int,int);
void SaveImagePpm(char*,float***,int,int);

//>Fourier
void fourn(float*,unsigned long*,int,int);
void FFTDD(float**,float**,int,int);
void IFFTDD(float**,float**,int,int);
void FFT1D(float*,float*,int);
void IFFT1D(float*,float*,int);
void Mod(float**,float**,float**,int,int);
void Mod(float*,float*,float*,int,int);
void Mult(float**,float,int,int);
void MultMatrix(float**,float**,float**,float**,float**,float**,int,int);
void SquareMatrix(float**,float**,float**,float**,int,int);
void CenterImg(float**,int,int);
void CenterVct(float*,int);

//>Lecture/Sauve Signal
float* LoadSignalDat(char*,int*); 
void SaveSignalDat(char*,float*,int);
void SaveSignalDatWav(char*,float*,int,int);

//>SpaceColor
void ImgRGBToImgHsl(float***,float***,int,int);
void RGBToHSV(float,float,float,float*,float*,float*);
void RGBToHSL(float,float,float,float*,float*,float*);
void RGBToHSL2(float,float,float,float*,float*,float*);

//>Aleat
float randomize(void);

//>Sonification
void RecalSignal(float*,int,float,float);
void ValHist1D(float**,int,int,float**,int,float*,int,int);
void GradHist1D(float***,int,int,float**,int,float*,int,int);
void ComputeHueMixture(float**,float**,int,int,float**,int,int);
void ComputeGrdMixture(float***,float**,int,int,float**,int,int);
void ComputeSatLuMixture(float***,float**,int,int,float**,int);
void ComputeGrdMoy(float***,float**,int,int,float*,int);
void ComputeOrientGradVar(float***,float**,int,int,float*,int);
float BilinearInterp(float*,int,float);
void MakeEvenSpectrum(float*,int);
void MakeOddSpectrum(float*,int);
void MkStat(float*,int);
void MakeSpectrogram(float**,int,int,int,int,float**,float**,float**,float*,float*,int,int);
float* BuildSonif3(int,float**,float**,float**,float*,int,int,int*);
float* BuildSonif4(int,float**,float**,float**,float*,float*,int,int,int*,bool octave=true,bool purity=true,bool distorsion=true,bool rythm=true);
