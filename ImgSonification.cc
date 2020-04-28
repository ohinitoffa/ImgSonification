//------------------------------------------------------
// module  :
// author  : Mignotte Max, Toffa Ohini
// date    :
// version : 1.0
// language: C++
// note    :
//------------------------------------------------------
//

//------------------------------------------------
// INCLUDED FILES --------------------------------
//------------------------------------------------
#include "FunctionSonif.h"
#include <signal.h>
#include <sys/time.h>

//------------------------------------------------
// CONSTANTES ------------------------------------
//------------------------------------------------
const int OCTAVE=1;
const int PURITY=1;
const int DISTORSION=1;
const int RYTHM=1;
const int ZOOM=1;
const int NBVISUSEQ=1;
const int SAVE=0;
const int QUIT=0;
const int TIME_DELAY=1;

//------------------------------------------------
// VARIABLES GLOBALES ----------------------------
//------------------------------------------------
char Name_Img[100];
int octave = OCTAVE;
int purity = PURITY;
int distorsion = DISTORSION;
int rythm = RYTHM;
int blackScreen;
int zoom;
int flag_quit;
int flag_save;
int test=0;

//------------------------------------------------
// FOR X-WINDOWS ---------------------------------
//------------------------------------------------
#include <X11/Xutil.h>
#include <X11/Xos.h>

Display   *display;
int	  screen_num;
int 	  depth;
Window	  root;
Visual*	  visual;
GC	  gc;
struct itimerval timerval;
bool canPlay = false;

//-----------------------------------------------------------------------
// OPEN_DISPLAY()
//-----------------------------------------------------------------------
int open_display()
{
    if ((display=XOpenDisplay(NULL))==NULL)
    { printf("Connection impossible\n");
        return(-1); }

    else
    { screen_num=DefaultScreen(display);
        visual=DefaultVisual(display,screen_num);
        depth=DefaultDepth(display,screen_num);
        root=RootWindow(display,screen_num);
        return 0; }
}

void timer_callback(int sign)
{
    canPlay = true;
}

void timer_init()
{
    signal(SIGALRM, &timer_callback);
    struct timeval interval = { 0, 0 };
    struct timeval value = { TIME_DELAY, 0 };
    timerval.it_interval = interval;
    timerval.it_value = value;
    //setitimer(ITIMER_REAL, &timerval, 0);
}

//-----------------------------------------------------------------------
// FABRIQUE_WINDOW()
// Cette fonction cr?e une fenetre X et l'affiche ? l'?cran.
//-----------------------------------------------------------------------
Window fabrique_window(char *nom_fen,int x,int y,int width,int height,int zoom)
{
    Window                 win;
    XSizeHints      size_hints;
    XWMHints          wm_hints;
    XClassHint     class_hints;
    XTextProperty  windowName, iconName;

    char *name=nom_fen;

    if(zoom<0) { width/=-zoom; height/=-zoom; }
    if(zoom>0) { width*=zoom;  height*=zoom;  }

    win=XCreateSimpleWindow(display,root,x,y,width,height,1,0,255);

    size_hints.flags=PPosition|PSize|PMinSize;
    size_hints.min_width=width;
    size_hints.min_height=height;

    XStringListToTextProperty(&name,1,&windowName);
    XStringListToTextProperty(&name,1,&iconName);
    wm_hints.initial_state=NormalState;
    wm_hints.input=True;
    wm_hints.flags=StateHint|InputHint;
    class_hints.res_name=nom_fen;
    class_hints.res_class=nom_fen;

    XSetWMProperties(display,win,&windowName,&iconName,
                     NULL,0,&size_hints,&wm_hints,&class_hints);

    gc=XCreateGC(display,win,0,NULL);

    XSelectInput(display,win,ExposureMask|KeyPressMask|ButtonPressMask|
                 ButtonReleaseMask|ButtonMotionMask|PointerMotionHintMask|
                 StructureNotifyMask);

    XMapWindow(display,win);
    return(win);
}

//--------------------------------------------------------------------------
// CREE_XIMAGENB()
// Cr?e une XImage ? partir d'un tableau de float
// L'image peut subir un zoom.
//--------------------------------------------------------------------------
XImage* cree_XimageNB(float** mat,int z,int length,int width)
{
    int lgth,wdth,lig,col,zoom_col,zoom_lig;
    float somme;
    unsigned char	 pix;
    unsigned char* dat;
    XImage* imageX;

    /*Zoom positiv*/
    /*------------*/
    if (z>0)
    {
        lgth=length*z;
        wdth=width*z;

        dat=(unsigned char*)malloc(lgth*(wdth*4)*sizeof(unsigned char));
        if (dat==NULL)
        { printf("Impossible d'allouer de la memoire.");
            exit(-1); }

        for(lig=0;lig<lgth;lig=lig+z) for(col=0;col<wdth;col=col+z)
        {
            pix=(unsigned char)mat[lig/z][col/z];
            for(zoom_lig=0;zoom_lig<z;zoom_lig++) for(zoom_col=0;zoom_col<z;zoom_col++)
            {
                dat[((lig+zoom_lig)*wdth*4)+((4*(col+zoom_col))+0)]=pix;
                dat[((lig+zoom_lig)*wdth*4)+((4*(col+zoom_col))+1)]=pix;
                dat[((lig+zoom_lig)*wdth*4)+((4*(col+zoom_col))+2)]=pix;
                dat[((lig+zoom_lig)*wdth*4)+((4*(col+zoom_col))+3)]=pix;
            }
        }
    } /*--------------------------------------------------------*/

    /*Zoom negatifv*/
    /*------------*/
    else
    {
        z=-z;
        lgth=(length/z);
        wdth=(width/z);

        dat=(unsigned char*)malloc(lgth*(wdth*4)*sizeof(unsigned char));
        if (dat==NULL)
        { printf("Impossible d'allouer de la memoire.");
            exit(-1); }

        for(lig=0;lig<(lgth*z);lig=lig+z) for(col=0;col<(wdth*z);col=col+z)
        {
            somme=0.0;
            for(zoom_lig=0;zoom_lig<z;zoom_lig++) for(zoom_col=0;zoom_col<z;zoom_col++)
            { somme+=mat[lig+zoom_lig][col+zoom_col]; }

            somme/=(z*z);
            dat[((lig/z)*wdth*4)+((4*(col/z))+0)]=(unsigned char)somme;
            dat[((lig/z)*wdth*4)+((4*(col/z))+1)]=(unsigned char)somme;
            dat[((lig/z)*wdth*4)+((4*(col/z))+2)]=(unsigned char)somme;
            dat[((lig/z)*wdth*4)+((4*(col/z))+3)]=(unsigned char)somme;
        }
    } /*--------------------------------------------------------*/

    imageX=XCreateImage(display,visual,depth,ZPixmap,0,(char*)dat,wdth,lgth,16,wdth*4);
    return (imageX);
}

//---------------------------------------------------------------------------
// CREE_XIMAGECOL()
// Cr?e une XImage ? partir d'un tableau 3D de float
// L'image peut subir un zoom.
//--------------------------------------------------------------------------
XImage* cree_XimageCol(float*** matRVB,int z,int length,int width)
{
    int i;
    int lgth,wdth,lig,col,zoom_col,zoom_lig;
    float somme;
    float sum[3];
    unsigned char	 pixR,pixV,pixB;
    unsigned char  pixN;
    unsigned char* dat;
    XImage* imageX;

    /*Zoom positif*/
    /*------------*/
    if (z>0)
    {
        lgth=length*z;
        wdth=width*z;

        dat=(unsigned char*)malloc(lgth*(wdth*4)*sizeof(unsigned char));
        if (dat==NULL)
        { printf("Impossible d'allouer de la memoire.");
            exit(-1); }

        for(lig=0;lig<lgth;lig=lig+z) for(col=0;col<wdth;col=col+z)
        {
            pixR=(unsigned char)matRVB[0][lig/z][col/z];
            pixV=(unsigned char)matRVB[1][lig/z][col/z];
            pixB=(unsigned char)matRVB[2][lig/z][col/z];
            somme=(1.0/3.0)*(pixR+pixV+pixB);
            pixN=(unsigned char)somme;

            for(zoom_lig=0;zoom_lig<z;zoom_lig++) for(zoom_col=0;zoom_col<z;zoom_col++)
            {
                dat[((lig+zoom_lig)*wdth*4)+((4*(col+zoom_col))+0)]=pixB;
                dat[((lig+zoom_lig)*wdth*4)+((4*(col+zoom_col))+1)]=pixV;
                dat[((lig+zoom_lig)*wdth*4)+((4*(col+zoom_col))+2)]=pixR;
                dat[((lig+zoom_lig)*wdth*4)+((4*(col+zoom_col))+3)]=pixN;
            }
        }
    } /*--------------------------------------------------------*/

    /*Zoom negatif*/
    /*------------*/
    else
    {
        z=-z;
        lgth=(length/z);
        wdth=(width/z);

        dat=(unsigned char*)malloc(lgth*(wdth*4)*sizeof(unsigned char));
        if (dat==NULL)
        { printf("Impossible d'allouer de la memoire.");
            exit(-1); }

        for(lig=0;lig<(lgth*z);lig=lig+z) for(col=0;col<(wdth*z);col=col+z)
        {
            sum[0]=sum[1]=sum[2]=0.0;

            for(i=0;i<3;i++)
                for(zoom_lig=0;zoom_lig<z;zoom_lig++) for(zoom_col=0;zoom_col<z;zoom_col++)
                    sum[i]+=matRVB[i][lig+zoom_lig][col+zoom_col];

            for(i=0;i<3;i++)  sum[i]/=(z*z);

            dat[((lig/z)*wdth*4)+((4*(col/z))+0)]=(unsigned char)sum[2];
            dat[((lig/z)*wdth*4)+((4*(col/z))+1)]=(unsigned char)sum[1];
            dat[((lig/z)*wdth*4)+((4*(col/z))+2)]=(unsigned char)sum[0];
            dat[((lig/z)*wdth*4)+((4*(col/z))+3)]=(unsigned char)sum[1];
        }
    } /*--------------------------------------------------------*/

    imageX=XCreateImage(display,visual,depth,ZPixmap,0,(char*)dat,wdth,lgth,16,wdth*4);
    return (imageX);
}

//---------------------------------------------------------
// Enregistre les arguments de la ligne de commande
//---------------------------------------------------------
void read_arguments(int argc,char** argv)
{
    int i;

    //Options
    if (argc<2)
    { printf("\n Usage %s [Img. in ppm format]",argv[0]);
        printf("\n Options : value by default indicated in []");
        printf("\n            -----------------------------------------");
        printf("\n               -o Octave   >[%d]",OCTAVE);
        printf("\n               -p Purity   >[%d]",PURITY);
        printf("\n               -d Distorsion    >[%d]",DISTORSION);
        printf("\n               -r Rythm    >[%d]",RYTHM);
        printf("\n               -h BlackScreen      >[%d]",0);
        printf("\n            -------------");
        printf("\n            -z zoom (%d)",ZOOM);
        printf("\n            -s save (Sans)");
        printf("\n            -q quit (Sans)");
        printf("\n            -------------");
        printf("\n");
        printf("\n  Example: %s 198023",argv[0]);
        printf("\n  Example: %s 12003",argv[0]);
        printf("\n  Example: %s 86000",argv[0]);
        printf("\n  Example: %s 134052",argv[0]);
        printf("\n  Example: %s 277095",argv[0]);
        printf("\n\n\n\n\n");
        exit(-1); }

    //Chargement Fichier
    if (argc>1) strcpy(Name_Img,argv[1]);


    //Boucle
    for(i=2;i<argc;i++)
    {
        switch(argv[i][1])
        {
        case 'o': octave=atoi(argv[++i]);      break;
        case 'p': purity=atoi(argv[++i]);      break;
        case 'd': distorsion=atoi(argv[++i]);      break;
        case 'r': rythm=atoi(argv[++i]);      break;
        case 'h': blackScreen=atoi(argv[++i]);   break;
        case 'z': zoom=atoi(argv[++i]);        break;
        case 's': flag_save=1;                 break;
        case 'q': flag_quit=1;                 break;
        case 't': test=atoi(argv[++i]);        break;
        }
    }
}


//---------------------------------------------------------
// MakeCross
//---------------------------------------------------------
void MakeCross(float*** Img,float*** Img_,int lgth,int wdth,int PosRow,int PosCol)
{
    int i,j,k,l;

    if(blackScreen)
    {
        for(k=0;k<TROIS;k++) for(i=0;i<lgth;i++) for(j=0;j<wdth;j++) Img[k][i][j]=0;
    }
    else
    {
        for(k=0;k<TROIS;k++) for(i=0;i<lgth;i++) for(j=0;j<wdth;j++) Img[k][i][j]=Img_[k][i][j];
    }

    for(k=-2;k<=2;k++) for(l=-2;l<=2;l++)
        if (((PosRow+k)>=0)&&((PosRow+k)<lgth)&&((PosCol+k)>=0)&&((PosCol+k)<wdth))
        { Img[0][PosRow+k][PosCol+l]=255;
            Img[1][PosRow+k][PosCol+l]=0;
            Img[2][PosRow+k][PosCol+l]=0; }
}

//----------------------------------------------------------
//----------------------------------------------------------
// PROGRAMME PRINCIPAL -------------------------------------
//----------------------------------------------------------
//----------------------------------------------------------
int main(int argc,char** argv)
{
    int k;
    XEvent ev;

    //>Pour Xwindow
    Window win_ppictureRGB=0;
    Window win_ppictureSEG=0;
    XImage *x_ppictureRGB;
    XImage *x_ppictureSEG;
    char nomfen_ppictureRGB[100];
    char nomfen_ppictureSEG[100];

    //>Initialisation
    blackScreen=0;
    zoom=ZOOM;
    flag_quit=QUIT;
    flag_save=SAVE;

    //>Read Arguments
    read_arguments(argc,argv);

    //>Display Options
    printf("\n\n");
    printf("\n SONIFICATION");
    printf("\n ----------------------");
    printf("\n Info.    : Image[%s]",Name_Img);
    printf("\n ---------------------");
    printf("\n Octave  : [%d]",octave);
    printf("\n Purity  : [%d]",purity);
    printf("\n Distorsion  : [%d]",distorsion);
    printf("\n Rythm   : [%d]",rythm);
    printf("\n Black Screen: [%d]",blackScreen);
    printf("\n ----------------");
    printf("\n Zoom     : [%d]",zoom);
    printf("\n ----------------");
    if (flag_quit)       printf("\n -> Avec Quit");
    if (flag_save)       printf("\n -> Avec Save");
    printf("\n\n");
    fflush(stdout);


    //------
    //------------------------------------------
    //->-[1] Read Imgs
    //------------------------------------------
    //------
    //>Load_Images
    int Length,Width;
    char ImgPathAndName[NCHAR];

    strcpy(ImgPathAndName,"BSDS300_IMAGES/");
    strcat(ImgPathAndName,Name_Img);
    strcat(ImgPathAndName,".ppm");
    GetLengthWidth(ImgPathAndName,&Length,&Width);
    float*** ImgRGB=fmatrix_allocate_3d(TROIS,Length,Width);
    float*** ImgRGB_=fmatrix_allocate_3d(TROIS,Length,Width);
    LoadImagePpm(ImgPathAndName,ImgRGB,Length,Width);
    CopyMat(ImgRGB,ImgRGB_,Length,Width);
    //SaveImagePpm((char*)"ImgRGB",ImgRGB,Length,Width);

    float** ImgSEG=fmatrix_allocate_2d(Length,Width);
    strcpy(ImgPathAndName,"BSDS300_SEGMENTATIONS/");
    strcat(ImgPathAndName,Name_Img);
    strcat(ImgPathAndName,"_1.pgm");
    if(test) LoadImagePgm("BSDS300_SEGMENTATIONS/square.pgm",ImgSEG,Length,Width);
    else LoadImagePgm(ImgPathAndName,ImgSEG,Length,Width);

    //>Conversion_HSL_
    float*** ImgHSL=fmatrix_allocate_3d(TROIS,Length,Width);
    ImgRGBToImgHsl(ImgRGB,ImgHSL,Length,Width);
    // Recal(ImgHSL,Length,Width); Le recalage pousse les sons dans les hautes fréquences
    //SaveImagePpm((char*)"ImgHSL",ImgHSL,Length,Width);

    //>SONIFICATION::Pre_Comnpute_================
    printf("\n PreCompute...");
    const int NBREGMAX=256;
    const int HUE_NBBINS=7;
    const int GRD_NBBINS=10;
    int CurRegionLabel,OldRegionLabel;
    float** HueMixture=fmatrix_allocate_2d(NBREGMAX,HUE_NBBINS);
    float** GrdMixture=fmatrix_allocate_2d(NBREGMAX,GRD_NBBINS);
    float** SatLuNbMixture=fmatrix_allocate_2d(NBREGMAX,3);
    float* GrdMoy=fmatrix_allocate_1d(NBREGMAX);
    float* OrientGrdVar=fmatrix_allocate_1d(NBREGMAX);
    float* AudioSignal;
    int LgthSig;
    ComputeHueMixture(ImgHSL[0],ImgSEG,Length,Width,HueMixture,NBREGMAX,HUE_NBBINS);
    ComputeGrdMixture(ImgRGB,ImgSEG,Length,Width,GrdMixture,NBREGMAX,GRD_NBBINS);
    ComputeSatLuMixture(ImgHSL,ImgSEG,Length,Width,SatLuNbMixture,NBREGMAX);
    ComputeGrdMoy(ImgRGB,ImgSEG,Length,Width,GrdMoy,NBREGMAX);
    ComputeOrientGradVar(ImgRGB,ImgSEG,Length,Width,OrientGrdVar,NBREGMAX);
    //=============================================

    //>DEBUG
    if (0) for(int k=0;k<NBREGMAX;k++)
    { printf("\n %d>",k);
        //for(int l=0;l<HUE_NBBINS;l++) printf("[%.2f]",HueMixture[k][l]);
        //for(int l=0;l<GRD_NBBINS;l++) printf("[%.2f]",GrdMixture[k][l]);
        //for(int l=0;l<3;l++) printf("[%.0f]",SatLuNbMixture[k][l]);
        //printf("[%.0f]",GrdMoy[k]);
        printf("[%.0f]",OrientGrdVar[k]);
    }

    //>DEBUG
    if (0)
    { MakeSpectrogram(ImgSEG,Length,Width,320,110,
                      HueMixture,GrdMixture,SatLuNbMixture,GrdMoy,OrientGrdVar,HUE_NBBINS,GRD_NBBINS);
        exit(-1); }

    //------
    //------------------------------------------
    //->-[2] VisuXWindows
    //------------------------------------------
    //------
    //>Start_Graphical_Session
    if (open_display()<0) printf(" Impossible to open a graphical session");
    sprintf(nomfen_ppictureRGB,"Image : %s",Name_Img);
    win_ppictureRGB=fabrique_window(nomfen_ppictureRGB,10,10,Width,Length,zoom);
    if(blackScreen)
    {
        for(k=0;k<TROIS;k++) for(int i=0;i<Length;i++) for(int j=0;j<Width;j++) ImgRGB[k][i][j]=0;
        x_ppictureRGB=cree_XimageCol(ImgRGB,zoom,Length,Width);
        for(int k=0;k<100;k++)
        { XPutImage(display,win_ppictureRGB,gc,x_ppictureRGB,0,0,0,0,x_ppictureRGB->width,x_ppictureRGB->height); }
    }
    else
    {
        win_ppictureSEG=fabrique_window(nomfen_ppictureSEG,Width+20,10,Width,Length,zoom);
        x_ppictureRGB=cree_XimageCol(ImgRGB,zoom,Length,Width);
        x_ppictureSEG=cree_XimageNB(ImgSEG,zoom,Length,Width);
        for(int k=0;k<100;k++)
        { XPutImage(display,win_ppictureRGB,gc,x_ppictureRGB,0,0,0,0,x_ppictureRGB->width,x_ppictureRGB->height);
            XPutImage(display,win_ppictureSEG,gc,x_ppictureSEG,0,0,0,0,x_ppictureSEG->width,x_ppictureSEG->height);  }
    }

    //>Red_Cross_On_Image
    int PosCurRow, PosCurCol;
    int STEP=5;
    PosCurRow=Length/2;
    PosCurCol=Width/2;
    CurRegionLabel=ImgSEG[PosCurRow][PosCurCol];
    printf("\n > Press Keys For [col:4...6][time:9...3]\n");
    printf("\r [Col::%d/%d][Row::%d/%d][Reg:%d]  ",PosCurCol,Width,PosCurRow,Length,CurRegionLabel);
    fflush(stdout);
    AudioSignal=NULL;

    timer_init();
    for(k=0;!flag_quit;k++)
    { XNextEvent(display,&ev);
        if (ev.type == KeyPress)
        { if  (ev.xkey.keycode == 0x53) PosCurCol-=STEP;
            if  (ev.xkey.keycode == 0x55) PosCurCol+=STEP;
            if  (ev.xkey.keycode == 0x50) PosCurRow-=STEP;
            if  (ev.xkey.keycode == 0x58) PosCurRow+=STEP;
            if  (ev.xkey.keycode == 0x09 ) break;
            PosCurCol=PosCurCol%Width;
            PosCurRow=PosCurRow%Length;
            if (PosCurCol<0) PosCurCol+=Width;
            if (PosCurRow<0) PosCurRow+=Length;
            CurRegionLabel=ImgSEG[PosCurRow][PosCurCol];
            printf("\r [Col::%d/%d][Row::%d/%d][Reg:%d]  ",PosCurCol,Width,PosCurRow,Length,CurRegionLabel);
            fflush(stdout);
            MakeCross(ImgRGB,ImgRGB_,Length,Width,PosCurRow,PosCurCol);
            x_ppictureRGB=cree_XimageCol(ImgRGB,zoom,Length,Width);
            XPutImage(display,win_ppictureRGB,gc,x_ppictureRGB,0,0,0,0,x_ppictureRGB->width,x_ppictureRGB->height);
            XDestroyImage(x_ppictureRGB);
            if  (ev.xkey.keycode == 0x54)
            { SaveSignalDat((char*)"AudioSig",AudioSignal,LgthSig);
                printf("<row=%d:col=%d>",PosCurRow,PosCurCol);
                printf("\n"); }
        }

        if (OldRegionLabel!=CurRegionLabel)
        {
            if (AudioSignal) free_fmatrix_1d(AudioSignal);
            //AudioSignal=BuildSonif3(CurRegionLabel,HueMixture,GrdMixture,SatLuNbMixture,GrdMoy,HUE_NBBINS,GRD_NBBINS,&LgthSig);
            //AudioSignal=BuildSonif4(CurRegionLabel,HueMixture,GrdMixture,SatLuNbMixture,GrdMoy,OrientGrdVar,HUE_NBBINS,GRD_NBBINS,&LgthSig);
            AudioSignal=BuildSonif4(CurRegionLabel,HueMixture,GrdMixture,SatLuNbMixture,GrdMoy,OrientGrdVar,HUE_NBBINS,GRD_NBBINS,&LgthSig,octave,purity,distorsion,rythm);
            SaveSignalDatWav((char*)"AudioSig",AudioSignal,LgthSig,LgthSig);
            OldRegionLabel=CurRegionLabel;
            canPlay = false;
            system("aplay -q AudioSig.wav&");
            setitimer(ITIMER_REAL, &timerval, 0);
        }
        else if(canPlay)
        {
            canPlay = false;
            system("aplay -q AudioSig.wav&");
            setitimer(ITIMER_REAL, &timerval, 0);
        }
    }

    //------
    //-------------
    // End
    //-------------
    //------
    //FreeMemory
    if (HueMixture) free_fmatrix_2d(HueMixture);
    if (GrdMixture) free_fmatrix_2d(GrdMixture);
    if (SatLuNbMixture) free_fmatrix_2d(SatLuNbMixture);
    if (AudioSignal) free_fmatrix_1d(AudioSignal);
    if (GrdMoy) free_fmatrix_1d(GrdMoy);
    if (OrientGrdVar) free_fmatrix_1d(OrientGrdVar);

    //ReturnWithoutProblem
    printf("\n\n\n C'est fini... \n\n\n");
    fflush(stdout);
    return 0;
}

