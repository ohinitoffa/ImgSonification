# Image Sonification
This repository countains the source code to reproduce the paper 
[A hierarchical Visual Feature-Based Approach For Image Sonification, O.K. Toffa and M. Mignotte, IEEE Transactions on Multimedia DOI: 10.1109/TMM.2020.2987710](https://ieeexplore.ieee.org/document/9070218)
<br>
This software run on linux platform and requires the software "aplay" and "sox" to be installed.

To compile
make -f Makefile

To visualize a sound signal use the script ViewSig.sh on the .dat file
for eg. ViewSig.sh SigTest.dat <br>
will generate SigTest.eps

For the comprehension of the sox software refer to 
http://www.iro.umontreal.ca/~mignotte/IFT3205/TpIFT3205_EffetsSonores_slides/slide004.html

\>\> ImgSonification


Usage   : MDSCCT [Img. in ppm format] <br>
 Options : value by default indicated in () <br>
            ----------------------------------------- <br>
            -o Octave <br>
            -p Purity <br>
            -d Distorsion <br>
            -r Rythm <br>
            -h BlackScreen <br>
            -z zoom <br>
			      -s save <br>
            - q quit <br>
			      <br>
            -------------

 Example : MDSCCT 198023

