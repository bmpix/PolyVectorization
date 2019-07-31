#ifndef _PARAMS_H_
#define _PARAMS_H_

//If you see gaps in the final vectorization, or some separate curves get connected, this is the constant to tune (0-255). However, an easier way is to adjust the contrast of your input image.
const double BACKGROUND_FOREGROUND_THRESHOLD = 90.0;

//The more regularization, the more orthogonal is the frame field => the sharper the junctions, but narrow junctions might become wabbly.
const double FRAME_FIELD_REGULARIZER_WEIGHT = 0.1;
//Adjust only if you have crazily noisy image
const double FRAME_FIELD_SMOOTHNESS_WEIGHT = 50.0;
//If you see if removing little branches in the vectorization, try tuning that, but just a bit
const double PRUNE_SHORT_BRANCHES_RATIO = 0.75;
//If you have a very low-res image and you want to preserve tiny holes (white areas), tune this.
const int MAX_NUMBER_OF_WHITE_PIXELS_IN_A_CONTRACTIBLE_LOOP = 4;


#endif
