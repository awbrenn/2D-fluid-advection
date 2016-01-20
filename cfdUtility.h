//
// Created by awbrenn on 1/20/16.
//

#ifndef ADVECTION_CFDUTILITY_H
#define ADVECTION_CFDUTILITY_H

void swapFloatPointers(float** a, float** b)
{
  float* temp = *a;
  *a = *b;
  *b = temp;
}

#endif //ADVECTION_CFDUTILITY_H
