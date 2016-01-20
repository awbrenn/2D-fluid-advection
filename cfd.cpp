//
// Created by awbrenn on 1/20/16.
//
#include <cmath>
#include "cfd.h"
#include "cfdUtility.h"

// dealing with negative results
int mod(int a, int b)
{
  return (a%b+b)%b;
}

cfd::cfd(const int nx, const int ny, const float dx)
{
  Nx = nx;
  Ny = ny;
  Dx = dx;
  density1 = new float[Nx*Ny]();
  density2 = new float[Nx*Ny]();
  velocity1 = new float[Nx*Ny*2]();
  velocity2 = new float[Nx*Ny*2]();
  color1 = new float[Nx*Ny*3]();
  color2 = new float[Nx*Ny*3]();
}


cfd::~cfd()
{
  delete density1;
  delete density2;
  delete velocity1;
  delete velocity2;
  delete color1;
  delete color2;
}


void cfd::bilinearlyInterpolate(float x, float y)
{
  // get index if sample
  const int i = mod((int) (x/Dx), Nx);
  const int j = mod((int) (y/Dx), Ny);

  // get weights of samples
  const float ax = std::abs(x/Dx - ((int)(x/Dx)));
  const float ay = std::abs(y/Dx - ((int)(y/Dx)));
  const float w1 = (1-ax) * (1-ay);
  const float w2 = ax * (1-ay);
  const float w3 = (1-ax) * ay;
  const float w4 = ax * ay;

  density2[dIndex(i,j)]    = density1[dIndex(i,j)]                           * w1 +
                             density1[dIndex((i + 1) % Nx, j)]               * w2 +
                             density1[dIndex(i, (j + 1) % Ny)]               * w3 +
                             density1[dIndex((i + 1) % Nx, (j + 1) % Ny)]    * w4;
  velocity2[vIndex(i,j,0)] = velocity1[vIndex(i,j,0)]                        * w1 +
                             velocity1[vIndex((i + 1) % Nx, j,0)]            * w2 +
                             velocity1[vIndex(i, (j + 1) % Ny,0)]            * w3 +
                             velocity1[vIndex((i + 1) % Nx, (j + 1) % Ny,0)] * w4;
  velocity2[vIndex(i,j,1)] = velocity1[vIndex(i,j,1)]                        * w1 +
                             velocity1[vIndex((i + 1) % Nx, j,1)]            * w2 +
                             velocity1[vIndex(i, (j + 1) % Ny,1)]            * w3 +
                             velocity1[vIndex((i + 1) % Nx, (j + 1) % Ny,1)] * w4;
  color2[cIndex(i,j,0)]    = color1[cIndex(i,j,0)]                           * w1 +
                             color1[cIndex((i + 1) % Nx, j,0)]               * w2 +
                             color1[cIndex(i, (j + 1) % Ny,0)]               * w3 +
                             color1[cIndex((i + 1) % Nx, (j + 1) % Ny,0)]    * w4;
  color2[cIndex(i,j,1)]    = color1[cIndex(i,j,1)]                           * w1 +
                             color1[cIndex((i + 1) % Nx, j,1)]               * w2 +
                             color1[cIndex(i, (j + 1) % Ny,1)]               * w3 +
                             color1[cIndex((i + 1) % Nx, (j + 1) % Ny,1)]    * w4;
  color2[cIndex(i,j,2)]    = color1[cIndex(i,j,2)]                           * w1 +
                             color1[cIndex((i + 1) % Nx, j,2)]               * w2 +
                             color1[cIndex(i, (j + 1) % Ny,2)]               * w3 +
                             color1[cIndex((i + 1) % Nx, (j + 1) % Ny,2)]    * w4;
}


void cfd::advect(const float dt)
{
  float x, y;

  for (int j=0; j<Ny; ++j)
  {
    for (int i=0; i<Nx; ++i)
    {
      x = i*Dx - velocity1[vIndex(i,j,0)]*dt;
      y = j*Dx - velocity1[vIndex(i,j,1)]*dt;
      bilinearlyInterpolate(x,y);
    }
  }

  swapFloatPointers(&density1, &density2);
  swapFloatPointers(&velocity1, &velocity2);
  swapFloatPointers(&color1, &color2);
}