//
// Created by awbrenn on 1/20/16.
//

#ifndef ADVECTION_CFD_H
#define ADVECTION_CFD_H


class cfd
{

  public:
    // constructors/destructors
    cfd(int nx, int ny, float dx);
    ~cfd();

    // public methods
    void bilinearlyInterpolate(const float x, const float y);
    void advect(const float dt = (float)(1.0/24.0));

    // getters
    int getNx()           const { return Nx; }
    int getNy()           const { return  Ny; }
    float getDx()         const { return Dx; }
    float getGravityX()   const { return gravityX;}
    float getGravityY()   const { return gravityY;}
    float* getDensity1()  const { return density1; }
    float* getDensity2()  const { return density2; }
    float* getVelocity1() const { return velocity1; }
    float* getVelocity2() const { return velocity2; }
    float* getColor1()    const { return color1; }
    float* getColor2()    const { return color2; }
    float* getDensitySourceField() { return densitySourceField; }

    // setters
    void setDensitySourceField(float* dsrc) { densitySourceField = dsrc; }

    // indexing
    int dIndex(int i, int j)        const { return i+Nx*j; }
    int vIndex(int i, int j, int c) const { return (i+Nx*j)*2+c; }
    int cIndex(int i, int j, int c) const { return (i+Nx*j)*3+c; }

  private:
    int     Nx, Ny;
    float   Dx;
    float   gravityX, gravityY;
    float   *density1, *density2;
    float   *velocity1, *velocity2;
    float   *color1, *color2;
    float   *densitySourceField;
};

#endif //ADVECTION_CFD_H
