//------------------------------------------------
//
//  cfd_paint
//
//
//-------------------------------------------------

//-------------------------------------------------
//
//  usage:
//
//  cfd_paint is an interactive paint program
//  in which the user paints density, color, and
//  or divergence sources that flow using
//  computational fluid dynamics and react with 
//  obstructions in the space.
//
//  There are two paint modes.  Typing 'o' puts the
//  program in obstruction painting mode. When you
//  hold down the left mouse button and paint, you
//  will see a black obstruction painted.  This 
//  obstruction may be any shape.
//
//  Typing 's' puts the program in source painting 
//  mode.  Now painting with the left mouse button
//  down injects density into the simulation.
//  The flow it produces evolves as you
//  continue to paint.  The flow bounces off any
//  obstructions that have been painted or are
//  subsequently painted.
//
//  Typing 'b' clears all obstructions, flow, density,
//  and color.
//
//  Typing '=' and '-' brightens and darkens the display.
//
//  Pressing the spacebar starts and stops the flow 
//  evolution. While the evolution is stopped, you
//  can continue painting obstructions.
//
//
//
//-------------------------------------------------




#include <cmath>
#include "CmdLineFind.h"

#include <GL/gl.h>   // OpenGL itself.
#include <GL/glu.h>  // GLU support library.
#include <GL/glut.h> // GLUT support library.

#include <iostream>
#include <OpenImageIO/imageio.h>
#include <omp.h>


using namespace std;
using namespace lux;
OIIO_NAMESPACE_USING

int iwidth, iheight, size;
float *display_map;
float *baseimage;

int paint_mode;
enum{ PAINT_OBSTRUCTION, PAINT_SOURCE, PAINT_DIVERGENCE, PAINT_COLOR };

bool toggle_animation_on_off;

float scaling_factor;

#define BRUSH_SIZE 11
float obstruction_brush[BRUSH_SIZE][BRUSH_SIZE];
float source_brush[BRUSH_SIZE][BRUSH_SIZE];

int xmouse_prev, ymouse_prev;

// dealing with negative results
int mod(int a, int b)
{
  return (a%b+b)%b;
}

void swapFloatPointers(float** a, float** b)
{
  float* temp = *a;
  *a = *b;
  *b = temp;
}

class cfd
{

  public:
    // constructors/destructors
    cfd(int nx, int ny, float dx);
    ~cfd();

    // public methods
    void bilinearlyInterpolate(float x, float y);
    void advect(const float dt);

    // getters
    int getNx()           const { return Nx; }
    int getNy()           const { return Ny; }
    float getDx()         const { return Dx; }
    float* getDensity1()  { return density1; }
    float* getDensity2()  { return density2; }
    float* getVelocity1() { return velocity1; }
    float* getVelocity2() { return velocity2; }
    float* getColor1()    { return color1; }
    float* getColor2()    { return color2; }

    // indexing
    int dIndex(int i, int j)        const { return i+Nx*j; }
    int vIndex(int i, int j, int c) const { return (i+Nx*j)*2+c; }
    int cIndex(int i, int j, int c) const { return (i+Nx*j)*3+c; }

  private:
    int     Nx, Ny;
    float   Dx;
    float   *density1;
    float   *density2;
    float   *velocity1;
    float   *velocity2;
    float   *color1;
    float   *color2;
};


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
  const float ax = abs(x/Dx - ((int)(x/Dx)));
  const float ay = abs(y/Dx - ((int)(y/Dx)));
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



////////  OpenImageIO reader

void readOIIOImage( const char* fname, float* img  )
{
  int xres, yres, channels;
  ImageInput *in = ImageInput::create (fname);
  if (! in) {return;}
  ImageSpec spec;
  in->open (fname, spec);
  xres = spec.width;
  yres = spec.height;
  channels = spec.nchannels;
  float* pixels = new float[xres*yres*channels];

  in->read_image (TypeDesc::FLOAT, pixels);
  long index = 0;
  for( int j=0;j<yres;j++)
  {
    for( int i=0;i<xres;i++ )
    {
      for( int c=0;c<channels;c++ )
      {
        img[ (i + xres*(yres - j - 1))*channels + c ] = pixels[index++];
      }
    }
  }

  in->close ();
  delete in;
}


//--------------------------------------------------------
//
//  Initialization routines
//
//  
// Initialize all of the fields to zero
void Initialize( float *data, int size, float value )
{
#pragma omp parallel for
   for(int i=0;i<size;i++ ) { data[i] = value; }
}


void InitializeBrushes()
{
  int brush_width = (BRUSH_SIZE-1)/2;
  for( int j=-brush_width;j<=brush_width;j++ )
  {
    int jj = j + brush_width;
    float jfactor =  (float(brush_width) - (float)fabs(j) )/float(brush_width);
    for( int i=-brush_width;i<=brush_width;i++ )
    {
      int ii = i + brush_width;
      float ifactor =  (float(brush_width) - (float)fabs(i) )/float(brush_width);
      float radius = (float) ((jfactor * jfactor + ifactor * ifactor) / 2.0);
      source_brush[ii][jj] = pow(radius,0.5);
      obstruction_brush[ii][jj] = (float)(1.0 - pow(radius, 1.0/4.0));
    }
  }
}


void setNbCores( int nb )
{
  omp_set_num_threads( nb );
}


//----------------------------------------------------

void ConvertToDisplay()
{
  for( int j=0;j<iheight;j++ )
  {
#pragma omp parallel for
    for(int i=0;i<iwidth;i++ )
    {
      int index = i + iwidth*j;
      float r,g,b;
      r = baseimage[index*3];
      g = baseimage[index*3+1];
      b = baseimage[index*3+2];
      display_map[3*index  ] = r * scaling_factor;
      display_map[3*index+1] = g * scaling_factor;
      display_map[3*index+2] = b * scaling_factor;
    }
  }
}


//------------------------------------------
//
//  Painting and display code
//

void resetScaleFactor( float amount )
{
   scaling_factor *= amount;
}



void DabSomePaint( int x, int y )
{
  int brush_width = (BRUSH_SIZE-1)/2;
  int xstart = x - brush_width;
  int ystart = y - brush_width;
  if( xstart < 0 ){ xstart = 0; }
  if( ystart < 0 ){ ystart = 0; }

  int xend = x + brush_width;
  int yend = y + brush_width;
  if( xend >= iwidth ){ xend = iwidth-1; }
  if( yend >= iheight ){ yend = iheight-1; }


  if( paint_mode == PAINT_OBSTRUCTION )
  {
    for(int ix=xstart;ix <= xend; ix++)
    {
      for( int iy=ystart;iy<=yend; iy++)
      {
        int index = ix + iwidth*(iheight-iy-1);
        baseimage[3*index] *= obstruction_brush[ix-xstart][iy-ystart];
        baseimage[3*index+1] *= obstruction_brush[ix-xstart][iy-ystart];
        baseimage[3*index+2] *= obstruction_brush[ix-xstart][iy-ystart];
      }
    }
  }
  else if( paint_mode == PAINT_SOURCE )
  {
    for(int ix=xstart;ix <= xend; ix++)
    {
      for( int iy=ystart;iy<=yend; iy++)
      {
        int index = ix + iwidth*(iheight-iy-1);
        baseimage[3*index] += source_brush[ix-xstart][iy-ystart];
        baseimage[3*index+1] += source_brush[ix-xstart][iy-ystart];
        baseimage[3*index+2] += source_brush[ix-xstart][iy-ystart];
      }
    }
  }

  return;
}


//----------------------------------------------------
//
//  GL and GLUT callbacks
//
//----------------------------------------------------



void cbDisplay( void )
{
  glClear(GL_COLOR_BUFFER_BIT );
  glDrawPixels( iwidth, iheight, GL_RGB, GL_FLOAT, display_map );
  glutSwapBuffers();
}

// animate and display new result
void cbIdle()
{
  ConvertToDisplay();
  glutPostRedisplay(); 
}


void cbOnKeyboard( unsigned char key, int x, int y )
{
  switch (key) 
  {
    case '-': case '_':
    resetScaleFactor( 0.9 );
    break;

    case '+': case '=':
    resetScaleFactor( 1.0/0.9 );
    break;

    case 'r':
    scaling_factor = 1.0;
    break;

    case ' ':
    toggle_animation_on_off = !toggle_animation_on_off;

    case 'o':
    paint_mode = PAINT_OBSTRUCTION;
    break;

    case 's':
    paint_mode = PAINT_SOURCE;
    break;

    default:
    break;
  }
}

void cbMouseDown( int button, int state, int x, int y )
{
  if( button != GLUT_LEFT_BUTTON ) { return; }
  if( state != GLUT_DOWN ) { return; }
  xmouse_prev = x;
  ymouse_prev = y;
  DabSomePaint( x, y );
}



void cbMouseMove( int x, int y )
{
  xmouse_prev = x;
  ymouse_prev = y;
  DabSomePaint( x, y ); 
}


void PrintUsage()
{
  cout << "cfd_paint keyboard choices\n";
  cout << "s       turns on painting source strength\n";
  cout << "o       turns on painting obstructions\n";
  cout << "+/-     increase/decrease brightness of display\n";
  cout << "r       resets brightness to default\n";
}


//---------------------------------------------------

int main(int argc, char** argv)
{
  CmdLineFind clf(argc, argv);

  iwidth = clf.find("-NX", 512, "Horizontal grid points");
  iheight = clf.find("-NY", iwidth, "Vertical grid points");

  setNbCores(4);

  string imagename = clf.find("-image", "grumpy.jpg", "Image to drive color");

  clf.usage("-h");
  clf.printFinds();
  PrintUsage();

  // initialize a few variables
  size = iwidth * iheight;
  scaling_factor = 1;
  toggle_animation_on_off = true;

  display_map = new float[3 * size];
  Initialize(display_map, 3 * size, 0.0);
  baseimage = new float[3 * size];
  Initialize(baseimage, 3 * size, 0.0);

  if (imagename != "") {
    readOIIOImage(imagename.c_str(), baseimage);
  }


  InitializeBrushes();

  paint_mode = PAINT_SOURCE;

  // GLUT routines
  glutInit(&argc, argv);

  glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
  glutInitWindowSize(iwidth, iheight);

  // Open a window 
  char title[] = "cfd Demo";
  glutCreateWindow(title);

  glClearColor(1, 1, 1, 1);

  glutDisplayFunc(&cbDisplay);
  glutIdleFunc(&cbIdle);
  glutKeyboardFunc(&cbOnKeyboard);
  glutMouseFunc(&cbMouseDown);
  glutMotionFunc(&cbMouseMove);

  glutMainLoop();
  return 1;
}
