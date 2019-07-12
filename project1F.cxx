#include <iostream>
#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <sstream>
#include <iomanip>
#include <string>
#define NORMALS

using std::cerr;
using std::endl;

double ceil_round(double f)
{
    return ceil(f-0.00001);
}

double floor_round(double f)
{
    return floor(f+0.00001);
}


class Matrix
{
  public:
    double          A[4][4];

    void            TransformPoint(const double *ptIn, double *ptOut);
    static Matrix   ComposeMatrices(const Matrix &, const Matrix &);
    void            Print(ostream &o);
};

void
Matrix::Print(ostream &o)
{
    for (int i = 0 ; i < 4 ; i++)
    {
        char str[256];
        sprintf(str, "(%.7f %.7f %.7f %.7f)\n", A[i][0], A[i][1], A[i][2], A[i][3]);
        o << str;
    }
}

Matrix
Matrix::ComposeMatrices(const Matrix &M1, const Matrix &M2)
{
    Matrix rv;
    for (int i = 0 ; i < 4 ; i++)
        for (int j = 0 ; j < 4 ; j++)
        {
            rv.A[i][j] = 0;
            for (int k = 0 ; k < 4 ; k++)
                rv.A[i][j] += M1.A[i][k]*M2.A[k][j];
        }

    return rv;
}

void
Matrix::TransformPoint(const double *ptIn, double *ptOut)
{
    ptOut[0] = ptIn[0]*A[0][0]
             + ptIn[1]*A[1][0]
             + ptIn[2]*A[2][0]
             + ptIn[3]*A[3][0];
    ptOut[1] = ptIn[0]*A[0][1]
             + ptIn[1]*A[1][1]
             + ptIn[2]*A[2][1]
             + ptIn[3]*A[3][1];
    ptOut[2] = ptIn[0]*A[0][2]
             + ptIn[1]*A[1][2]
             + ptIn[2]*A[2][2]
             + ptIn[3]*A[3][2];
    ptOut[3] = ptIn[0]*A[0][3]
             + ptIn[1]*A[1][3]
             + ptIn[2]*A[2][3]
             + ptIn[3]*A[3][3];
}

class Camera
{
  public:
    double          near, far;
    double          angle;
    double          position[3];
    double          focus[3];
    double          up[3];

    Matrix          ViewTransform(void);
    Matrix          CameraTransform(void);
    Matrix          DeviceTransform(void);
};


Matrix Camera::ViewTransform(void){
	Matrix m;
	m.A[0][0] = 1/tan(angle/2), m.A[0][1] = 0, m.A[0][2] = 0, m.A[0][3] = 0;
	m.A[1][0] = 0, m.A[1][1] = 1/tan(angle/2), m.A[1][2] = 0, m.A[1][3] = 0;
	
	m.A[2][0] = 0, m.A[2][1] = 0, m.A[2][2] = (far+near)/(far-near), m.A[2][3] = -1;
	m.A[3][0] = 0, m.A[3][1] = 0, m.A[3][2] = 2*(far*near)/(far-near), m.A[3][3] = 0;

	return m;
}

Matrix
Camera::CameraTransform(void){
	double O[3], U[3], V[3], W[3];

	O[0] = position[0], O[1] = position[1], O[2] = position[2];

	W[0] = O[0] - focus[0], W[1] = O[1] - focus[1], W[2] = O[2] - focus[2];
	U[0] = up[1] * W[2] - up[2] * W[1], U[1] = W[0] * up[2] - up[0] * W[2], U[2] = up[0] * W[1] - up[1] * W[0];
	
	V[0] = W[1] * U[2] - W[2] * U[1], V[1] = U[0] * W[2] - W[0] * U[2], V[2] = W[0] * U[1] - W[1] * U[0];

	double Unorm, Vnorm, Wnorm;
	Unorm = sqrt(U[0]*U[0] + U[1]*U[1] + U[2]*U[2]);
	Vnorm = sqrt(V[0]*V[0] + V[1]*V[1] + V[2]*V[2]);
	Wnorm = sqrt(W[0]*W[0] + W[1]*W[1] + W[2]*W[2]);
	for (int i = 0; i < 3; i++){
		U[i] = U[i] / Unorm;	
		V[i] = V[i] / Vnorm;	
		W[i] = W[i] / Wnorm;	
	}

	Matrix m;
	double t[3];
	t[0] = 0 - O[0], t[1] = 0 - O[1], t[2] = 0 - O[2];

	m.A[0][0] = U[0], m.A[0][1] = V[0], m.A[0][2] = W[0], m.A[0][3] = 0;
	m.A[1][0] = U[1], m.A[1][1] = V[1], m.A[1][2] = W[1], m.A[1][3] = 0;
	m.A[2][0] = U[2], m.A[2][1] = V[2], m.A[2][2] = W[2], m.A[2][3] = 0;
	
	m.A[3][0] = U[0]*t[0] + U[1]*t[1] + U[2]*t[2], m.A[3][1] = V[0]*t[0] + V[1]*t[1] + V[2]*t[2], m.A[3][2] = W[0]*t[0] + W[1]*t[1] + W[2]*t[2], m.A[3][3] = 1;
	
	return m;
}

Matrix
Camera::DeviceTransform(void){
	Matrix m;

	m.A[0][0] = 500.0, m.A[0][1] = 0, m.A[0][2] = 0, m.A[0][3] = 0;
	m.A[1][0] = 0, m.A[1][1] = 500.0, m.A[1][2] = 0, m.A[1][3] = 0;
	m.A[2][0] = 0, m.A[2][1] = 0, m.A[2][2] = 1.0, m.A[2][3] = 0;
	m.A[3][0] = 500.0, m.A[3][1] = 500.0, m.A[3][2] = 0, m.A[3][3] = 1.0;

	return m;
}

double SineParameterize(int curFrame, int nFrames, int ramp)
{
    int nNonRamp = nFrames-2*ramp;
    double height = 1./(nNonRamp + 4*ramp/M_PI);
    if (curFrame < ramp)
    {
        double factor = 2*height*ramp/M_PI;
        double eval = cos(M_PI/2*((double)curFrame)/ramp);
        return (1.-eval)*factor;
    }
    else if (curFrame > nFrames-ramp)
    {
        int amount_left = nFrames-curFrame;
        double factor = 2*height*ramp/M_PI;
        double eval =cos(M_PI/2*((double)amount_left/ramp));
        return 1. - (1-eval)*factor;
    }
    double amount_in_quad = ((double)curFrame-ramp);
    double quad_part = amount_in_quad*height;
    double curve_part = height*(2*ramp)/M_PI;
    return quad_part+curve_part;
}

Camera
GetCamera(int frame, int nframes)
{
    double t = SineParameterize(frame, nframes, nframes/10);
    Camera c;
    c.near = 5;
    c.far = 200;
    c.angle = M_PI/6;
    c.position[0] = 40*sin(2*M_PI*t);
    c.position[1] = 40*cos(2*M_PI*t);
    c.position[2] = 40;
    c.focus[0] = 0;
    c.focus[1] = 0;
    c.focus[2] = 0;
    c.up[0] = 0;
    c.up[1] = 1;
    c.up[2] = 0;
    return c;
}

vtkImageData *
NewImage(int width, int height)
{
    vtkImageData *img = vtkImageData::New();
    img->SetDimensions(width, height, 1);
    img->AllocateScalars(VTK_UNSIGNED_CHAR, 3);

    return img;
}

void
WriteImage(vtkImageData *img, const char *filename)
{
   std::string full_filename = filename;
   full_filename += ".png";
   vtkPNGWriter *writer = vtkPNGWriter::New();
   writer->SetInputData(img);
   writer->SetFileName(full_filename.c_str());
   writer->Write();
   writer->Delete();
}

class Screen
{
  public:
      unsigned char   *buffer;
      double *depthbuffer;
      int *idbuffer;
      int width, height;
      void set_pixel(int x, int y, unsigned char colors[3], int id){
	  if (x >= 0 && y >= 0 && x < width && y < height){
          	buffer[3 * (x + (y * width))] = colors[0];
          	buffer[3 * (x + (y * width)) + 1] = colors[1];
          	buffer[3 * (x + (y * width)) + 2] = colors[2];
		idbuffer[x + (y * width)] = id;
	  }
      }
      double get_depthbuffer(int x, int y){
	  if (x >= 0 && y >= 0 && x < width && y < height){
	  	return depthbuffer[x + (y * width)];
      	  	}
	  else{
		return 0;
	  }
      }
      void set_depthbuffer(int x, int y, double depth){
	  if (x >= 0 && y >= 0 && x < width && y < height){
      	  	depthbuffer[x + (y * width)] = depth;
      		}
     }
     void print(){
	for (int i = 0; i < width; i++){
		for (int j = 0; j < height; j++){
			std::cerr << "Index: " << (j + (i*height)) << " (col=" << j << ", row=" << i << ") comes from triangle " << idbuffer[j + (i*height)] << " and has color " << (int)buffer[3 * (j + (i*height))] << ", " << (int)buffer[3 * (j + (i*height))+1] << ", " << (int)buffer[3 * (j + (i*height))+2] << ", and z-value " << depthbuffer[j+(i*height)] << endl;
		}
	}
     }
};
struct LightingParameters
{
    LightingParameters(void)
    {
         lightDir[0] = -0.6;
         lightDir[1] = 0;
         lightDir[2] = -0.8;
         Ka = 0.3;
         Kd = 0.7;
         Ks = 2.3;
         alpha = 2.5;
    };

    double lightDir[3]; // The direction of the light source
    double Ka;           // The coefficient for ambient lighting.
    double Kd;           // The coefficient for diffuse lighting.
    double Ks;           // The coefficient for specular lighting.
    double alpha;        // The exponent term for specular lighting.
};


class Triangle : public Screen	
{
  public:
      bool		isGoingUp;
      int	     id;
      double 	     world_X[3], world_Y[3], world_Z[3];
      double         X[3], Y[3], Z[3];
      double         LB[3], RB[3], P[3];
      double         LB_colors[3], RB_colors[3], P_colors[3];
      double 	     slope1, slope2, bigslope;
      double	     left_b, right_b, big_b;
      double         leftEnd, midEnd, rightEnd, rowMin, rowMax;
      double         colors[3][3];
      double 	     normals[3][3];
      double 	     shading[3];
      void printTriangle(void);
      void splitTriangle(Screen, Camera);
      void setPoints(double[3], double[3], double[3], double[3], double[3], double[3], double, double, double, Screen, Camera);
      void calculateSlope();
      void fillTriangle(Screen, Camera);
      void printPoints();
      double calculateShading(LightingParameters &, double *, double *);
};

LightingParameters lp;
double Triangle::calculateShading(LightingParameters &lb, double *viewDirection, double *normal)
{
	double spec, def, amb, total;
	double R[3], V[3], light[3];
	double R_norm, V_norm, lighting_norm, RV_dot;
	
	amb = lp.Ka;
	def = lp.lightDir[0] * normal[0] + lp.lightDir[1] * normal[1] + lp.lightDir[2] * normal[2];
	
	R[0] = (2.0 * def * normal[0]) - lp.lightDir[0];
	R[1] = (2.0 * def * normal[1]) - lp.lightDir[1];
	R[2] = (2.0 * def * normal[2]) - lp.lightDir[2];
	def = abs(def * lp.Kd);

	R_norm = sqrt(R[0]*R[0] + R[1]*R[1] + R[2]*R[2]);
	V_norm = sqrt(viewDirection[0]*viewDirection[0] + viewDirection[1]*viewDirection[1] + viewDirection[2]*viewDirection[2]);

	for (int i = 0; i < 3; i++){
		R[i] = R[i] / R_norm;
		V[i] = viewDirection[i] / V_norm;
	}
	RV_dot = R[0] * V[0] + R[1] * V[1] + R[2] * V[2];
	if (RV_dot < 0.0){spec = 0.0;} 
	else{spec = (lp.Ks * std::pow(RV_dot, lp.alpha));}
	total = amb + def + spec;
	return total;
}

void Triangle::printPoints()
{
	std::cerr << std::setprecision(16) << "\nTriangle: (" << X[0] << ", " << Y[1] << "), (" << X[1] << ", " << Y[1] << "), (" << X[2] << ", " << Y[2] << ")\n" << endl;
}

void Triangle::printTriangle()
{
	if (isGoingUp == 1){
	std::cerr <<  "Rasterizing GoingUpTriangle.\nTriangle:\n" << "LB(" <<  LB[0] << " , " <<  LB[1] << " , "  << LB[2] << "), color = (" << LB_colors[0] << ", " << LB_colors[1] << ", " << LB_colors[2] << ")\nRB(" << RB[0] << " , " << RB[1] << " , " << RB[2]  << "), colors = (" << RB_colors[0] << ", " << RB_colors[1] << ", " << RB_colors[2] << ")\nP(" << P[0] << " , " << P[1] << " , " << P[2] << "), color = (" << P_colors[0] << ", " << P_colors[1] << ", " << P_colors[2] << ")" << endl;
	} else {

	std::cerr <<  "Rasterizing GoingDownTriangle.\nTriange:\n" << "LB(" <<  LB[0] << " , " <<  LB[1] << " , "  << LB[2] << "), color = (" << LB_colors[0] << ", " << LB_colors[1] << ", " << LB_colors[2] << ")\nRB(" << RB[0] << " , " << RB[1] << " , " << RB[2]  << "), colors = (" << RB_colors[0] << ", " << RB_colors[1] << ", " <<  RB_colors[2] << ")\nP(" << P[0] << " , " << P[1] << " , " << P[2] << "), color = (" << P_colors[0] << ", " << P_colors[1] << ", " << P_colors[2] << ")" << endl;
	}
}

void Triangle::setPoints(double midPoint[3], double P1[3], double P2[3], double MP_colors[3], double P1_colors[3], double P2_colors[3], double MP_shading, double P1_shading, double P2_shading, Screen screen, Camera c)
{
	double top[3], bottom[3], base[0];
	double top_colors[3];
	double bottom_colors[3];
	double base_colors[3];

	double top_shading;
	double bottom_shading;
	double base_shading;

	double t;
	double temp_color_t;
	double temp_shading_t;

	Triangle t1, t2;
	t1.isGoingUp = 1;
	t2.isGoingUp = 0;

	t1.id = id;
	t2.id = id;

	if (P1[1] - P2[1] == 0 || P1[0] - P2[0] == 0){
		bigslope = 0;
		big_b = 0;
	} else {
		bigslope = (P1[1] - P2[1]) / (P1[0] - P2[0]); 
		big_b = P2[1] - (bigslope * P2[0]);
	}

	if (P1[1] > P2[1]){
		top[0] = P1[0], top[1] = P1[1], top[2] = P1[2];
		top_colors[0] = P1_colors[0], top_colors[1] = P1_colors[1], top_colors[2] = P1_colors[2];
		top_shading = P1_shading;	
		
		bottom[0] = P2[0], bottom[1] = P2[1], bottom[2] = P2[2];
		bottom_colors[0] = P2_colors[0], bottom_colors[1] = P2_colors[1], bottom_colors[2] = P2_colors[2];
		bottom_shading = P2_shading;	
	} else {
		top[0] = P2[0], top[1] = P2[1], top[2] = P2[2];
		top_colors[0] = P2_colors[0], top_colors[1] = P2_colors[1], top_colors[2] = P2_colors[2];
		top_shading = P2_shading;	
		
		bottom[0] = P1[0], bottom[1] = P1[1], bottom[2] = P1[2];
		bottom_colors[0] = P1_colors[0], bottom_colors[1] = P1_colors[1], bottom_colors[2] = P1_colors[2];
		bottom_shading = P1_shading;	
	}

	t1.P[0] = top[0], t1.P[1] = top[1], t1.P[2] = top[2];
	t1.P_colors[0] = top_colors[0], t1.P_colors[1] = top_colors[1], t1.P_colors[2] = top_colors[2];
	t1.shading[0] = top_shading;

	t2.P[0] = bottom[0], t2.P[1] = bottom[1], t2.P[2] = bottom[2];
	t2.P_colors[0] = bottom_colors[0], t2.P_colors[1] = bottom_colors[1], t2.P_colors[2] = bottom_colors[2];
	t2.shading[0] = bottom_shading;
	if (bottom[0] < top[0]){
		if (bigslope != 0){
			base[0] = (midPoint[1] - big_b)/bigslope;
		} else {
			base[0] = P1[0];
		}
		
		if (base[0] < midPoint[0]){
			t1.LB[0] = base[0], t1.LB[1] = midPoint[1];
			t1.RB[0] = midPoint[0], t1.RB[1] = midPoint[1], t1.RB[2] = midPoint[2];
			
			if (bottom[0] == top[0]){
				t = ((top[2] - bottom[2])/(top[1] - bottom[1])); 
				t1.LB[2] = bottom[2] + (t * (midPoint[1] - bottom[1]));	
			
				temp_color_t = ((top_colors[0] - bottom_colors[0]) / (top[1] - bottom[1]));
				t1.LB_colors[0] =  bottom_colors[0] + (temp_color_t * (midPoint[1] - bottom[1]));
				temp_color_t = ((top_colors[1] - bottom_colors[1]) / (top[1] - bottom[1]));
				t1.LB_colors[1] =  bottom_colors[1] + (temp_color_t * (midPoint[1] - bottom[1]));
				temp_color_t = ((top_colors[2] - bottom_colors[2]) / (top[1] - bottom[1]));
				t1.LB_colors[2] =  bottom_colors[2] + (temp_color_t * (midPoint[1] - bottom[1]));
				
				temp_shading_t = ((top_shading - bottom_shading) / (top[1] - bottom[1]));
				t1.shading[1] = bottom_shading + (temp_shading_t * (MP_shading - bottom_shading));
			} else {
				t = (base[0] - top[0]) / (bottom[0] - top[0]);
				t1.LB[2] = top[2] + (t * (bottom[2] - top[2]));
				t1.LB_colors[0] = top_colors[0] + (t *(bottom_colors[0] - top_colors[0]));
				t1.LB_colors[1] = top_colors[1] + (t *(bottom_colors[1] - top_colors[1]));
				t1.LB_colors[2] = top_colors[2] + (t *(bottom_colors[2] - top_colors[2]));
			
				t1.shading[1] = top_shading + (t * (bottom_shading - top_shading));
			}

			t1.RB_colors[0] = MP_colors[0], t1.RB_colors[1] = MP_colors[1], t1.RB_colors[2] = MP_colors[2];
			t1.shading[2] = MP_shading;

			t2.LB[0] = t1.LB[0], t2.LB[1] = t1.LB[1], t2.LB[2] = t1.LB[2];
			t2.RB[0] = t1.RB[0], t2.RB[1] = t1.RB[1], t2.RB[2] = t1.RB[2];
			
			t2.LB_colors[0] = t1.LB_colors[0], t2.LB_colors[1] = t1.LB_colors[1], t2.LB_colors[2] = t1.LB_colors[2];
			t2.shading[1] = t1.shading[1];

			t2.RB_colors[0] = MP_colors[0], t2.RB_colors[1] = MP_colors[1], t2.RB_colors[2] = MP_colors[2];
			t2.shading[2] = t1.shading[2];
		} else {
			t1.RB[0] = base[0], t1.RB[1] = midPoint[1];
			t1.LB[0] = midPoint[0], t1.LB[1] = midPoint[1], t1.LB[2] = midPoint[2];
			if (bottom[0] == top[0]){
				t = ((top[2] - bottom[2])/(top[1] - bottom[1])); 
				t1.RB[2] = bottom[2] + (t * (midPoint[1] - bottom[1]));	
			
				temp_color_t = ((top_colors[0] - bottom_colors[0]) / (top[1] - bottom[1]));
				t1.RB_colors[0] =  bottom_colors[0] + (temp_color_t * (midPoint[1] - bottom[1]));
				temp_color_t = ((top_colors[1] - bottom_colors[1]) / (top[1] - bottom[1]));
				t1.RB_colors[1] =  bottom_colors[1] + (temp_color_t * (midPoint[1] - bottom[1]));
				temp_color_t = ((top_colors[2] - bottom_colors[2]) / (top[1] - bottom[1]));
				t1.RB_colors[2] =  bottom_colors[2] + (temp_color_t * (midPoint[1] - bottom[1]));
					
				temp_shading_t = ((top_shading - bottom_shading) / (top[1] - bottom[1]));
				t1.shading[2] = bottom_shading + (temp_shading_t * (MP_shading - bottom_shading));
			} else {
				t = (base[0] - top[0]) / (bottom[0] - top[0]);
				t1.RB[2] = top[2] + (t * (bottom[2] - top[2]));
				t1.RB_colors[0] = top_colors[0] + (t *(bottom_colors[0] - top_colors[0]));
				t1.RB_colors[1] = top_colors[1] + (t *(bottom_colors[1] - top_colors[1]));
				t1.RB_colors[2] = top_colors[2] + (t *(bottom_colors[2] - top_colors[2]));
		
				t1.shading[2] = top_shading + (t * (bottom_shading - top_shading));
			}
			
			t1.LB_colors[0] = MP_colors[0], t1.LB_colors[1] = MP_colors[1], t1.LB_colors[2] = MP_colors[2];
			t1.shading[1] = MP_shading;

			t2.LB[0] = t1.LB[0], t2.LB[1] = t1.LB[1], t2.LB[2] = t1.LB[2];
			t2.RB[0] = t1.RB[0], t2.RB[1] = t1.RB[1], t2.RB[2] = t1.RB[2];
			
			t2.LB_colors[0] = MP_colors[0], t2.LB_colors[1] = MP_colors[1], t2.LB_colors[2] = MP_colors[2];
			t2.shading[1] = t1.shading[1];

			t2.RB_colors[0] = t1.RB_colors[0], t2.RB_colors[1] = t1.RB_colors[1], t2.RB_colors[2] = t1.RB_colors[2];
			t2.shading[2] = t1.shading[2];
		}
	} else {
		if (bigslope != 0){
			base[0] = (midPoint[1] - big_b)/bigslope;
		} else {
			base[0] = P1[0];
		}
		
		if (base[0] < midPoint[0]){
			t1.LB[0] = base[0], t1.LB[1] = midPoint[1];
			t1.RB[0] = midPoint[0], t1.RB[1] = midPoint[1], t1.RB[2] = midPoint[2];
			if (bottom[0] == top[0]){
				t = ((top[2] - bottom[2])/(top[1] - bottom[1])); 
				t1.LB[2] = bottom[2] + (t * (midPoint[1] - bottom[1]));	
			
				temp_color_t = ((top_colors[0] - bottom_colors[0]) / (top[1] - bottom[1]));
				t1.LB_colors[0] =  bottom_colors[0] + (temp_color_t * (midPoint[1] - bottom[1]));
				temp_color_t = ((top_colors[1] - bottom_colors[1]) / (top[1] - bottom[1]));
				t1.LB_colors[1] =  bottom_colors[1] + (temp_color_t * (midPoint[1] - bottom[1]));
				temp_color_t = ((top_colors[2] - bottom_colors[2]) / (top[1] - bottom[1]));
				t1.LB_colors[2] =  bottom_colors[2] + (temp_color_t * (midPoint[1] - bottom[1]));
				
				temp_shading_t = ((top_shading - bottom_shading) / (top[1] - bottom[1]));
				t1.shading[1] = bottom_shading + (temp_shading_t * (MP_shading - bottom_shading));
			} else {
				t = (base[0] - top[0]) / (bottom[0] - top[0]);
				t1.LB[2] = top[2] + (t * (bottom[2] - top[2]));
				t1.LB_colors[0] = top_colors[0] + (t *(bottom_colors[0] - top_colors[0]));
				t1.LB_colors[1] = top_colors[1] + (t *(bottom_colors[1] - top_colors[1]));
				t1.LB_colors[2] = top_colors[2] + (t *(bottom_colors[2] - top_colors[2]));
			
				t1.shading[1] = top_shading + (t * (bottom_shading - top_shading));
			}
			
			t1.RB_colors[0] = MP_colors[0], t1.RB_colors[1] = MP_colors[1], t1.RB_colors[2] = MP_colors[2];
			t1.shading[2] = MP_shading;

			t2.RB[0] = t1.RB[0], t2.RB[1] = t1.RB[1], t2.RB[2] = t1.RB[2];
			t2.LB[0] = t1.LB[0], t2.LB[1] = t1.LB[1], t2.LB[2] = t1.LB[2];
			
			t2.RB_colors[0] = MP_colors[0], t2.RB_colors[1] = MP_colors[1], t2.RB_colors[2] = MP_colors[2];
			t2.shading[2] = t1.shading[2];
			t2.LB_colors[0] = t1.LB_colors[0], t2.LB_colors[1] = t1.LB_colors[1], t2.LB_colors[2] = t1.LB_colors[2];
			t2.shading[1] = t1.shading[1];
		} else {
			t1.RB[0] = base[0], t1.RB[1] = midPoint[1];
			t1.LB[0] = midPoint[0], t1.LB[1] = midPoint[1], t1.LB[2] = midPoint[2];
			if (bottom[0] == top[0]){
				t = ((top[2] - bottom[2])/(top[1] - bottom[1])); 
				t1.RB[2] = bottom[2] + (t * (midPoint[1] - bottom[1]));	
			
				temp_color_t = ((top_colors[0] - bottom_colors[0]) / (top[1] - bottom[1]));
				t1.RB_colors[0] =  bottom_colors[0] + (temp_color_t * (midPoint[1] - bottom[1]));
				temp_color_t = ((top_colors[1] - bottom_colors[1]) / (top[1] - bottom[1]));
				t1.RB_colors[1] =  bottom_colors[1] + (temp_color_t * (midPoint[1] - bottom[1]));
				temp_color_t = ((top_colors[2] - bottom_colors[2]) / (top[1] - bottom[1]));
				t1.RB_colors[2] =  bottom_colors[2] + (temp_color_t * (midPoint[1] - bottom[1]));
				
				temp_shading_t = ((top_shading - bottom_shading) / (top[1] - bottom[1]));
				t1.shading[2] = bottom_shading + (temp_shading_t * (MP_shading - bottom_shading));
			} else {
				t = (base[0] - top[0]) / (bottom[0] - top[0]);
				t1.RB[2] = top[2] + (t * (bottom[2] - top[2]));
				t1.RB_colors[0] = top_colors[0] + (t *(bottom_colors[0] - top_colors[0]));
				t1.RB_colors[1] = top_colors[1] + (t *(bottom_colors[1] - top_colors[1]));
				t1.RB_colors[2] = top_colors[2] + (t *(bottom_colors[2] - top_colors[2]));
			
				t1.shading[2] = top_shading + (t * (bottom_shading - top_shading));
			}		
			
			t1.LB_colors[0] = MP_colors[0], t1.LB_colors[1] = MP_colors[1], t1.LB_colors[2] = MP_colors[2];
			t1.shading[1] = MP_shading;

			t2.RB[0] = t1.RB[0], t2.RB[1] = t1.RB[1], t2.RB[2] = t1.RB[2];
			t2.LB[0] = t1.LB[0], t2.LB[1] = t1.LB[1], t2.LB[2] = t1.LB[2];
			
			t2.LB_colors[0] = MP_colors[0], t2.LB_colors[1] = MP_colors[1], t2.LB_colors[2] = MP_colors[2];
			t2.RB_colors[0] = t1.RB_colors[0], t2.RB_colors[1] = t1.RB_colors[1], t2.RB_colors[2] = t1.RB_colors[2];
			t2.shading[1] = t1.shading[1];
			t2.shading[2] = t1.shading[2];
		}

	}
	t1.calculateSlope(); 
	t2.calculateSlope();
	t1.fillTriangle(screen, c);
	t2.fillTriangle(screen, c);
}

void Triangle::splitTriangle(Screen screen, Camera c)
{
	double viewDirection[3];
	if ((Y[0] > Y[1] && Y[0] < Y[2]) || (Y[0] < Y[1] && Y[0] > Y[2])){
		double MP[3], P1[3], P2[3];
		double MP_colors[3], P1_colors[3], P2_colors[3];
		double MP_shading, P1_shading, P2_shading;
		
		viewDirection[0] = c.position[0] - world_X[0], viewDirection[1] = c.position[1] - world_Y[0], viewDirection[2] = c.position[2] - world_Z[0];
		MP[0] = X[0], MP[1] = Y[0], MP[2] = Z[0];
		MP_colors[0] = colors[0][0], MP_colors[1] = colors[0][1], MP_colors[2] = colors[0][2];
		MP_shading = calculateShading(lp, viewDirection, normals[0]);

		viewDirection[0] = c.position[0] - world_X[1], viewDirection[1] = c.position[1] - world_Y[1], viewDirection[2] = c.position[2] - world_Z[1];
		P1[0] = X[1], P1[1] = Y[1], P1[2] = Z[1];
		P1_colors[0] = colors[1][0], P1_colors[1] = colors[1][1], P1_colors[2] = colors[1][2];
		P1_shading = calculateShading(lp, viewDirection, normals[1]);
		
		viewDirection[0] = c.position[0] - world_X[2], viewDirection[1] = c.position[1] - world_Y[2], viewDirection[2] = c.position[2] - world_Z[2];
		P2[0] = X[2], P2[1] = Y[2], P2[2] = Z[2];
		P2_colors[0] = colors[2][0], P2_colors[1] = colors[2][1], P2_colors[2] = colors[2][2];
		P2_shading = calculateShading(lp, viewDirection, normals[2]);

		setPoints(MP, P1, P2, MP_colors, P1_colors, P2_colors, MP_shading, P1_shading, P2_shading, screen, c);
	}
	else if ((Y[1] > Y[0] && Y[1] < Y[2]) || (Y[1] < Y[0] && Y[1] > Y[2])){
		double MP[3], P1[3], P2[3];
		double MP_colors[3], P1_colors[3], P2_colors[3];
		double MP_shading, P1_shading, P2_shading;
		double viewDirection[3];
		
		viewDirection[0] = c.position[0] - world_X[1], viewDirection[1] = c.position[1] - world_Y[1], viewDirection[2] = c.position[2] - world_Z[1];
		MP[0] = X[1], MP[1] = Y[1], MP[2] = Z[1];
		MP_colors[0] = colors[1][0], MP_colors[1] = colors[1][1], MP_colors[2] = colors[1][2];
		MP_shading = calculateShading(lp, viewDirection, normals[1]);
		
		viewDirection[0] = c.position[0] - world_X[0], viewDirection[1] = c.position[1] - world_Y[0], viewDirection[2] = c.position[2] - world_Z[0];
		P1[0] = X[0], P1[1] = Y[0], P1[2] = Z[0];
		P1_colors[0] = colors[0][0], P1_colors[1] = colors[0][1], P1_colors[2] = colors[0][2];
		P1_shading = calculateShading(lp, viewDirection, normals[0]);
		
		viewDirection[0] = c.position[0] - world_X[2], viewDirection[1] = c.position[1] - world_Y[2], viewDirection[2] = c.position[2] - world_Z[2];
		P2[0] = X[2], P2[1] = Y[2], P2[2] = Z[2];
		P2_colors[0] = colors[2][0], P2_colors[1] = colors[2][1], P2_colors[2] = colors[2][2];
		P2_shading = calculateShading(lp, viewDirection, normals[2]);


		setPoints(MP, P1, P2, MP_colors, P1_colors, P2_colors, MP_shading, P1_shading, P2_shading, screen, c);
	}
	else if ((Y[2] > Y[1] && Y[2] < Y[0]) || (Y[2] < Y[1] && Y[2] > Y[0])){
		double MP[3], P1[3], P2[3];
		double MP_colors[3], P1_colors[3], P2_colors[3];
		double MP_shading, P1_shading, P2_shading;
		double viewDirection[3];
		
		viewDirection[0] = c.position[0] - world_X[2], viewDirection[1] = c.position[1] - world_Y[2], viewDirection[2] = c.position[2] - world_Z[2];
		MP[0] = X[2], MP[1] = Y[2], MP[2] = Z[2];
		MP_colors[0] = colors[2][0], MP_colors[1] = colors[2][1], MP_colors[2] = colors[2][2];
		MP_shading = calculateShading(lp, viewDirection, normals[2]);

		viewDirection[0] = c.position[0] - world_X[1], viewDirection[1] = c.position[1] - world_Y[1], viewDirection[2] = c.position[2] - world_Z[1];
		P1[0] = X[1], P1[1] = Y[1], P1[2] = Z[1];
		P1_colors[0] = colors[1][0], P1_colors[1] = colors[1][1], P1_colors[2] = colors[1][2];
		P1_shading = calculateShading(lp, viewDirection, normals[1]);
		
		viewDirection[0] = c.position[0] - world_X[0], viewDirection[1] = c.position[1] - world_Y[0], viewDirection[2] = c.position[2] - world_Z[0];
		P2[0] = X[0], P2[1] = Y[0], P2[2] = Z[0];
		P2_colors[0] = colors[0][0], P2_colors[1] = colors[0][1], P2_colors[2] = colors[0][2];
		P2_shading = calculateShading(lp, viewDirection, normals[0]);

		setPoints(MP, P1, P2, MP_colors, P1_colors, P2_colors, MP_shading, P1_shading, P2_shading, screen, c);
	} else {
		if (Y[0] == Y[1]){
			P[0] = X[2], P[1] = Y[2], P[2] = Z[2];
			P_colors[0] = colors[2][0], P_colors[1] = colors[2][1], P_colors[2] = colors[2][2];
			viewDirection[0] = c.position[0] - world_X[2], viewDirection[1] = c.position[1] - world_Y[2], viewDirection[2] = c.position[2] - world_Z[2];
			shading[0] = calculateShading(lp, viewDirection, normals[2]);
			if (X[0] < X[1]){
				LB[0] = X[0], LB[1] = Y[0], LB[2] = Z[0];
				LB_colors[0] = colors[0][0], LB_colors[1] = colors[0][1], LB_colors[2] = colors[0][2];
				shading[1] = calculateShading(lp, viewDirection, normals[0]);
				RB[0] = X[1], RB[1] = Y[1], RB[2] = Z[1];
				RB_colors[0] = colors[1][0], RB_colors[1] = colors[1][1], RB_colors[2] = colors[1][2];
				shading[2] = calculateShading(lp, viewDirection, normals[1]);
			} else {
				LB[0] = X[1], LB[1] = Y[1], LB[2] = Z[1];
				LB_colors[0] = colors[1][0], LB_colors[1] = colors[1][1], LB_colors[2] = colors[1][2];
				shading[1] = calculateShading(lp, viewDirection, normals[1]);
					
				RB[0] = X[0], RB[1] = Y[0], RB[2] = Z[0];
				RB_colors[0] = colors[0][0], RB_colors[1] = colors[0][1], RB_colors[2] = colors[0][2];
				shading[2] = calculateShading(lp, viewDirection, normals[0]);
			}
			if (P[1] < LB[1]){
				isGoingUp = 0;
			} else {
				isGoingUp = 1;
			}
		} else if (Y[0] == Y[2]){
			P[0] = X[1], P[1] = Y[1], P[2] = Z[1];
			P_colors[0] = colors[1][0], P_colors[1] = colors[1][1], P_colors[2] = colors[1][2];
			viewDirection[0] = c.position[0] - world_X[1], viewDirection[1] = c.position[1] - world_Y[1], viewDirection[2] = c.position[2] - world_Z[1];
			shading[0] = calculateShading(lp, viewDirection, normals[1]);
			if (X[0] < X[2]){
				LB[0] = X[0], LB[1] = Y[0], LB[2] = Z[0];
				LB_colors[0] = colors[0][0], LB_colors[1] = colors[0][1], LB_colors[2] = colors[0][2];
				shading[1] = calculateShading(lp, viewDirection, normals[0]);
				
				RB[0] = X[2], RB[1] = Y[2], RB[2] = Z[2];
				RB_colors[0] = colors[2][0], RB_colors[1] = colors[2][1], RB_colors[2] = colors[2][2];
				shading[2] = calculateShading(lp, viewDirection, normals[2]);
			} else {
				LB[0] = X[2], LB[1] = Y[2], LB[2] = Z[2];
				LB_colors[0] = colors[2][0], LB_colors[1] = colors[2][1], LB_colors[2] = colors[2][2];
				shading[1] = calculateShading(lp, viewDirection, normals[2]);
				
				RB[0] = X[0], RB[1] = Y[0], RB[2] = Z[0];
				RB_colors[0] = colors[0][0], RB_colors[1] = colors[0][1], RB_colors[2] = colors[0][2];
				shading[2] = calculateShading(lp, viewDirection, normals[0]);
			}
			if (P[1] < LB[1]){
				isGoingUp = 0;
			} else {
				isGoingUp = 1;
			}
		} else if (Y[1] == Y[2]){
			P[0] = X[0], P[1] = Y[0], P[2] = Z[0];
			P_colors[0] = colors[0][0], P_colors[1] = colors[0][1], P_colors[2] = colors[0][2];
			viewDirection[0] = c.position[0] - world_X[0], viewDirection[1] = c.position[1] - world_Y[0], viewDirection[2] = c.position[2] - world_Z[0];
			shading[0] = calculateShading(lp, viewDirection, normals[0]);
			if (X[1] < X[2]){
				LB[0] = X[1], LB[1] = Y[1], LB[2] = Z[1];
				LB_colors[0] = colors[1][0], LB_colors[1] = colors[1][1], LB_colors[2] = colors[1][2];
				shading[1] = calculateShading(lp, viewDirection, normals[1]);
				
				RB[0] = X[2], RB[1] = Y[2], RB[2] = Z[2];
				RB_colors[0] = colors[2][0], RB_colors[1] = colors[2][1], RB_colors[2] = colors[2][2];
				shading[2] = calculateShading(lp, viewDirection, normals[2]);
			} else {
				LB[0] = X[2], LB[1] = Y[2], LB[2] = Z[2];
				LB_colors[0] = colors[2][0], LB_colors[1] = colors[2][1], LB_colors[2] = colors[2][2];
				shading[1] = calculateShading(lp, viewDirection, normals[2]);
				
				RB[0] = X[1], RB[1] = Y[1], RB[2] = Z[1];
				RB_colors[0] = colors[1][0], RB_colors[1] = colors[1][1], RB_colors[2] = colors[1][2];
				shading[2] = calculateShading(lp, viewDirection, normals[1]);
			}
			if (P[1] < LB[1]){
				isGoingUp = 0;
			} else {
				isGoingUp = 1;
			}
		}
		calculateSlope();
		fillTriangle(screen, c);
	}
}

void Triangle::calculateSlope()
{
	if (LB[0] - P[0] == 0 || LB[1] - P[1] == 0){
		slope1 = 0;
		left_b = 0;
	} else {
		slope1 = (LB[1] - P[1]) / (LB[0] - P[0]);
		left_b = LB[1] - (slope1 * LB[0]);
	}
	
	if (RB[0] - P[0] == 0 || RB[1] - P[1] == 0){
		slope2 = 0;
		right_b = 0;
	} else { 
		slope2 = (RB[1] - P[1]) / (RB[0] - P[0]);
		right_b = RB[1] - (slope2 * RB[0]);
	}
}

void Triangle::fillTriangle(Screen screen, Camera c)
{	
	if (isGoingUp == 1)
	{
		rowMin = ceil_round(LB[1]), rowMax = floor_round(P[1]);
		double t;
		double left_t, right_t;
		double leftZ, rightZ;
		double left_color[3], right_color[3];
		double left_shading, right_shading;
		double temp_color_t;
		double temp_shading_t;
		for (float r = rowMin; r <= rowMax; r++){
			if (slope1 != 0){
				leftEnd = (r - left_b)/slope1;
			} else {
				leftEnd = LB[0];
			}
			if (slope2 != 0){
				rightEnd = (r - right_b)/slope2;
			} else {
				rightEnd = RB[0];
			}
				
			if (LB[0] - P[0] == 0){
				left_t = ((LB[2] - P[2])/(LB[1] - P[1])); 
				leftZ =  P[2]+(left_t * (r - P[1]));
			
				temp_color_t = ((LB_colors[0] - P_colors[0]) / (LB[1] - P[1]));
				left_color[0] = P_colors[0] + (temp_color_t * (r - P[1]));
				temp_color_t = ((LB_colors[1] - P_colors[1]) / (LB[1] - P[1]));
				left_color[1] = P_colors[1] + (temp_color_t * (r - P[1]));
				temp_color_t = ((LB_colors[2] - P_colors[2]) / (LB[1] - P[1]));
				left_color[2] = P_colors[2] + (temp_color_t * (r - P[1]));
		
				temp_shading_t = ((shading[1] - shading[0]) / (LB[1] - P[1]));
				left_shading = shading[0] + (temp_shading_t * (r - P[1]));
			} else{
				left_t = (leftEnd - P[0]) / (LB[0] - P[0]);
				leftZ = P[2] + (left_t * (LB[2] - P[2]));
		
				left_color[0] = P_colors[0] + (left_t *(LB_colors[0] - P_colors[0]));
				left_color[1] = P_colors[1] + (left_t *(LB_colors[1] - P_colors[1]));
				left_color[2] = P_colors[2] + (left_t *(LB_colors[2] - P_colors[2]));
			
				left_shading = shading[0] + (left_t * (shading[1] - shading[0]));
			}
			
			if (P[0] - RB[0] == 0){
				
				right_t = ((RB[2] - P[2])/(RB[1] - P[1]));
				rightZ =  P[2]+(right_t * (r - P[1]));

				temp_color_t = ((RB_colors[0] - P_colors[0]) / (RB[1] - P[1]));
				right_color[0] = P_colors[0] + (temp_color_t * (r - P[1]));
				temp_color_t = ((RB_colors[1] - P_colors[1]) / (RB[1] - P[1]));
				right_color[1] = P_colors[1] + (temp_color_t * (r - P[1]));
				temp_color_t = ((RB_colors[2] - P_colors[2]) / (RB[1] - P[1]));
				right_color[2] = P_colors[2] + (temp_color_t * (r - P[1]));
				
				temp_shading_t = ((shading[2] - shading[0]) / (RB[1] - P[1]));
				right_shading = shading[0] + (temp_shading_t * (r - P[1]));
			} else{
				right_t = (rightEnd - RB[0]) / (P[0] - RB[0]);
				rightZ = RB[2] + (right_t * (P[2] - RB[2]));
				right_color[0] = RB_colors[0] + (right_t *(P_colors[0] - RB_colors[0]));
				right_color[1] = RB_colors[1] + (right_t *(P_colors[1] - RB_colors[1]));
				right_color[2] = RB_colors[2] + (right_t *(P_colors[2] - RB_colors[2]));
				
				right_shading = shading[2] + (right_t * (shading[0] - shading[2]));
			}		

			for (int c = ceil_round(leftEnd); c <= floor_round(rightEnd); c++){
				double color_d[3];
				unsigned char color_uc[3];
				double currentZ;
				double current_t;
				double current_shading;
				double color_t;
				double shading_t;

				current_t = ((leftZ - rightZ)/(leftEnd - rightEnd));
				currentZ = rightZ + (current_t * (c - rightEnd));
			
				color_t = ((left_color[0] - right_color[0]) / (leftEnd - rightEnd));  
				color_d[0] = right_color[0] + (color_t * (c - rightEnd));
				
				color_t = ((left_color[1] - right_color[1]) / (leftEnd - rightEnd));  
				color_d[1] = right_color[1] + (color_t * (c - rightEnd));
				
				color_t = ((left_color[2] - right_color[2]) / (leftEnd - rightEnd));  
				color_d[2] = right_color[2] + (color_t * (c - rightEnd));

				shading_t = ((left_shading - right_shading) / (leftEnd - rightEnd));
				current_shading = right_shading + (shading_t * (c - rightEnd));

				color_d[0] = color_d[0] * current_shading;
				color_d[1] = color_d[1] * current_shading;
				color_d[2] = color_d[2] * current_shading;

				if (color_d[0] > 1.0){ color_d[0] = 1.0;}
				if (color_d[1] > 1.0){ color_d[1] = 1.0;}
				if (color_d[2] > 1.0){ color_d[2] = 1.0;}
				
				color_uc[0] = ceil_round(color_d[0] * 255.0);
				color_uc[1] = ceil_round(color_d[1] * 255.0);
				color_uc[2] = ceil_round(color_d[2] * 255.0);
				
				if (currentZ > screen.get_depthbuffer(c, r)){
					screen.set_pixel(c, r, color_uc, id);
					screen.set_depthbuffer(c, r, currentZ);
				}
			}
		}
	} else
	{
		int xValue, yValue;
		xValue = 0, yValue = 1;
		
		double left_t, right_t;
		double leftZ, rightZ;
		double left_color[3], right_color[3];
		double left_shading, right_shading;
		double temp_color_t;
		double temp_shading_t;
		rowMin = ceil_round(P[yValue]), rowMax = floor_round(LB[yValue]);
		
		for (float r = rowMax; r >= rowMin; r--){
			if (slope1 != 0){
				leftEnd = (r - left_b)/slope1;
			} else {
				leftEnd = LB[0];
			}
			if (slope2 != 0){
				rightEnd = (r - right_b)/slope2;
			} else {
				rightEnd = RB[0];
			}
			if (LB[0] - P[0] == 0){
				left_t = ((LB[2] - P[2])/(LB[1] - P[1])); 
				leftZ =  P[2]+(left_t * (r - P[1]));

				temp_color_t = ((LB_colors[0] - P_colors[0]) / (LB[1] - P[1]));
				left_color[0] = P_colors[0] + (temp_color_t * (r - P[1]));
				temp_color_t = ((LB_colors[1] - P_colors[1]) / (LB[1] - P[1]));
				left_color[1] = P_colors[1] + (temp_color_t * (r - P[1]));
				temp_color_t = ((LB_colors[2] - P_colors[2]) / (LB[1] - P[1]));
				left_color[2] = P_colors[2] + (temp_color_t * (r - P[1]));
				
				temp_shading_t = ((shading[1] - shading[0]) / (LB[1] - P[1]));
				left_shading = shading[0] + (temp_shading_t * (r - P[1]));
			} else{
				left_t = (leftEnd - P[0]) / (LB[0] - P[0]);
				leftZ = P[2] + (left_t * (LB[2] - P[2]));
				
				left_color[0] = P_colors[0] + (left_t *(LB_colors[0] - P_colors[0]));
				left_color[1] = P_colors[1] + (left_t *(LB_colors[1] - P_colors[1]));
				left_color[2] = P_colors[2] + (left_t *(LB_colors[2] - P_colors[2]));
				
				left_shading = shading[0] + (left_t * (shading[1] - shading[0]));
			}
			
			if (P[0] - RB[0] == 0){
				right_t = ((RB[2] - P[2])/(RB[1] - P[1]));
				rightZ =  P[2]+(right_t * (r - P[1]));

				temp_color_t = ((RB_colors[0] - P_colors[0]) / (RB[1] - P[1]));
				right_color[0] = P_colors[0] + (temp_color_t * (r - P[1]));
				temp_color_t = ((RB_colors[1] - P_colors[1]) / (RB[1] - P[1]));
				right_color[1] = P_colors[1] + (temp_color_t * (r - P[1]));
				temp_color_t = ((RB_colors[2] - P_colors[2]) / (RB[1] - P[1]));
				right_color[2] = P_colors[2] + (temp_color_t * (r - P[1]));
				
				temp_shading_t = ((shading[2] - shading[0]) / (RB[1] - P[1]));
				right_shading = shading[0] + (temp_shading_t * (r - P[1]));
			} else{
				right_t = (rightEnd - RB[0]) / (P[0] - RB[0]);
				rightZ = RB[2] + (right_t * (P[2] - RB[2]));
				
				right_color[0] = RB_colors[0] + (right_t *(P_colors[0] - RB_colors[0]));
				right_color[1] = RB_colors[1] + (right_t *(P_colors[1] - RB_colors[1]));
				right_color[2] = RB_colors[2] + (right_t *(P_colors[2] - RB_colors[2]));
				
				right_shading = shading[2] + (right_t * (shading[0] - shading[2]));
			}		
						
			for (int c = ceil_round(leftEnd); c <= floor_round(rightEnd); c++){
				
				double color_d[3];
				unsigned char color_uc[3];
				double currentZ;
				double current_t;
				double current_shading;
				double color_t;
				double shading_t;
				
				current_t = ((rightZ - leftZ)/(rightEnd - leftEnd));
				currentZ = leftZ + (current_t * (c - leftEnd));
				
				color_t = ((right_color[0] - left_color[0]) / (rightEnd - leftEnd));  
				color_d[0] = left_color[0] + (color_t * (c - leftEnd));
				
				color_t = ((right_color[1] - left_color[1]) / (rightEnd - leftEnd));  
				color_d[1] = left_color[1] + (color_t * (c - leftEnd));
				
				color_t = ((right_color[2] - left_color[2]) / (rightEnd - leftEnd));  
				color_d[2] = left_color[2] + (color_t * (c - leftEnd));
				
				shading_t = ((right_shading - left_shading) / (rightEnd - leftEnd));
				current_shading = left_shading + (shading_t * (c - leftEnd));

				color_d[0] = color_d[0] * current_shading;
				color_d[1] = color_d[1] * current_shading;
				color_d[2] = color_d[2] * current_shading;

				if (color_d[0] > 1.0){ color_d[0] = 1.0;}
				if (color_d[1] > 1.0){ color_d[1] = 1.0;}
				if (color_d[2] > 1.0){ color_d[2] = 1.0;}
				
				color_uc[0] = ceil_round(color_d[0] * 255.0);
				color_uc[1] = ceil_round(color_d[1] * 255.0);
				color_uc[2] = ceil_round(color_d[2] * 255.0);

				if (currentZ > screen.get_depthbuffer(c, r)){
					screen.set_pixel(c, r, color_uc, id);
					screen.set_depthbuffer(c, r, currentZ);
				}
			}
		}
	}
};

std::vector<Triangle>
GetTriangles(void)
{
    vtkPolyDataReader *rdr = vtkPolyDataReader::New();
    rdr->SetFileName("proj1f_geometry.vtk");
    cerr << "Reading" << endl;
    rdr->Update();
    cerr << "Done reading" << endl;
    if (rdr->GetOutput()->GetNumberOfCells() == 0)
    {
        cerr << "Unable to open file!!" << endl;
        exit(EXIT_FAILURE);
    }
    vtkPolyData *pd = rdr->GetOutput();

    int numTris = pd->GetNumberOfCells();
    vtkPoints *pts = pd->GetPoints();
    vtkCellArray *cells = pd->GetPolys();
    vtkDoubleArray *var = (vtkDoubleArray *) pd->GetPointData()->GetArray("hardyglobal");
    double *color_ptr = var->GetPointer(0);
    //vtkFloatArray *var = (vtkFloatArray *) pd->GetPointData()->GetArray("hardyglobal");
    //float *color_ptr = var->GetPointer(0);
    vtkFloatArray *n = (vtkFloatArray *) pd->GetPointData()->GetNormals();
    float *normals = n->GetPointer(0);
    std::vector<Triangle> tris(numTris);
    vtkIdType npts;
    vtkIdType *ptIds;
    int idx;
    for (idx = 0, cells->InitTraversal() ; cells->GetNextCell(npts, ptIds) ; idx++)
    {
        if (npts != 3)
        {
            cerr << "Non-triangles!! ???" << endl;
            exit(EXIT_FAILURE);
        }
        double *pt = NULL;
        pt = pts->GetPoint(ptIds[0]);
        tris[idx].X[0] = pt[0];
        tris[idx].Y[0] = pt[1];
        tris[idx].Z[0] = pt[2];
#ifdef NORMALS
        tris[idx].normals[0][0] = normals[3*ptIds[0]+0];
        tris[idx].normals[0][1] = normals[3*ptIds[0]+1];
        tris[idx].normals[0][2] = normals[3*ptIds[0]+2];
#endif
        pt = pts->GetPoint(ptIds[1]);
        tris[idx].X[1] = pt[0];
        tris[idx].Y[1] = pt[1];
        tris[idx].Z[1] = pt[2];
#ifdef NORMALS
        tris[idx].normals[1][0] = normals[3*ptIds[1]+0];
        tris[idx].normals[1][1] = normals[3*ptIds[1]+1];
        tris[idx].normals[1][2] = normals[3*ptIds[1]+2];
#endif
        pt = pts->GetPoint(ptIds[2]);
        tris[idx].X[2] = pt[0];
        tris[idx].Y[2] = pt[1];
        tris[idx].Z[2] = pt[2];
#ifdef NORMALS
        tris[idx].normals[2][0] = normals[3*ptIds[2]+0];
        tris[idx].normals[2][1] = normals[3*ptIds[2]+1];
        tris[idx].normals[2][2] = normals[3*ptIds[2]+2];
#endif

        // 1->2 interpolate between light blue, dark blue
        // 2->2.5 interpolate between dark blue, cyan
        // 2.5->3 interpolate between cyan, green
        // 3->3.5 interpolate between green, yellow
        // 3.5->4 interpolate between yellow, orange
        // 4->5 interpolate between orange, brick
        // 5->6 interpolate between brick, salmon
        double mins[7] = { 1, 2, 2.5, 3, 3.5, 4, 5 };
        double maxs[7] = { 2, 2.5, 3, 3.5, 4, 5, 6 };
        unsigned char RGB[8][3] = { { 71, 71, 219 }, 
                                    { 0, 0, 91 },
                                    { 0, 255, 255 },
                                    { 0, 128, 0 },
                                    { 255, 255, 0 },
                                    { 255, 96, 0 },
                                    { 107, 0, 0 },
                                    { 224, 76, 76 } 
                                  };
        for (int j = 0 ; j < 3 ; j++)
        {
            float val = color_ptr[ptIds[j]];
            int r;
            for (r = 0 ; r < 7 ; r++)
            {
                if (mins[r] <= val && val < maxs[r])
                    break;
            }
            if (r == 7)
            {
                cerr << "Could not interpolate color for " << val << endl;
                exit(EXIT_FAILURE);
            }
            double proportion = (val-mins[r]) / (maxs[r]-mins[r]);
            tris[idx].colors[j][0] = (RGB[r][0]+proportion*(RGB[r+1][0]-RGB[r][0]))/255.0;
            tris[idx].colors[j][1] = (RGB[r][1]+proportion*(RGB[r+1][1]-RGB[r][1]))/255.0;
            tris[idx].colors[j][2] = (RGB[r][2]+proportion*(RGB[r+1][2]-RGB[r][2]))/255.0;
        }
    }

    return tris;
}

int main()
{
	vtkImageData *image = NewImage(1000, 1000);
	unsigned char *buffer = 
	  (unsigned char *) image->GetScalarPointer(0,0,0);
	double depthbuffer[1000*1000];
	int *idbuffer = new int[1000*1000];
	int npixels = 1000*1000;
	for (int i = 0 ; i < npixels*3 ; i++){
		 buffer[i] = 0;
	}
	for (int i = 0; i < npixels; i++){
		depthbuffer[i] = -1.0;
	}
	for (int i = 0; i < npixels; i++){
		idbuffer[i] = -1;
	}
	
	std::vector<Triangle> triangles = GetTriangles();
	
	Screen screen;
	screen.buffer = buffer;
	screen.depthbuffer = depthbuffer;
	screen.idbuffer = idbuffer;
	screen.width = 1000;
	screen.height = 1000;
	for (int i = 0; i < 1; i++){
		for (int i = 0 ; i < npixels*3 ; i++){
			 buffer[i] = 0;
		}
		for (int i = 0; i < npixels; i++){
			depthbuffer[i] = -1.0;
		}
		for (int i = 0; i < npixels; i++){
			idbuffer[i] = -1;
		}
		screen.buffer = buffer;
		screen.depthbuffer = depthbuffer;
		screen.idbuffer = idbuffer;

		Matrix comp, total, camera_transform, view_transform, device_transform;

		Camera c = GetCamera(i, 1000);
		double in[4], out[4];
		
		camera_transform = c.CameraTransform(), view_transform = c.ViewTransform(), device_transform = c.DeviceTransform();
		
		comp = comp.ComposeMatrices(camera_transform, view_transform);
		total = total.ComposeMatrices(comp, device_transform);

		for (int i = 0; i < triangles.size(); i++){
			Triangle t_mod;
			Triangle t = triangles[i];
			
			t_mod.world_X[0] = t.X[0], t_mod.world_X[1] = t.X[1], t_mod.world_X[2] = t.X[2]; 
			t_mod.world_Y[0] = t.Y[0], t_mod.world_Y[1] = t.Y[1], t_mod.world_Y[2] = t.Y[2]; 
			t_mod.world_Z[0] = t.Z[0], t_mod.world_Z[1] = t.Z[1], t_mod.world_Z[2] = t.Z[2]; 

			for (int j = 0; j <= 2; j++){
				in[0] = t.X[j];
				in[1] = t.Y[j];
				in[2] = t.Z[j];
				in[3] = 1;

				total.TransformPoint(in, out);
				
				out[0] = out[0] / out[3];
				out[1] = out[1] / out[3];
				out[2] = out[2] / out[3];
				t_mod.X[j] = out[0];
				t_mod.Y[j] = out[1];
				t_mod.Z[j] = out[2];
			}

			t_mod.colors[0][0] = t.colors[0][0], t_mod.colors[0][1] = t.colors[0][1], t_mod.colors[0][2] = t.colors[0][2];
			t_mod.colors[1][0] = t.colors[1][0], t_mod.colors[1][1] = t.colors[1][1], t_mod.colors[1][2] = t.colors[1][2];
			t_mod.colors[2][0] = t.colors[2][0], t_mod.colors[2][1] = t.colors[2][1], t_mod.colors[2][2] = t.colors[2][2];
			t_mod.normals[0][0] = t.normals[0][0], t_mod.normals[0][1] = t.normals[0][1], t_mod.normals[0][2] = t.normals[0][2];
			t_mod.normals[1][0] = t.normals[1][0], t_mod.normals[1][1] = t.normals[1][1], t_mod.normals[1][2] = t.normals[1][2];
			t_mod.normals[2][0] = t.normals[2][0], t_mod.normals[2][1] = t.normals[2][1], t_mod.normals[2][2] = t.normals[2][2];
			t_mod.id = i;
			t_mod.splitTriangle(screen, c);
		}
		
		std::stringstream ss;
		ss << "frame";
		ss << std::setw(3) << std::setfill('0') << i;
		const std::string s = ss.str();
		const char* tmp = s.c_str();
		WriteImage(image, tmp);
	}
}
