// cameras/realistic.cpp*
#include "stdafx.h"
#include "cameras/realistic.h"
#include "paramset.h"
#include "sampler.h"
#include "montecarlo.h"
#include "filters/box.h"
#include "film/image.h"
#include "samplers/stratified.h"
#include "intersection.h"
#include "renderer.h"

#include <fstream>
#include <algorithm>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>

#define _USE_MATH_DEFINES

using namespace std;


void DebugPoint(const Point& P) {
	cout << "(" << P.x << ", " << P.y << ", " << P.z << ")" << endl;
}

void DebugVector(const Vector& V) {
	cout << "<" << V.x << ", " << V.y << ", " << V.z << ">" << endl;
}

void DebugRay(const Ray& R) {
	cout << "(" << R.o.x << ", " << R.o.y << ", " << R.o.z << ") -> <" << R.d.x << ", " << R.d.y << ", " << R.d.z << ">" << endl;
}

LensInterface::LensInterface(Transform &c2o, float radius, float n_i, float n_t, float aperture) {
	this->c2o = c2o;
	this->radius = radius;
	this->n_i = n_i;
	this->n_t = n_t;
	this->aperture = aperture;
}

LensInterface::~LensInterface() {
}

bool LensInterface::Intersect(const Ray &r, float *tHit) const {

    Point phit;
    // Transform _Ray_ to object space
    Ray ray;
    (c2o)(r, &ray);

	float thit;
	if (radius != 0) {

    	// Compute quadratic sphere coefficients
	    float A = ray.d.x*ray.d.x + ray.d.y*ray.d.y + ray.d.z*ray.d.z;
	    float B = 2 * (ray.d.x*ray.o.x + ray.d.y*ray.o.y + ray.d.z*ray.o.z);
	    float C = ray.o.x*ray.o.x + ray.o.y*ray.o.y +
	              ray.o.z*ray.o.z - radius*radius;

	    // Solve quadratic equation for _t_ values
	    float t0, t1;
	    if (!Quadratic(A, B, C, &t0, &t1))
	        return false;

	    // Compute intersection distance along ray
	    if (t0 > ray.maxt || t1 < ray.mint)
	        return false;

		if (radius < 0) {
		    thit = t0;	
			if (thit < ray.mint) return false;
		}
		else if (radius > 0){
			thit = t1;
			if (thit > ray.maxt) return false;
		}
	}
	else {
		thit = (0.f - ray.o.z)/ray.d.z;
		if (thit > ray.maxt || thit < ray.mint) return false;
	}
	phit = ray(thit);

	// Do not include intersections outside the aperature
	if (sqrt(pow(phit.x, 2) + pow(phit.y, 2)) > aperture/2.f)
		return false;

    // Update _tHit_ for quadric intersection
    *tHit = thit;
    return true;
}

bool LensInterface::RefractRay(const Ray &ray, Ray &rayOut) const {
	
	// Compute ray-sphere intersection
	float tHit;
	if (!Intersect(ray, &tHit)) return false;

	// Compute refracted ray
	Point P = ray(tHit);
	Vector I = Normalize(ray.d);
	Vector N;
	Point center;
	(Inverse(c2o))(Point(0.f, 0.f, 0.f), &center);
	if (radius < 0) N = Normalize(P - center);
	else if (radius > 0) N = Normalize(center - P);
	else N = Vector(0.f, 0.f, -1.f);
	
	// Snell's Law
	float mu = n_i/n_t;
	float coeff = 1 - pow(mu, 2)*(1- pow(Dot(I, N), 2));
	if (coeff < 0) return false; // Total internal reflection
	
	float gamma = -mu*Dot(I, N) - sqrt(coeff);
	Vector T = Normalize(mu*I + gamma*N);
	rayOut.o = P;
	rayOut.d = T;
	
	return true;
}

RealisticCamera *CreateRealisticCamera(const ParamSet &params,
        const AnimatedTransform &cam2world, Film *film) {
	   // Extract common camera parameters from \use{ParamSet}
	   float hither = params.FindOneFloat("hither", -1);
	   float yon = params.FindOneFloat("yon", -1);
	   float shutteropen = params.FindOneFloat("shutteropen", -1);
	   float shutterclose = params.FindOneFloat("shutterclose", -1);

	   // Realistic camera-specific parameters
	   string specfile = params.FindOneString("specfile", "");
	   float filmdistance = params.FindOneFloat("filmdistance", 70.0); // about 70 mm default to film
	   float fstop = params.FindOneFloat("aperture_diameter", 1.0);
	   float filmdiag = params.FindOneFloat("filmdiag", 35.0);
	   assert(hither != -1 && yon != -1 && shutteropen != -1 &&
	      shutterclose != -1 && filmdistance!= -1);
	   if (specfile == "") {
	       Severe( "No lens spec file supplied!\n" );
	   }
	   return new RealisticCamera(cam2world, hither, yon,
	      shutteropen, shutterclose, filmdistance, fstop,
	      specfile, filmdiag, film);
}

RealisticCamera::RealisticCamera(const AnimatedTransform &cam2world,
                                 float hither, float yon,
                                 float sopen, float sclose,
                                 float filmdistance, float aperture_diameter_,
                                 const string &specfile,
                                 float filmdiag,
								 Film *f)
                                 : Camera(cam2world, sopen, sclose, f),
								   ShutterOpen(sopen),
								   ShutterClose(sclose),
								   film(f)
{
	ifstream file(specfile.c_str());
	string line;
	float distance = 0.f;
	float n_t = 1.f;
	
	// Build lens vector
	while(getline(file, line))
	{
		if (line[0] == '#') continue;
		
	    stringstream linestream(line);
		float radius;
		float thickness;
		float n_i;
		float aperture;

		linestream >> radius >> thickness >> n_i >> aperture;
		if (n_i == 0) n_i = 1.f;
		if (radius == 0) aperture = aperture_diameter_;
		
		Transform c2o = Translate(Vector(0.f, 0.f, distance+radius));
		LensInterface *lens = new LensInterface(c2o, radius, n_i, n_t, aperture);
		lenses.push_back(lens);
		distance += thickness;
		n_t = n_i;
	}
	
	// Set attributes
	distance += filmdistance;
	filmZ = -distance;
	filmDist = filmdistance;
	filmDiag = filmdiag;
}


RealisticCamera::~RealisticCamera()
{

}

	
float RealisticCamera::GenerateRay(const CameraSample &sample, Ray *ray) const
{
	
	// Sample film plane
	float scale = filmDiag/sqrt(pow(film->xResolution,2) + pow(film->yResolution,2));
	float x = film->xResolution*scale/2 - sample.imageX*scale;
	float y = sample.imageY*scale - film->yResolution*scale/2;
	Point P = Point(x, y, filmZ);
	
	// Sample back lens for ray direction
	float xp = 2*sample.lensU - 1;
	float yp = 2*sample.lensV - 1;
	float theta = xp/yp;
	float r = yp * lenses[lenses.size()-1]->Aperture()/2.f;
	Point Pp = Point(0.f, 0.f, filmZ+filmDist) + r*Vector(cos(theta), sin(theta), 0.f);
	
	// Generate ray
	Ray rayIn = Ray(P, Normalize(Pp - P), 0.f);
	Ray rayOut;

	// Trace ray through lenses
	for (int i = lenses.size() - 1; i >= 0; i--) {

		LensInterface *lens = lenses[i];
		float thit;
		if (!lens->Intersect(rayIn, &thit)) return 0.f;
		if (!lens->RefractRay(rayIn, rayOut)) return 0.f;
		rayIn = rayOut;
	}
	
	// Transform ray to world space coordinates
	(CameraToWorld)(rayOut, ray);
	ray->d = Normalize(ray->d);
	
	// Return weight
	float A = M_PI*pow(lenses[lenses.size()-1]->Aperture()/2.f, 2);
	return A*pow(Dot(Normalize(Pp - P), Vector(0.f, 0.f, 1.f)), 4) / pow(filmDist, 2);
}
