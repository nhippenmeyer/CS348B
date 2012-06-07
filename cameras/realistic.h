#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_CAMERAS_REALISTIC_H
#define PBRT_CAMERAS_REALISTIC_H

#include "pbrt.h"
#include "camera.h"
#include "film.h"


class LensInterface {
public:
	LensInterface(Transform &c2o, float radius, float n_i, float n_t, float aperature);
	~LensInterface();
    bool Intersect(const Ray &ray, float *tHit) const;
	bool RefractRay(const Ray &ray, Ray &rayOut) const;
	float Aperture() {return aperture;}
private:
	Transform c2o;
	float radius;
	float n_i, n_t;
	float aperture;
};

class RealisticCamera : public Camera {
public:
   RealisticCamera(const AnimatedTransform &cam2world,
      float hither, float yon, float sopen,
      float sclose, float filmdistance, float aperture_diameter,
      const string &specfile,
      float filmdiag,
	  Film *film);
	~RealisticCamera();
	float GenerateRay(const CameraSample &sample, Ray *ray) const;
	void GenerateCameraRay(const CameraSample &sample, Ray *ray) const;

private:
	float ShutterOpen;
	float ShutterClose;
	Film * film;
	vector<LensInterface* > lenses;
	float filmZ;
	float filmDist;
	float filmDiag;
};

RealisticCamera *CreateRealisticCamera(const ParamSet &params,
        const AnimatedTransform &cam2world, Film *film);

#endif
