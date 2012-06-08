#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_CAMERAS_REALISTIC_H
#define PBRT_CAMERAS_REALISTIC_H

#include "pbrt.h"
#include "camera.h"
#include "film.h"

struct LensElement{
    float radius;
    float separation;
    float n;
    float aperture;
};

class RealisticCamera : public Camera {
public:
   RealisticCamera(const AnimatedTransform &cam2world,
      float hither, float yon, float sopen,
      float sclose, float filmdistance, float aperture_diameter,
      const string &specfile,
      float filmdiag,
	  Film *film, bool preproc);
   ~RealisticCamera();
   float GenerateRay(const CameraSample &sample, Ray *) const;
   void GenerateCameraRay(const CameraSample &sample, Ray *) const;
    void runLensFlare(const Scene * scene, const Renderer * renderer) const;
private:
   float ShutterOpen;
   float ShutterClose;
   float filmDiag;
   float filmDistance;
   Film * film;
   vector<LensElement> lensEls;
   bool IntersectLensEl(const Ray &r, float *tHit, float radius, float dist,  Vector & normalVec) const


;
};

RealisticCamera *CreateRealisticCamera(const ParamSet &params,
        const AnimatedTransform &cam2world, Film *film);

#endif
