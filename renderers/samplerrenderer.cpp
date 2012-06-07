
/*
    pbrt source code Copyright(c) 1998-2010 Matt Pharr and Greg Humphreys.

    This file is part of pbrt.

    pbrt is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.  Note that the text contents of
    the book "Physically Based Rendering" are *not* licensed under the
    GNU GPL.

    pbrt is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

 */


// renderers/samplerrenderer.cpp*
#include "stdafx.h"
#include "renderers/samplerrenderer.h"
#include "scene.h"
#include "film.h"
#include "volume.h"
#include "sampler.h"
#include "integrator.h"
#include "progressreporter.h"
#include "camera.h"
#include "intersection.h"
#include <iostream>
#include <fstream>
#include <iomanip> 

using namespace std;


static uint32_t hash(char *key, uint32_t len)
{
    uint32_t   hash, i;
    for (hash=0, i=0; i<len; ++i) {
        hash += key[i];
        hash += (hash << 10);
        hash ^= (hash >> 6);
    }
    hash += (hash << 3);
    hash ^= (hash >> 11);
    hash += (hash << 15);
    return hash;
} 

// SamplerRendererTask Definitions
void SamplerRendererTask::Run() {
	
    PBRT_STARTED_RENDERTASK(taskNum);
    // Get sub-_Sampler_ for _SamplerRendererTask_
    Sampler *sampler = mainSampler->GetSubSampler(taskNum, taskCount);
    if (!sampler)
    {
        reporter.Update();
        PBRT_FINISHED_RENDERTASK(taskNum);
        return;
    }

    // Declare local variables used for rendering loop
    MemoryArena arena;
    RNG rng(taskNum);

    // Allocate space for samples and intersections
    int maxSamples = sampler->MaximumSampleCount();
    Sample *samples = origSample->Duplicate(maxSamples);
    RayDifferential *rays = new RayDifferential[maxSamples];
    Spectrum *Ls = new Spectrum[maxSamples];
    Spectrum *Ts = new Spectrum[maxSamples];
    Intersection *isects = new Intersection[maxSamples];

	FILE * lightfieldBin;
	ofstream lightfield;

	// if (preprocess) {
	// 	// Create lightfield file
	// 	lightfieldBin = fopen("lightfield.bin", "a");
	//   	lightfield.open("lightfield.txt", ios::out | ios::app );
	// } 

    // Get samples from _Sampler_ and update image
    int sampleCount;
    while ((sampleCount = sampler->GetMoreSamples(samples, rng)) > 0) {
        // Generate camera rays and compute radiance along rays
        for (int i = 0; i < sampleCount; ++i) {
            // Find camera ray for _sample[i]_
            PBRT_STARTED_GENERATING_CAMERA_RAY(&samples[i]);
			float rayWeight = camera->GenerateRayDifferential(samples[i], &rays[i]);
            rays[i].ScaleDifferentials(1.f / sqrtf(sampler->samplesPerPixel));
            PBRT_FINISHED_GENERATING_CAMERA_RAY(&samples[i], &rays[i], rayWeight);

			Ray* cameraRay = new Ray();
			camera->GenerateCameraRay(samples[i], cameraRay);

			if (preprocess) {

            	// Evaluate radiance along camera ray
	            PBRT_STARTED_CAMERA_RAY_INTEGRATION(&rays[i], &samples[i]);
	            if (visualizeObjectIds) {
	                if (rayWeight > 0.f && scene->Intersect(rays[i], &isects[i])) {
	                    // random shading based on shape id...
	                    uint32_t ids[2] = { isects[i].shapeId, isects[i].primitiveId };
	                    uint32_t h = hash((char *)ids, sizeof(ids));
	                    float rgb[3] = { (h & 0xff), (h >> 8) & 0xff, (h >> 16) & 0xff };
	                    Ls[i] = Spectrum::FromRGB(rgb);
	                    Ls[i] /= 255.f;
						camera->lightfield->AddRayToField(*cameraRay, rgb);

	                }
	                else
	                    Ls[i] = 0.f;
	            }
	            else {
	            	if (rayWeight > 0.f) 
		                Ls[i] = rayWeight * renderer->Li(scene, rays[i], &samples[i], rng,
		                                                 arena, &isects[i], &Ts[i]);
		            else {
		                Ls[i] = 0.f;
		                Ts[i] = 1.f;
		            }
				}
				
				/*
				short x = (short)(samples[i].imageX);
				short y = (short)(samples[i].imageY);                                                                                                                                                                                                             
				short u = (short)(samples[i].lensU * 256);
				short v = (short)(samples[i].lensV * 256);
				short r = (short)(Ls[i].GetCoeff(0) * 256);
				short g = (short)(Ls[i].GetCoeff(1) * 256);
				short b = (short)(Ls[i].GetCoeff(2) * 256);
				fwrite((void*)(&x), sizeof(x), 1, lightfieldBin);
				fwrite((void*)(&y), sizeof(y), 1, lightfieldBin);
				fwrite((void*)(&u), sizeof(u), 1, lightfieldBin);
				fwrite((void*)(&v), sizeof(v), 1, lightfieldBin);
				lightfield << x << "\t" << y << "\t" << u << "\t" << v;
				if (r > 0 && g > 0 && b > 0) {
					fwrite((void*)(&r), sizeof(r), 1, lightfieldBin);
					fwrite((void*)(&g), sizeof(g), 1, lightfieldBin);
					fwrite((void*)(&b), sizeof(b), 1, lightfieldBin);
					lightfield << "\t" << r << "\t" << g << "\t" << b;
				}
				lightfield << endl;
				*/
			} else {
				/*
				float rgb[3];
				rgb[0] = 0.f;
				rgb[1] = 0.f;
				rgb[2] = 1.f;
				Ls[i] = Spectrum::FromRGB(rgb);
				string line;
			  	ifstream lightfield("lightfield.txt");
		  		if (lightfield.is_open()) {
				   while (lightfield.good()) {
				      	getline (lightfield,line);
					  	char * pch;
						float rgb[3];
						rgb[0] = 0.f;
						rgb[1] = 0.f;
						rgb[2] = 0.f;
						// pch = strtok((char *)line.c_str(), ",");
						// rgb[0] = atof(pch);
						// for (int j = 1; j < 3; j++) {
						// 	pch = strtok(NULL, ",");
						// 	rgb[j] = atof(pch);
						// }
						Ls[i] = Spectrum::FromRGB(rgb);
				    }
				    lightfield.close();
				}
				*/
				float* rgb = new float[3];
				if (!camera->lightfield->Intersect(*cameraRay, rgb))
					rgb[0] = rgb[1] = rgb[2] = 0.f;
	            Ls[i] = Spectrum::FromRGB(rgb);
	            Ls[i] /= 255.f;
			}
			
            // Issue warning if unexpected radiance value returned
            if (Ls[i].HasNaNs()) {
                Error("Not-a-number radiance value returned "
                      "for image sample.  Setting to black.");
                Ls[i] = Spectrum(0.f);
            }
            else if (Ls[i].y() < -1e-5) {
                Error("Negative luminance value, %f, returned"
                      "for image sample.  Setting to black.", Ls[i].y());
                Ls[i] = Spectrum(0.f);
            }
            else if (isinf(Ls[i].y())) {
                Error("Infinite luminance value returned"
                      "for image sample.  Setting to black.");
                Ls[i] = Spectrum(0.f);
            }
            PBRT_FINISHED_CAMERA_RAY_INTEGRATION(&rays[i], &samples[i], &Ls[i]);
        }

        // Report sample results to _Sampler_, add contributions to image
        if (sampler->ReportResults(samples, rays, Ls, isects, sampleCount))
        {
            for (int i = 0; i < sampleCount; ++i)
            {
                PBRT_STARTED_ADDING_IMAGE_SAMPLE(&samples[i], &rays[i], &Ls[i], &Ts[i]);
                camera->film->AddSample(samples[i], Ls[i]);
                PBRT_FINISHED_ADDING_IMAGE_SAMPLE();
            }
        }

        // Free _MemoryArena_ memory from computing image sample values
        arena.FreeAll();
    }
	
	// if (preprocess) {
	// 	fclose(lightfieldBin);
	// 	lightfield.close();
	// }

    // Clean up after _SamplerRendererTask_ is done with its image region
    camera->film->UpdateDisplay(sampler->xPixelStart,
        sampler->yPixelStart, sampler->xPixelEnd+1, sampler->yPixelEnd+1);
    delete sampler;
    delete[] samples;
    delete[] rays;
    delete[] Ls;
    delete[] Ts;
    delete[] isects;
    reporter.Update();
    PBRT_FINISHED_RENDERTASK(taskNum);
}



// SamplerRenderer Method Definitions
SamplerRenderer::SamplerRenderer(Sampler *s, Camera *c,
                                 SurfaceIntegrator *si, VolumeIntegrator *vi,
                                 bool visIds) {
    sampler = s;
    camera = c;
    surfaceIntegrator = si;
    volumeIntegrator = vi;
    visualizeObjectIds = visIds;
}


SamplerRenderer::~SamplerRenderer() {
    delete sampler;
    delete camera;
    delete surfaceIntegrator;
    delete volumeIntegrator;
}


void SamplerRenderer::Render(const Scene *scene) {
    PBRT_FINISHED_PARSING();
    // Allow integrators to do preprocessing for the scene
    PBRT_STARTED_PREPROCESSING();
    surfaceIntegrator->Preprocess(scene, camera, this);
    volumeIntegrator->Preprocess(scene, camera, this);
    PBRT_FINISHED_PREPROCESSING();
    PBRT_STARTED_RENDERING();
    // Allocate and initialize _sample_
    Sample *sample = new Sample(sampler, surfaceIntegrator,
                                volumeIntegrator, scene);

    // Create and launch _SamplerRendererTask_s for rendering image

    // Compute number of _SamplerRendererTask_s to create for rendering
    int nPixels = camera->film->xResolution * camera->film->yResolution;
    int nTasks = max(32 * NumSystemCores(), nPixels / (16*16));
    nTasks = RoundUpPow2(nTasks);
    ProgressReporter reporter(nTasks, "Rendering");

    vector<Task *> preprocessTasks;
    for (int i = 0; i < nTasks; ++i)
        preprocessTasks.push_back(new SamplerRendererTask(scene, this, camera,
                                                      	  reporter, sampler, sample, 
                                                      	  visualizeObjectIds, 
                                                      	  nTasks-1-i, nTasks, true));
    EnqueueTasks(preprocessTasks);
    WaitForAllTasks();
    for (uint32_t i = 0; i < preprocessTasks.size(); ++i)
        delete preprocessTasks[i];
    reporter.Done();

	vector<Task *> renderTasks;
    for (int i = 0; i < nTasks; ++i)
        renderTasks.push_back(new SamplerRendererTask(scene, this, camera,
                                                      reporter, sampler, sample, 
                                                      visualizeObjectIds, 
                                                      nTasks-1-i, nTasks, false));
    EnqueueTasks(renderTasks);
    WaitForAllTasks();
    for (uint32_t i = 0; i < renderTasks.size(); ++i)
        delete renderTasks[i];

    PBRT_FINISHED_RENDERING();
    // Clean up after rendering and store final image
    delete sample;
    camera->film->WriteImage();
}


Spectrum SamplerRenderer::Li(const Scene *scene,
        const RayDifferential &ray, const Sample *sample, RNG &rng,
        MemoryArena &arena, Intersection *isect, Spectrum *T) const {
    Assert(ray.time == sample->time);
    Assert(!ray.HasNaNs());
    // Allocate local variables for _isect_ and _T_ if needed
    Spectrum localT;
    if (!T) T = &localT;
    Intersection localIsect;
    if (!isect) isect = &localIsect;
    Spectrum Li = 0.f;
    if (scene->Intersect(ray, isect))
        Li = surfaceIntegrator->Li(scene, this, ray, *isect, sample,
                                   rng, arena);
    else {
        // Handle ray that doesn't intersect any geometry
        for (uint32_t i = 0; i < scene->lights.size(); ++i)
           Li += scene->lights[i]->Le(ray);
    }
    Spectrum Lvi = volumeIntegrator->Li(scene, this, ray, sample, rng,
                                        T, arena);
    return *T * Li + Lvi;
}


Spectrum SamplerRenderer::Transmittance(const Scene *scene,
        const RayDifferential &ray, const Sample *sample, RNG &rng,
        MemoryArena &arena) const {
    return volumeIntegrator->Transmittance(scene, this, ray, sample,
                                           rng, arena);
}


