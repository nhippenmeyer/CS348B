
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


// materials/matte.cpp*
#include "stdafx.h"
#include "materials/matte.h"
#include "paramset.h"
#include "reflection.h"
#include "diffgeom.h"
#include "texture.h"

// MatteMaterial Method Definitions
BSDF *MatteMaterial::GetBSDF(const DifferentialGeometry &dgGeom,
                             const DifferentialGeometry &dgShading,
                             MemoryArena &arena) const {
    // Allocate _BSDF_, possibly doing bump mapping with _bumpMap_
    DifferentialGeometry dgs;
    if (bumpMap)
        Bump(bumpMap, dgGeom, dgShading, &dgs);
    else
        dgs = dgShading;
    BSDF *bsdf = BSDF_ALLOC(arena, BSDF)(dgs, dgGeom.nn);

    // Evaluate textures for _MatteMaterial_ material and allocate BRDF
    Spectrum r = Kd->Evaluate(dgs).Clamp();
    float sig = Clamp(sigma->Evaluate(dgs), 0.f, 90.f);
    if (sig == 0.)
        bsdf->Add(BSDF_ALLOC(arena, Lambertian)(r));
    else
        bsdf->Add(BSDF_ALLOC(arena, OrenNayar)(r, sig));
    return bsdf;
}


MatteMaterial *CreateMatteMaterial(const Transform &xform,
        const TextureParams &mp) {
    Reference<Texture<Spectrum> > Kd = mp.GetSpectrumTexture("Kd", Spectrum(0.5f));
    Reference<Texture<float> > sigma = mp.GetFloatTexture("sigma", 0.f);
    Reference<Texture<float> > bumpMap = mp.GetFloatTexture("bumpmap", 0.f);
    return new MatteMaterial(Kd, sigma, bumpMap);
}


