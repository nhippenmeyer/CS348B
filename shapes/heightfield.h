
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

#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_SHAPES_HEIGHTFIELD_H
#define PBRT_SHAPES_HEIGHTFIELD_H

// shapes/heightfield.h*
#include "shape.h"
#include "shapes/trianglemesh.h"
#include <vector>
using std::vector;

// Heightfield Declarations
class Heightfield : public Shape {
public:
    // Heightfield Public Methods
    Heightfield(const Transform *o2, const Transform *w2o, bool ro, int nu, int nv, const float *zs);
    ~Heightfield();
    bool CanIntersect() const;
    bool Intersect(const Ray &ray, float *tHit, float *rayEpsilon,
                   DifferentialGeometry *dg) const;
    bool IntersectP(const Ray &ray) const;
    void Refine(vector<Reference<Shape> > &refined) const;
    BBox ObjectBound() const;
    virtual void GetShadingGeometry(const Transform &obj2world,
            const DifferentialGeometry &dg,
            DifferentialGeometry *dgShading) const;
private:
	class QuadNode {
	public:
		QuadNode() {
			isLeaf = false;
			NEchild = SEchild = SWchild = NWchild = NULL;
		}
		bool isLeaf;
		BBox worldBound;
		vector<Reference<Shape> > triangles;
		QuadNode* NEchild;
		QuadNode* SEchild;
		QuadNode* SWchild;
		QuadNode* NWchild;	
	};
	void ConstructQuadTree(QuadNode* current, vector<TriangleMesh* > meshes);
	bool IntersectAccel(QuadNode* node, const Ray &ray, float *tHit, float *rayEpsilon,
				   DifferentialGeometry *dg) const;
	bool IntersectAccelP(QuadNode* node, const Ray &ray) const;
    // Heightfield Private Data
    float *z;
    int nx, ny;
	vector<TriangleMesh* > tmesh;
	QuadNode *quadTree;
};


Heightfield *CreateHeightfieldShape(const Transform *o2w, const Transform *w2o,
        bool reverseOrientation, const ParamSet &params);

#endif // PBRT_SHAPES_HEIGHTFIELD_H
