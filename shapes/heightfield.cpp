
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


// shapes/heightfield.cpp*
#include "stdafx.h"
#include "shapes/heightfield.h"
#include "paramset.h"
#include <iostream>
using std::cout;
using std::endl;


// Heightfield Method Definitions
Heightfield::Heightfield(const Transform *o2w, const Transform *w2o,
        bool ro, int nx, int ny, const float *zs)
    : Shape(o2w, w2o, ro) {
    this->nx = nx;
    this->ny = ny;
    z = new float[nx*ny];
    memcpy(z, zs, nx*ny*sizeof(float));
	
    tmesh.reserve((nx-1)*(ny-1));

	// Construct triangular meshes
    int x, y;
    for (y = 0; y < ny-1; ++y) {
        for (x = 0; x < nx-1; ++x) {

		    int *verts = new int[6];
		    Point *P = new Point[4];
			Normal *N = new Normal[4];
		    float *uvs = new float[8];

		    int pos = 0;
			int dy, dx;
			for (dy = 0; dy <= 1; ++dy) {
				for (dx = 0; dx <= 1; ++dx) {
					P[pos].x = uvs[2*pos]   = (float)(x+dx) / (float)(nx-1);
		            P[pos].y = uvs[2*pos+1] = (float)(y+dy) / (float)(ny-1);
		            P[pos].z = z[(x+dx)+(y+dy)*nx];
		            ++pos;
				}
			}
			
			// Compute average normals at each vertex
			Vector V0 = Cross(P[1]-P[0], P[3]-P[0]) + Cross(P[3]-P[0], P[2]-P[0]);
			Vector V1 = Cross(P[3]-P[1], P[0]-P[1]);
			Vector V2 = Cross(P[0]-P[2], P[3]-P[2]);
			Vector V3 = Cross(P[2]-P[3], P[0]-P[3]) + Cross(P[0]-P[3], P[1]-P[3]);
			
			Point *Pp = new Point[10];
			if (x > 0) {
				Pp[9] = Point((float)(x-1)/(float)(nx-1), (float)(y+1)/(float)(ny-1), z[x-1+(y+1)*nx]);
				Pp[0] = Point((float)(x-1)/(float)(nx-1), (float)y/(float)(ny-1), z[x-1+y*nx]);
				V2 += Cross(Pp[9]-P[2], Pp[0]-P[2]);
				V2 += Cross(Pp[0]-P[2], P[0]-P[2]);
				V0 += Cross(P[2]-P[0], Pp[0]-P[0]);
			}
			if (y > 0) {
				Pp[2] = Point((float)x/(float)(nx-1), (float)(y-1)/(float)(ny-1), z[x+(y-1)*nx]);
				Pp[3] = Point((float)(x+1)/(float)(nx-1), (float)(y-1)/(float)(ny-1), z[x+1+(y-1)*nx]);
				V0 += Cross(Pp[2]-P[0], P[1]-P[0]);
				V1 += Cross(P[0]-P[1], Pp[2]-P[1]);
				V1 += Cross(Pp[2]-P[1], Pp[3]-P[1]);
			}
			if (x > 0 && y > 0) {
				Pp[1] = Point((float)(x-1)/(float)(nx-1), (float)(y-1)/(float)(ny-1), z[x-1+(y-1)*nx]);
				V0 += Cross(Pp[0]-P[0], Pp[1]-P[0]);
				V0 += Cross(Pp[1]-P[0], Pp[2]-P[0]);
			}
			if (x < nx-2) {
				Pp[4] = Point((float)(x+2)/(float)(nx-1), (float)y/(float)(ny-1), z[x+2+y*nx]);
				Pp[5] = Point((float)(x+2)/(float)(nx-1), (float)(y+1)/(float)(ny-1), z[x+2+(y+1)*nx]);
				V1 += Cross(Pp[4]-P[1], Pp[5]-P[1]);
				V1 += Cross(Pp[5]-P[1], P[3]-P[1]);
				V3 += Cross(P[1]-P[3], Pp[5]-P[3]);
			}
			if (x < nx-2 && y > 0) {
				V1 += Cross(Pp[3]-P[1], Pp[4]-P[1]);
			}
			if (y < ny-2) {
				Pp[7] = Point((float)(x+1)/(float)(nx-1), (float)(y+2)/(float)(ny-1), z[x+1+(y+2)*nx]);
				Pp[8] = Point((float)x/(float)(nx-1), (float)(y+2)/(float)(ny-1), z[x+(y+2)*nx]);
				V3 += Cross(Pp[7]-P[3], P[2]-P[3]);
				V2 += Cross(P[3]-P[2], Pp[7]-P[2]);
				V2 += Cross(Pp[7]-P[2], Pp[8]-P[2]);
			}
			if (x < nx-2 && y < ny-2) {
				Pp[6] = Point((float)(x+2)/(float)(nx-1), (float)(y+2)/(float)(ny-1), z[x+2+(y+2)*nx]);
				V3 += Cross(Pp[5]-P[3], Pp[6]-P[3]);
				V3 += Cross(Pp[6]-P[3], Pp[7]-P[3]);
			}
			if (x > 0 && y < ny-2) {
				V2 += Cross(Pp[8]-P[2], Pp[9]-P[2]);
			}

			N[0] = Normal(V0);
			N[1] = Normal(V1);
			N[2] = Normal(V2);
			N[3] = Normal(V3);

			// Each mesh contains 2 connected triangles
		    int *vp = verts;
			*vp++ = 0;
            *vp++ = 1;
            *vp++ = 3;
            *vp++ = 0;
            *vp++ = 3;
            *vp++ = 2;

		    ParamSet paramSet;
		    paramSet.AddInt("indices", verts, 6);
		    paramSet.AddFloat("uv", uvs, 8);
		    paramSet.AddPoint("P", P, 4);
			paramSet.AddNormal("N", N, 4);
		    tmesh.push_back(CreateTriangleMeshShape(ObjectToWorld, WorldToObject, ReverseOrientation, paramSet));

		    delete[] P;
		    delete[] uvs;
		    delete[] verts;
        }
    }

	// Construct acceleration structure
	quadTree = new QuadNode();
	ConstructQuadTree(quadTree, tmesh);
}


Heightfield::~Heightfield() {
    delete[] z;
}


void Heightfield::GetShadingGeometry(const Transform &obj2world,
        const DifferentialGeometry &dg,
        DifferentialGeometry *dgShading) const {
	
	dg.shape->GetShadingGeometry(obj2world, dg, dgShading);
}


void Heightfield::ConstructQuadTree(QuadNode* current, vector<TriangleMesh* > meshes) {
			
	// Base case, construct leaf node		
	if (meshes.size() == 1) {
		meshes[0]->Refine(current->triangles);
		current->isLeaf = true;
		current->worldBound = meshes[0]->WorldBound();
		return;
	}
	
	// Compute bounding boxes for current set of meshes
	BBox objectBound = meshes[0]->ObjectBound();
	BBox worldBound = meshes[0]->WorldBound();
	vector<TriangleMesh* >::const_iterator it;
	for (it=meshes.begin(); it < meshes.end(); ++it) {
		objectBound = Union(objectBound, (*it)->ObjectBound());
		worldBound = Union(worldBound, (*it)->WorldBound());
	}
	float x = (objectBound.pMax.x-objectBound.pMin.x)/2.f + objectBound.pMin.x;
	float y = (objectBound.pMax.y-objectBound.pMin.y)/2.f + objectBound.pMin.y;
	current->worldBound = worldBound;
	
	// Partition meshes into 4 distinct containers
	float eps = 0.0001;
	vector<TriangleMesh* > NEmeshes, SEmeshes, SWmeshes, NWmeshes;
	for (it=meshes.begin(); it < meshes.end(); ++it) {
		BBox bound = (*it)->ObjectBound();
		if (bound.pMax.x > x+eps && bound.pMax.y > y+eps)
			NEmeshes.push_back(*it);
		else if (bound.pMax.x > x+eps && bound.pMin.y < y-eps)
			SEmeshes.push_back(*it);
		else if (bound.pMin.x < x-eps && bound.pMin.y < y-eps)
			SWmeshes.push_back(*it);
		else if (bound.pMin.x < x-eps && bound.pMax.y > y+eps)
			NWmeshes.push_back(*it);
	}
	
	// Construct child nodes with corresponding meshes
	if (NEmeshes.size() > 0) {
		QuadNode* NEchild = new QuadNode();
		ConstructQuadTree(NEchild, NEmeshes);
		current->NEchild = NEchild;
	}                                                
	if (SEmeshes.size() > 0) {    
		QuadNode* SEchild = new QuadNode();          
		ConstructQuadTree(SEchild, SEmeshes);
		current->SEchild = SEchild;
	}                                                
	if (SWmeshes.size() > 0) { 
		QuadNode* SWchild = new QuadNode();          
		ConstructQuadTree(SWchild, SWmeshes);
		current->SWchild = SWchild;
	}                                                
	if (NWmeshes.size() > 0) {    
		QuadNode* NWchild = new QuadNode();          
		ConstructQuadTree(NWchild, NWmeshes);
		current->NWchild = NWchild;
	}
}


BBox Heightfield::ObjectBound() const {
    float minz = z[0], maxz = z[0];
    for (int i = 1; i < nx*ny; ++i) {
        if (z[i] < minz) minz = z[i];
        if (z[i] > maxz) maxz = z[i];
    }
    return BBox(Point(0,0,minz), Point(1,1,maxz));
}


bool Heightfield::CanIntersect() const {
    return true;
}


bool Heightfield::IntersectAccel(QuadNode* node, const Ray &ray, float *tHit, float *rayEpsilon,
               DifferentialGeometry *dg) const {
			
	bool hit = false;
	
	// Base case, check if ray intersects with either triangle
	if (node->isLeaf) {
		vector<Reference<Shape> >::const_iterator triIt;
		for (triIt=node->triangles.begin(); triIt < node->triangles.end(); ++triIt) {
			if ((*triIt)->IntersectP(ray)) {
				float tHitTemp;
				float epsTemp;
				DifferentialGeometry *dgTemp = new DifferentialGeometry();
				(*triIt)->Intersect(ray, &tHitTemp, &epsTemp, dgTemp);
				if (tHitTemp < *tHit) {
					*tHit = tHitTemp;
					*rayEpsilon = epsTemp;
					*dg = *dgTemp;
					hit = true;
				}
				delete dgTemp;
			}
		}
		return hit;
	}
	
	// If ray intersects current bounding box, check for ray intersections on each child node
	if (node->worldBound.IntersectP(ray)) {
		if (node->NEchild && IntersectAccel(node->NEchild, ray, tHit, rayEpsilon, dg)) hit = true;
		if (node->SEchild && IntersectAccel(node->SEchild, ray, tHit, rayEpsilon, dg)) hit = true;
		if (node->SWchild && IntersectAccel(node->SWchild, ray, tHit, rayEpsilon, dg)) hit = true;
		if (node->NWchild && IntersectAccel(node->NWchild, ray, tHit, rayEpsilon, dg)) hit = true;
	}
	return hit;
}


bool Heightfield::IntersectAccelP(QuadNode* node, const Ray &ray) const {
	
	// Base case, return true if ray intersects with either triangle
	if (node->isLeaf) {
		vector<Reference<Shape> >::const_iterator triIt;
		for (triIt=node->triangles.begin(); triIt < node->triangles.end(); ++triIt) {
			if ((*triIt)->IntersectP(ray)) {
				return true;
			}
		}
		return false;
	}
	
	// If ray intersects current bounding box, check for ray intersections on each child node
	if (node->worldBound.IntersectP(ray)) {
		if (node->NEchild && IntersectAccelP(node->NEchild, ray)) return true;
		if (node->SEchild && IntersectAccelP(node->SEchild, ray)) return true;
		if (node->SWchild && IntersectAccelP(node->SWchild, ray)) return true;
		if (node->NWchild && IntersectAccelP(node->NWchild, ray)) return true;
	}
	return false;
}


bool Heightfield::Intersect(const Ray &ray, float *tHit, float *rayEpsilon,
               DifferentialGeometry *dg) const {
	
	*tHit = ray.maxt;
	return IntersectAccel(quadTree, ray, tHit, rayEpsilon, dg);
}


bool Heightfield::IntersectP(const Ray &ray) const {

	return IntersectAccelP(quadTree, ray);
}


void Heightfield::Refine(vector<Reference<Shape> > &refined) const {
    int ntris = 2*(nx-1)*(ny-1);
    refined.reserve(ntris);
    int *verts = new int[3*ntris];
    Point *P = new Point[nx*ny];
    float *uvs = new float[2*nx*ny];
    int nverts = nx*ny;
    int x, y;
    // Compute heightfield vertex positions
    int pos = 0;
    for (y = 0; y < ny; ++y) {
        for (x = 0; x < nx; ++x) {
            P[pos].x = uvs[2*pos]   = (float)x / (float)(nx-1);
            P[pos].y = uvs[2*pos+1] = (float)y / (float)(ny-1);
            P[pos].z = z[pos];
            ++pos;
        }
    }

    // Fill in heightfield vertex offset array
    int *vp = verts;
    for (y = 0; y < ny-1; ++y) {
        for (x = 0; x < nx-1; ++x) {
#define VERT(x,y) ((x)+(y)*nx)
            *vp++ = VERT(x, y);
            *vp++ = VERT(x+1, y);
            *vp++ = VERT(x+1, y+1);
    
            *vp++ = VERT(x, y);
            *vp++ = VERT(x+1, y+1);
            *vp++ = VERT(x, y+1);
        }
#undef VERT
    }
    ParamSet paramSet;
    paramSet.AddInt("indices", verts, 3*ntris);
    paramSet.AddFloat("uv", uvs, 2 * nverts);
    paramSet.AddPoint("P", P, nverts);
    refined.push_back(CreateTriangleMeshShape(ObjectToWorld, WorldToObject, ReverseOrientation, paramSet));
    delete[] P;
    delete[] uvs;
    delete[] verts;
}


Heightfield *CreateHeightfieldShape(const Transform *o2w, const Transform *w2o,
        bool reverseOrientation, const ParamSet &params) {
    int nu = params.FindOneInt("nu", -1);
    int nv = params.FindOneInt("nv", -1);
    int nitems;
    const float *Pz = params.FindFloat("Pz", &nitems);
    Assert(nitems == nu*nv);
    Assert(nu != -1 && nv != -1 && Pz != NULL);
    return new Heightfield(o2w, w2o, reverseOrientation, nu, nv, Pz);
}


