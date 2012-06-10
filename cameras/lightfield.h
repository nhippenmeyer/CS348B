#include "core/geometry.h"

class Plane {
public:
	Plane(float z, float width, float height, int xRes, int yRes);
	bool Intersect(Ray& ray, int* coords);
	float z, width, height;
	int xRes, yRes;
};

class UVPlane : public Plane {
public:
	UVPlane(float z, float width, float height, int xRes, int yRes);
};

class STPlane : public Plane {
public:
	STPlane(float z, float width, float height, int xRes, int yRes);
	short **r, **g, **b;
};

class LightField {
	
public:
	LightField(float zClose, float widthClose, float heightClose, int xResClose, int yResClose,
	 		   float zFar, float widthFar, float heightFar, int xResFar, int yResFar);
    LightField(char * fileName);  //loads data for LightField in from file
	bool Intersect(Ray& r, float* rgb);	
	void AddRayToField(Ray& r, float* rgb);
	bool writeToFile();
    bool loadFromFile();
private:
	UVPlane* uvPlane;
	STPlane** stPlanes;
};
