#pragma once
#include <embree3/rtcore.h>
#include <embree3/rtcore_ray.h>
#include <embree3/rtcore_common.h>
#include "Model.h"
#include <glm/glm.hpp>
#include <map>

#define COUNTOF(x) (sizeof(x)/sizeof(*x))

namespace pathtracer
{
///////////////////////////////////////////////////////////////////////////
// Add a model to the embree scene
///////////////////////////////////////////////////////////////////////////
void addModel(helpers::Model* model, const glm::mat4& model_matrix);
// move model
void moveModel(const helpers::Model* model, const glm::mat4& model_matrix);

///////////////////////////////////////////////////////////////////////////
// Build an acceleration structure for the scene
///////////////////////////////////////////////////////////////////////////
void buildBVH();

///////////////////////////////////////////////////////////////////////////
// This struct is what an embree Ray must look like. It contains the
// information about the ray to be shot and (after intersect() has been
// called) the geometry the ray hit.
///////////////////////////////////////////////////////////////////////////
struct RTC_ALIGN(16) Ray
{
	Ray(const glm::vec3& origin = glm::vec3(0.0f),
	    const glm::vec3& direction = glm::vec3(0.0f),
	    float near = 0.0f,
	    float far = FLT_MAX,
        uint32_t id_n = 0,
        uint32_t ff = 0)
	    : o(origin), d(direction), tnear(near), tfar(far), id(id_n), flags(ff)
        {}
	// Ray data
	glm::vec3 o;
    float tnear = 0.0f;
    
	glm::vec3 d;
	float time = 0.0f;
    
	float tfar = FLT_MAX;
    uint32_t mask = 0xFFFFFFFF;
    uint32_t id;
    uint32_t flags;
	
	
};

struct RTC_ALIGN(16) Hit
{
    Hit(const glm::vec3& normal = glm::vec3(0.0f))
        : n(normal)
        
    {
        geomID = RTC_INVALID_GEOMETRY_ID;
		primID = RTC_INVALID_GEOMETRY_ID;
        for(uint32_t index = 0; index < COUNTOF(instID); index++) { instID[index] = RTC_INVALID_GEOMETRY_ID; }
	}
    // Hit Data
	glm::vec3 n;
	float u, v;
    
	
    unsigned int primID = RTC_INVALID_GEOMETRY_ID;
	unsigned int geomID = RTC_INVALID_GEOMETRY_ID;
    unsigned int instID[RTC_MAX_INSTANCE_LEVEL_COUNT];
};


struct RTC_ALIGN(16) Ray_t
{
    Ray_t(const glm::vec3& origin = glm::vec3(0.0f),
	     const glm::vec3& direction = glm::vec3(0.0f))
    {
        ray = Ray(origin, direction);
        hit = Hit();
    }
    struct Ray ray;
    struct Hit hit;
};

///////////////////////////////////////////////////////////////////////////
// This struct describes an intersection, as extracted from the Embree
// ray.
///////////////////////////////////////////////////////////////////////////
struct Intersection
{
	glm::vec3 position;
	glm::vec3 geometry_normal;
	glm::vec3 shading_normal;
	glm::vec3 wo;
	const helpers::Material* material;
	const helpers::Model* model;
};
Intersection getIntersection(const Ray_t&);

///////////////////////////////////////////////////////////////////////////
// Test a ray against the scene and find the closest intersection
///////////////////////////////////////////////////////////////////////////
bool intersect(Ray_t&);

///////////////////////////////////////////////////////////////////////////
// Test whether a ray is intersected by the scene (do not return an
// intersection).
///////////////////////////////////////////////////////////////////////////
bool occluded(const Ray&);
} // namespace pathtracer
