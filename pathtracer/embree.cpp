#include "embree.h"
#include <iostream>
#include <map>


using namespace std;
using namespace glm;

namespace pathtracer
{
///////////////////////////////////////////////////////////////////////////
// Global variables
///////////////////////////////////////////////////////////////////////////
RTCDevice embree_device;
RTCScene embree_scene;

///////////////////////////////////////////////////////////////////////////
// Build an acceleration structure for the scene
///////////////////////////////////////////////////////////////////////////
void buildBVH()
{

	//cout << "Embree building BVH..." << flush;
	rtcCommitScene(embree_scene);
	//cout << "done.\n";
}

///////////////////////////////////////////////////////////////////////////
// Called when there is an embree error
///////////////////////////////////////////////////////////////////////////
void embreeErrorHandler(void* unused, enum RTCError code, const char* str)
{
	cout << "Embree ERROR: " << str << endl;
	exit(1);
}

///////////////////////////////////////////////////////////////////////////
// Used to map an Embree geometry ID to our scene Meshes and Materials
///////////////////////////////////////////////////////////////////////////
map<uint32_t, const helpers::Model*> map_geom_ID_to_model;
map<uint32_t, const helpers::Mesh*> map_geom_ID_to_mesh;
map<uint32_t, RTCGeometry> map_geom_ID_to_geom;


///////////////////////////////////////////////////////////////////////////
// Add a model to the embree scene
///////////////////////////////////////////////////////////////////////////
void addModel(helpers::Model* model, const mat4& model_matrix)
{
	///////////////////////////////////////////////////////////////////////
	// Lazy initialize embree on first use
	///////////////////////////////////////////////////////////////////////
	cout << "Initializing embree..." << flush;
	static bool embree_is_initialized = false;
	if(!embree_is_initialized)
	{
		embree_is_initialized = true;
		embree_device = rtcNewDevice(NULL);
		rtcSetDeviceErrorFunction(embree_device, embreeErrorHandler, NULL);
        embree_scene = rtcNewScene(embree_device);
        rtcSetSceneFlags(embree_scene, RTC_SCENE_FLAG_DYNAMIC);
		//embree_scene = rtcDeviceNewScene(embree_device, RTC_SCENE_FLAG_DYNAMIC, RTC_INTERSECT1);
	}
	cout << "done.\n";
	
	///////////////////////////////////////////////////////////////////////
	// Transform and add each mesh in the model as a geometry in embree,
	// and create mappings so that we can connect an embree geom_ID to a
	// Material.
	///////////////////////////////////////////////////////////////////////
	cout << "Adding " << model->m_name << " to embree scene..." << flush;
    
	for(auto& mesh : model->m_meshes)
	{

        RTCGeometry geometry = rtcNewGeometry(embree_device, RTC_GEOMETRY_TYPE_TRIANGLE);
        
		// Transform vertices
        vec4* embree_vertices = (vec4*)rtcSetNewGeometryBuffer(geometry,
                                                               RTC_BUFFER_TYPE_VERTEX,
                                                               0,
                                                               RTC_FORMAT_FLOAT3,
                                                               sizeof(vec4),
                                                               mesh.m_number_of_vertices
        );
		for(uint32_t i = 0; i < mesh.m_number_of_vertices; i++)
		{
			embree_vertices[i] = model_matrix * vec4(model->m_positions[mesh.m_start_index + i], 1.0f);
		}
		
		// Add triangle indices
		int* embree_tri_idxs = (int*)rtcSetNewGeometryBuffer(geometry,
                                                             RTC_BUFFER_TYPE_INDEX,
                                                             0,
                                                             RTC_FORMAT_UINT3,
                                                             sizeof(int)*3,
                                                             mesh.m_number_of_vertices
        );
		for(uint32_t i = 0; i < mesh.m_number_of_vertices; i++)
		{
			embree_tri_idxs[i] = i;
		}
		
		// Commit
		rtcCommitGeometry(geometry);
        unsigned int geom_ID = rtcAttachGeometry(embree_scene, geometry);
		rtcReleaseGeometry(geometry);
        map_geom_ID_to_mesh[geom_ID] = &mesh;
		map_geom_ID_to_model[geom_ID] = model;
		mesh.geom = geometry;

	}
	cout << "done.\n";
}

void moveModel(const helpers::Model* model, const mat4& model_matrix) {
	for (auto& mesh : model->m_meshes) {
        RTCGeometry geometry = mesh.geom;
		vec4* embree_vertices = (vec4*)rtcGetGeometryBufferData(geometry, RTC_BUFFER_TYPE_VERTEX, 0);
		for (uint32_t i = 0; i < mesh.m_number_of_vertices; i++)
		{
			embree_vertices[i] = model_matrix * vec4(model->m_positions[mesh.m_start_index + i], 1.0f);
		}
		
		rtcUpdateGeometryBuffer(geometry, RTC_BUFFER_TYPE_VERTEX, 0);
		rtcCommitGeometry(geometry);
		buildBVH(); // Need to rebuild bvh each time something is moved.
	}
	
}


///////////////////////////////////////////////////////////////////////////
// Extract an intersection from an embree ray.
///////////////////////////////////////////////////////////////////////////
Intersection getIntersection(const Ray_t& rh)
{
    Hit h = rh.hit;
    Ray r = rh.ray;
	const helpers::Model* model = map_geom_ID_to_model[h.geomID];
	const helpers::Mesh* mesh = map_geom_ID_to_mesh[h.geomID];
	Intersection i;
	i.material = &(model->m_materials[mesh->m_material_idx]);
	i.model = model; // new
	vec3 n0 = model->m_normals[((mesh->m_start_index / 3) + h.primID) * 3 + 0];
	vec3 n1 = model->m_normals[((mesh->m_start_index / 3) + h.primID) * 3 + 1];
	vec3 n2 = model->m_normals[((mesh->m_start_index / 3) + h.primID) * 3 + 2];
	float w = 1.0f - (h.u + h.v);
	i.shading_normal = normalize(w * n0 + h.u * n1 + h.v * n2);
	i.geometry_normal = -normalize(h.n);
	i.position = r.o + r.tfar * r.d;
	i.wo = normalize(-r.d);
	return i;
}

///////////////////////////////////////////////////////////////////////////
// Test a ray against the scene and find the closest intersection
///////////////////////////////////////////////////////////////////////////
void filterFunction(const struct RTCFilterFunctionNArguments* args) {
    
}
bool intersect(Ray_t& rh)
{
    RTCIntersectContext context;
    rtcInitIntersectContext(&context);
	rtcIntersect1(embree_scene, &context, (RTCRayHit*)&rh);
	return rh.hit.geomID != RTC_INVALID_GEOMETRY_ID;
}

///////////////////////////////////////////////////////////////////////////
// Test whether a ray is intersected by the scene (do not return an
// intersection).
///////////////////////////////////////////////////////////////////////////
bool occluded(Ray& r)
{
	rtcOccluded1(embree_scene, NULL, (RTCRay*)&r);
    /* TODO: this condition has not been tested since 
     * only intersection testing is currently used. tfar should be -inf
     * if occluded according to docs
     */
    return r.tfar < 0.0f;
	
}
} // namespace pathtracer
