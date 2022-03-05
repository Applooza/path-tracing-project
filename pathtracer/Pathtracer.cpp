#include "Pathtracer.h"
#include <memory>
#include <iostream>
#include <map>
#include <algorithm>
#include "material.h"
#include "embree.h"
#include "sampling.h"

#define MY_FAV_MIN 0.00001f
#define MIN_THROUGHPUT 0.00001f // simply tested a bit, might not be optimal value

using namespace std;
using namespace glm;

namespace pathtracer
{




	///////////////////////////////////////////////////////////////////////////////
	// Global variables
	///////////////////////////////////////////////////////////////////////////////
	Settings settings;
	Environment environment;
	Image rendered_image;
	


	AreaLight::AreaLight(float i, float r, const vec3& pos) {
		intensity_multiplier = i;
		radius = r;
		area = r * r * M_PI;
		color = vec3(1.f);
		position = pos;
		normal = vec3(0.f, -1.f, 0.f); // always assume normal is down
	}



	//vector<AreaLight*> lights;
	map<const helpers::Model*, AreaLight*> lights;

	
	
	// create a new light with some standard data, fix me
	AreaLight* createLight(const helpers::Model* m) {
		float sx, sz;
		concentricSampleDisk(&sx, &sz);
		sx *= 10.f;
		sz *= 10.f;
		AreaLight* newLight = new AreaLight(50000.f, 3.f, vec3(sx, 40.f, sz));
		lights.insert(make_pair(m, newLight));
		return newLight;
	}

	



	///////////////////////////////////////////////////////////////////////////
	// Restart rendering of image
	///////////////////////////////////////////////////////////////////////////
	void restart()
	{
		// No need to clear image,
		
		rendered_image.number_of_samples = 0;
	}

	///////////////////////////////////////////////////////////////////////////
	// On window resize, window size is passed in, actual size of pathtraced
	// image may be smaller (if we're subsampling for speed)
	///////////////////////////////////////////////////////////////////////////
	void resize(int w, int h)
	{
		rendered_image.width = w / settings.subsampling;
		rendered_image.height = h / settings.subsampling;
		rendered_image.data.resize(rendered_image.width * rendered_image.height);
		restart();
	}

	///////////////////////////////////////////////////////////////////////////
	// Return the radiance from a certain direction wi from the environment
	// map.
	///////////////////////////////////////////////////////////////////////////
	vec3 Lenvironment(const vec3& wi)
	{
		const float theta = acos(std::max(-1.0f, std::min(1.0f, wi.y)));
		float phi = atan(wi.z, wi.x);
		if (phi < 0.0f)
			phi = phi + 2.0f * M_PI;
		vec2 lookup = vec2(phi / (2.0 * M_PI), theta / M_PI);
		return environment.multiplier * environment.map.sample(lookup.x, lookup.y);
	}

	////////////////////////////////////////////////////////////////////////
	// Weight functions for multiple importance sampling
	////////////////////////////////////////////////////////////////////////
	inline float BalanceHeuristic(int nf, float fPdf, int ng, float gPdf) {
		return (nf * fPdf) / (nf * fPdf + ng * gPdf);
	}

	inline float PowerHeuristic(uint numf, float fPdf, uint numg, float gPdf) {
		float f = numf * fPdf;
		float g = numg * gPdf;

		return (f * f) / (f * f + g * g);
	}

	// Given a direction, return the probability density for this direction based on the brdf
	float sample_pdf(BRDF& mat, Intersection& h, vec3& wi, int nSamples) {
		float light_scatter = 0.f;
		float light_scatter_sum = 0.f;
		for (int i = 0; i < nSamples; i++) {
			mat.get_pdf(wi, h.wo, h.shading_normal, light_scatter); 
			if(light_scatter >= MY_FAV_MIN)
				light_scatter_sum += light_scatter/nSamples;
		}
		
		return light_scatter_sum;
	}
	
	////////////////////////////////////////////////////////////////////////////
	// Sample a light and set direction and pdf etc
	////////////////////////////////////////////////////////////////////////////
	vec3 sample_light(const AreaLight* lp, BRDF& mat, vec3& wi, const Intersection& hit, float& p, float& r_square) {
		float sx, sz;
		concentricSampleDisk(&sx, &sz);
		vec3 light_sample_point = lp->position + lp->radius * vec3(sx, 0.0f, sz); // this assumes positive y is world up
		r_square = pow(length(light_sample_point - hit.position), 2.f);
		float falloff_factor = 1.0f / r_square;
		vec3 Li = lp->intensity_multiplier * lp->color * falloff_factor; 
		wi = normalize(light_sample_point - hit.position); // almost same as "dir" but without bias
		float dot_l = abs(dot(hit.shading_normal, wi));
		//float dot_l = std::max(0.f, dot(wi, hit.shading_normal));
		p = r_square / (abs(dot(lp->normal, -wi)) * lp->area); // pdf w.r.t solid angle
		//p = (1 / lp.area);	// pdf w.r.t area
		//float g = dot_l * abs(dot(lp.normal, -wi)) / r_square;
		return mat.f(wi, hit.wo, hit.shading_normal) * Li * dot_l;

	}

	// helper function
	bool myAll(const vec3& v, float t) { return (v.x <= t) && (v.y <= t) && (v.z <= t); }

	///////////////////////////////////////////////////////////////////////////
	// Calculate the radiance going from one point (r.hitPosition()) in one
	// direction (-r.d), through path tracing.
	///////////////////////////////////////////////////////////////////////////
	vec3 Li(Ray_t& primary_ray)
	{
		vec3 L = vec3(0.0f);
		vec3 path_throughput = vec3(1.0);
		Ray_t current_ray = primary_ray;


		int bounces;
		for (bounces = 0; bounces < settings.max_bounces; bounces++) {
			///////////////////////////////////////////////////////////////////
			// Get the intersection information from the ray
			///////////////////////////////////////////////////////////////////
			Intersection hit = getIntersection(current_ray);
			vec3 hit_bias = hit.position + hit.shading_normal * EPSILON;
			///////////////////////////////////////////////////////////////////
			// Create a Material tree for evaluating brdfs and calculating
			// sample directions.
			///////////////////////////////////////////////////////////////////
			Diffuse diffuse(hit.material->m_color);
			BlinnPhong dielectric(hit.material->m_shininess, hit.material->m_fresnel, &diffuse);
			BlinnPhongMetal metal(hit.material->m_color, hit.material->m_shininess,
				hit.material->m_fresnel);
			LinearBlend metal_blend(hit.material->m_metalness, &metal, &dielectric);
			LinearBlend reflectivity_blend(hit.material->m_reflectivity, &metal_blend, &diffuse);
			BRDF& mat = reflectivity_blend;

			///////////////////////////////////////////////
			// Emission
			//////////////////////////////////////////////
			L += path_throughput * hit.material->m_emission * hit.material->m_color;


			// path continues along wi
			vec3 wi(0.0f);

			//////////////////////////////////////////////
			// Multiple Importance sampling
			//////////////////////////////////////////////
			// if we intersected a light source the path is done
			if (!hit.model->isLightSource) {
				if (settings.MIS_enabled) {
					//////////////////////////////////////////////
					// Area light sampling
					//////////////////////////////////////////////
					{
						vec3 light_contribution(0.f);
						vec3 wi_l(0.f);
						float light_pdf, r_square;
						// sample each light
						for (auto lightP : lights) {
							// sample a number of times
							for (int i = 0; i < settings.light_samples; i++) {
								vec3 f = sample_light(lightP.second, mat, wi_l, hit, light_pdf, r_square);
								Ray_t shadow_ray(hit_bias, wi_l);
								// with light models we need to actually intersect instead of using occluded() :/
								bool isLight = false;
								if (intersect(shadow_ray)) {
									isLight = getIntersection(shadow_ray).model->isLightSource;
								}
								if (isLight) {
									float light_scatter = sample_pdf(mat, hit, wi_l, 4); // average bsdf pdf of sampling this direction (solid angle)
									float weight = PowerHeuristic(settings.light_samples, light_pdf,
										settings.brdf_samples, light_scatter);
									if (light_scatter >= MY_FAV_MIN)
										// I suppose light_pdf should be checked as well
										light_contribution += f * path_throughput * weight / light_pdf; //Li and cos_term accounted for in f
								}
							}
						}
						if (length(light_contribution) != 0.f) //gets rid of black pixels
							L += light_contribution / (float)settings.light_samples;
					}
					////////////////////////////////////////
					// BRDF Sampling
					////////////////////////////////////////
					float pdf{}, cosine_term{}, weight{};
					vec3 brdf_sample;
					vec3 light_contribution(0.f);
					float lPdf_brdf = 0.0f;
					
					// sample brdf multiple times
					for (int i = 0; i < settings.brdf_samples; i++) {
						brdf_sample = mat.sample_wi(wi, hit.wo, hit.shading_normal, pdf);
						//cosine_term = std::max(0.f, dot(hit.shading_normal, wi));
						cosine_term = abs(dot(hit.shading_normal, wi));
							
						Ray_t ssr(hit_bias, wi);
						bool isLight = false;
						if (!intersect(ssr)) {
							continue;
						}
						Intersection lightHit = getIntersection(ssr);
						if (!lightHit.model->isLightSource || (pdf <= MY_FAV_MIN))
							continue; // we did not hit a light
							
						float r_square = pow(length(lightHit.position - hit.position), 2.f);
						AreaLight* light = lights.at(lightHit.model);
						lPdf_brdf = r_square / (abs(dot(light->normal, -wi)) * light->area);
						
						
						weight = PowerHeuristic(settings.brdf_samples, pdf, settings.light_samples, lPdf_brdf);
						float falloff_factor = 1 / r_square;
						vec3 Li = light->intensity_multiplier * light->color * falloff_factor;
						light_contribution += brdf_sample * cosine_term * Li * path_throughput * weight / pdf;
						
						

					}
					if (pdf <= MY_FAV_MIN) // checking again here because don't want to return L while iterating above
						return L;

					if (length(light_contribution) != 0.f)
						L += light_contribution / (float)settings.brdf_samples;


					
					path_throughput = path_throughput * (brdf_sample * cosine_term) / pdf;
					if (myAll(path_throughput, MIN_THROUGHPUT)) return L;

				}
				/////////////////////////////////////////////
				// NO MIS 
				/////////////////////////////////////////////
				else {
					// sample all lights
					for (auto lightP : lights) {
						float sx, sz;
						concentricSampleDisk(&sx, &sz);
						{
							AreaLight* light = lightP.second;
							vec3 light_sample_point = light->position + light->radius * vec3(sx, 0.0f, sz); // this assumes positive y is world up
							vec3 dir = normalize(light_sample_point - hit_bias);

							Ray_t shadow_ray(hit_bias, dir);
							bool isLight = false;
							if (intersect(shadow_ray)) {
								isLight = getIntersection(shadow_ray).model->isLightSource;
							}
							if (isLight) {
								const float r_square = pow(length(light_sample_point - hit.position), 2.f);
								const float falloff_factor = 1.0f / r_square;
								vec3 Li = light->intensity_multiplier * light->color * falloff_factor;
								vec3 wi_l = normalize(light_sample_point - hit.position);
								float g = abs(dot(light->normal, -wi_l) * dot(hit.shading_normal, wi_l)) / r_square;
								float light_pdf = 1 / light->area;
								if(light_pdf >= MY_FAV_MIN)
									L += mat.f(wi_l, hit.wo, hit.shading_normal) * Li * g * path_throughput / light_pdf;
							}
						}
					}
					float pdf;
					vec3 sample = mat.sample_wi(wi, hit.wo, hit.shading_normal, pdf);
					if (pdf <= MY_FAV_MIN)
						return L;

					//float cosine_term = std::max(0.f, dot(wi, hit.shading_normal));
					float cosine_term = abs(dot(wi, hit.shading_normal));
					path_throughput = path_throughput * (sample * cosine_term) / pdf;
					if (myAll(path_throughput, MIN_THROUGHPUT)) return L;
				}


				///////////////////////////////
				// CREATE NEXT RAY ON PATH
				///////////////////////////////
				current_ray = Ray_t(hit_bias, wi);
				if (!intersect(current_ray)) {

					return L + path_throughput * Lenvironment(current_ray.ray.d);
				}
				
			}
			// we hit a light source so return L
			else 
				return L;

		}


		// Return the final outgoing radiance for the primary ray
		return L;
	}

	///////////////////////////////////////////////////////////////////////////
	// Used to homogenize points transformed with projection matrices
	///////////////////////////////////////////////////////////////////////////
	inline static glm::vec3 homogenize(const glm::vec4& p)
	{
		return glm::vec3(p * (1.f / p.w));
	}

	
	///////////////////////////////////////////////////////////////////////////
	// sample distance to object in front of camera to set focal distance
	///////////////////////////////////////////////////////////////////////////
	void traceSetFocal(const glm::mat4& V, const glm::vec3& dir) {
		vec3 camera_pos = vec3(glm::inverse(V) * vec4(0.0f, 0.0f, 0.0f, 1.0f));
		Ray_t sampleRay(camera_pos, dir);
		if (intersect(sampleRay)) {
			Intersection hit = getIntersection(sampleRay);
			float dist = length(hit.position - camera_pos);
			//dist = std::max(0.f, std::min(20.f, dist)); // clamp between 0 and 20.f
			settings.focald = dist;
		}
	}

	///////////////////////////////////////////////////////////////////////////
	// Trace one path per pixel and accumulate the result in an image
	///////////////////////////////////////////////////////////////////////////
	void tracePaths(const glm::mat4& V, const glm::mat4& P)
	{
		// Stop here if we have as many samples as we want
		if ((int(rendered_image.number_of_samples) > settings.max_paths_per_pixel)
			&& (settings.max_paths_per_pixel != 0))
		{
			return;
		}
		vec3 camera_pos = vec3(glm::inverse(V) * vec4(0.0f, 0.0f, 0.0f, 1.0f));

		// Trace one path per pixel (the omp parallel stuf magically distributes the
		// pathtracing on all cores of your CPU).
		int num_rays = 0;
		vector<vec4> local_image(rendered_image.width * rendered_image.height, vec4(0.0f));

		float focal_dist = settings.focald;
		float aperture_radius = settings.lensr;

#pragma omp parallel for
		for (int y = 0; y < rendered_image.height; y++)
		{
			for (int x = 0; x < rendered_image.width; x++)
			{
				vec3 color;
                
				// the current pixel on a virtual screen. [0,1]
				vec2 screenCoord = vec2(float(x) / float(rendered_image.width),
					float(y) / float(rendered_image.height));

				// Anti aliasing. Get random float, map to [-1,1], divide by width/height to get coordinate within pixel
				float jitterX = (randf() * 2.0f - 1.0f) / rendered_image.width;
				float jitterY = (randf() * 2.0f - 1.0f) / rendered_image.height;

				// Calculate direction, remaps to [-1,1]
				vec4 viewCoord = vec4((screenCoord.x + jitterX) * 2.0f - 1.0f, (screenCoord.y + jitterY) * 2.0f - 1.0f, 1.0f, 1.0f);
				vec3 p = homogenize(inverse(P * V) * viewCoord);

                Ray_t primaryRay(camera_pos, normalize(p - camera_pos));
				
				if (settings.thin_lens)
				{
					// Sample point on lense
					float angle = randf() * 2.0f * M_PI;
					float radius = sqrt(randf());
					vec2 offset = vec2(cos(angle), sin(angle)) * radius * aperture_radius;

					// Compute point on plane of focus
					vec3 local_ray = normalize(vec3(viewCoord.x, viewCoord.y, viewCoord.z));
					float ft = focal_dist / local_ray.z;
					vec3 pFocus = primaryRay.ray.o + primaryRay.ray.d * ft;

					// Update ray for effect of lense
					primaryRay.ray.o = vec3(glm::inverse(V) * vec4(offset.x, offset.y, 0.0f, 1.0f));;
					primaryRay.ray.d = normalize(pFocus - primaryRay.ray.o);
				}

                
                
				// Intersect ray with scene
				if (intersect(primaryRay))
				{
					// If it hit something, evaluate the radiance from that point
					color = Li(primaryRay);
				}
				else
				{
					// Otherwise evaluate environment
					color = Lenvironment(primaryRay.ray.d);
				}
				// Accumulate the obtained radiance to the pixels color
				float n = float(rendered_image.number_of_samples);
				rendered_image.data[y * rendered_image.width + x] =
					rendered_image.data[y * rendered_image.width + x] * (n / (n + 1.0f))
					+ (1.0f / (n + 1.0f)) * color;
			}
		}
		rendered_image.number_of_samples += 1;
	}
}; // namespace pathtracer

