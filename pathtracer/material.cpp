#include "material.h"
#include "sampling.h"

namespace pathtracer
{
///////////////////////////////////////////////////////////////////////////
// A Lambertian (diffuse) material
///////////////////////////////////////////////////////////////////////////
vec3 Diffuse::f(const vec3& wi, const vec3& wo, const vec3& n)
{
	if(dot(wi, n) <= 0.0f) 
		return vec3(0.0f);
	if(!sameHemisphere(wi, wo, n))
		return vec3(0.0f);
	return (1.0f / M_PI) * color;
}

vec3 Diffuse::sample_wi(vec3& wi, const vec3& wo, const vec3& n, float& p)
{
	vec3 tangent = normalize(perpendicular(n));
	vec3 bitangent = normalize(cross(tangent, n));
	vec3 sample = cosineSampleHemisphere();
	wi = normalize(sample.x * tangent + sample.y * bitangent + sample.z * n);
	if(dot(wi, n) <= 0.0f)
		p = 0.0f;
	else
		p = max(0.0f, dot(n, wi)) / M_PI;
	return f(wi, wo, n);
}

float Diffuse::get_pdf(const vec3& wi, const vec3& wo, const vec3& n, float& p)
{
	if (dot(wo, n) <= 0.0f) return 0.0f;
	if (dot(wi, n) <= 0.0f)
		p = 0.0f;
	else
		p = max(0.0f, dot(n, wi)) / M_PI;
	return 0.f;
}
///////////////////////////////////////////////////////////////////////////
// A Blinn Phong Dielectric Microfacet BRFD
///////////////////////////////////////////////////////////////////////////
vec3 BlinnPhong::refraction_brdf(const vec3& wi, const vec3& wo, const vec3& n)
{
	

	if (refraction_layer == NULL) {
		return vec3(0.0f);
	}
	else {
		vec3 wh = normalize(wi + wo);
		float fwi = R0 + (1.0f - R0) * pow(1.0f - dot(wh, wi), 5.0f);
		return (1.0f - fwi) * refraction_layer->f(wi, wo, n);
	}
	
}
vec3 BlinnPhong::reflection_brdf(const vec3& wi, const vec3& wo, const vec3& n)
{

	float won = dot(wo, n);
	float win = dot(wi, n);
	

	if (won < 0 || win < 0 || !sameHemisphere(wi, wo, n)) return vec3(0.0f);


	vec3 wh = normalize(wi + wo);


	float fwi = R0 + (1.0f - R0) * pow(1.0f - dot(wh, wi), 5.0f);

	// R0 already defined 
	float s = shininess;

	float whn = dot(n, wh);
	float wowh = dot(wo, wh);

	float dwh = ((s + 2) / (2 * M_PI)) * pow(whn, s);
	


	float t1 = 2.0f * (whn * won) / wowh;
	float t2 = 2.0f * (whn * win) / wowh;
	float gwiwo = min(1.0f, min(t1, t2));

	float denom = 4.0f * won * win;

	if (denom <= 0.00005f) return vec3(0.0f);

	float brdf = (fwi * dwh * gwiwo) / denom;

	return vec3(brdf);
}

vec3 BlinnPhong::f(const vec3& wi, const vec3& wo, const vec3& n)
{
	if (dot(wi, n) <= 0.0f || !sameHemisphere(wi, wo, n)) return vec3(0.0f);

	return reflection_brdf(wi, wo, n) + refraction_brdf(wi, wo, n);

	
}



float BlinnPhong::get_pdf(const vec3& wi, const vec3& wo, const vec3& n, float& p) {

	vec3 wh = normalize(wi + wo);

	if (dot(wo, n) <= 0.0f) return 0.0f;
	
	if (randf() < 0.5f) {
		float pwh = (shininess + 1) * pow(dot(n, wh), shininess) / (2.0f * M_PI);
		p = pwh / (4.0f * dot(wo, wh)); 
		p *= 0.5f;
	}
	else {
		if (refraction_layer == NULL) return 0.f;
		refraction_layer->get_pdf(wi, wo, n, p);
		p *= 0.5f;
	}
	
	return 0.f;
	
}




vec3 BlinnPhong::sample_wi(vec3& wi, const vec3& wo, const vec3& n, float& p)
{

	if (dot(wo, n) <= 0.0f) return vec3(0.0f);
	vec3 tangent = normalize(perpendicular(n));
	vec3 bitangent = normalize(cross(tangent, n));
	float phi = 2.0f * M_PI * randf();
	float cos_theta = pow(randf(), 1.0f / (shininess + 1));
	float sin_theta = sqrt(max(0.0f, 1.0f - cos_theta * cos_theta));
	vec3 wh = normalize(sin_theta * cos(phi) * tangent +
		sin_theta * sin(phi) * bitangent +
		cos_theta * n);
	
	
	if (randf() < 0.5f) {
		float pwh = (shininess + 1) * pow(dot(n, wh), shininess) / (2.0f * M_PI);
		wi = reflect(-wo, wh);			// forgot to negate first, big confusion
		p = pwh / (4.0f * dot(wo, wh));
		p *= 0.5f;
		
		return reflection_brdf(wi, wo, n);
	}
	else {
		if (refraction_layer == NULL) return vec3(0.0f);
		vec3 brdf = refraction_layer->sample_wi(wi, wo, n, p);
		p *= 0.5f;
		float f = R0 + ((1.0f - R0) * pow(1.0f - abs(dot(wh, wi)), 5.0f));
		return (1 - f) * brdf;
		

	}
	

}

///////////////////////////////////////////////////////////////////////////
// A Blinn Phong Metal Microfacet BRFD (extends the BlinnPhong class)
///////////////////////////////////////////////////////////////////////////
vec3 BlinnPhongMetal::refraction_brdf(const vec3& wi, const vec3& wo, const vec3& n)
{
	return vec3(0.0f);
}
vec3 BlinnPhongMetal::reflection_brdf(const vec3& wi, const vec3& wo, const vec3& n)
{
	return BlinnPhong::reflection_brdf(wi, wo, n) * color;
};

///////////////////////////////////////////////////////////////////////////
// A Linear Blend between two BRDFs
///////////////////////////////////////////////////////////////////////////
vec3 LinearBlend::f(const vec3& wi, const vec3& wo, const vec3& n)
{
	// using weight, should be right
	if (dot(wi, n) <= 0.0f || !sameHemisphere(wi, wo, n)) return vec3(0.0f);

	return w * bsdf0->f(wi, wo, n) + (1-w) * bsdf1->f(wi, wo, n);
}

vec3 LinearBlend::sample_wi(vec3& wi, const vec3& wo, const vec3& n, float& p)
{

	if (randf() < w) {
		vec3 ret = bsdf0->sample_wi(wi, wo, n, p);
		p *= w; 
		return ret;
	}
	else {
		vec3 ret = bsdf1->sample_wi(wi, wo, n, p);
		p *= (1-w);
		return ret;
	}


}

float LinearBlend::get_pdf(const vec3& wi, const vec3& wo, const vec3& n, float& p)
{
	
	if (randf() < w) {
		float ret = bsdf0->get_pdf(wi, wo, n, p);
		p *= w; 
		return ret;
	}
	else {
		float ret = bsdf1->get_pdf(wi, wo, n, p);
		p *= (1-w);
		return ret;
	}
}


///////////////////////////////////////////////////////////////////////////
// A perfect specular refraction.
///////////////////////////////////////////////////////////////////////////
} // namespace pathtracer