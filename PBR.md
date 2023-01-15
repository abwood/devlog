# Notes on Writing glTF PBR

Over a recent holiday break, I decided to rewrite my PBR shaders to be a bit more maintainable, readable, and accurate. The previous iteration was rooted in the early days of glTF 2.0, based off the original glTF Reference Viewer ... quickly renamed to glTF Sample Viewer as it was clear that the old implementation was good but should never be considered an authoritative "reference". This rewrite was not possible without the following resources:

 * Appendix B of the glTF 2.0 Specification
 * glTF Sample Viewer Source
 * Enterprise PBR (cross referencing the above two links)

## Core glTF 2.0 PBR (Metallic-Roughness)

The key starting point for this glTF PBR rewrite was to begin in Appendix B of the [glTF specification](https://registry.khronos.org/glTF/specs/2.0/glTF-2.0.html#appendix-b-brdf-implementation), implementing core glTF 2.0 closest to the intent of the Khronos PBR working group. Understand that this particular implementation aims to convey spirit of the material rather than promoting optimizations and PBR witchraft. Perfect for my use.

Rather than recite Appendix B back here, I'll simply fill in the blanks. Here are the nuts and bolts of my approach, any bugs are my own. We'll start with how punctual lighting (spot, point, direction lights) is handled on the PBR material and then weave in IBL.

I'll start at a high level implementation modeled after Appendix B, and drill down into each of these building blocks individually. First we'll need to calculate the specular and diffuse surface reflections, which will then be the core inputs into dielectric and metallic brdfs. The final color is selected based on the metallic value defined on the material.

```
    float specularBrdf = specular_brdf(alphaRoughness, NdotH, NdotL, NdotV, LdotH, VdotH);
    vec3 f_specular = NdotL * vec3(specularBrdf);
    vec3 f_diffuse = NdotL * diffuse_brdf(baseColor.rgb);

    vec3 dielectric_brdf = 
        fresnel_mix(
            pbrInputs.ior,
            f_diffuse,          // base
            f_specular,         // layer
            VdotH);

    vec3 metal_brdf = 
        conductor_fresnel(
            baseColor.rgb,      // f0
            f_specular,         // bsdf
            VdotH);

    vec3 color = mix(dielectric_brdf, metal_brdf, metallic);
```

Note the addition of NdotL factor to f_specular and f_diffuse above, which ensures that these terms are only visible on the surfaces that face the light source. Special care must be taken to ensure that these dot products are clamped to 0.0, or sometimes to a small 0.0001 for cases when the dot product is used in the denominator.

Diffuse and specular brdfs follow exactly as defined in Appendix B.
```
vec3 diffuse_brdf(vec3 color)
{
    return (1.0 / c_Pi) * color;
}

float specular_brdf(float alphaRoughness, float NdotL, float NdotV, float LdotH, float VdotH)
{
    float G = geometricOcclusion(alphaRoughness, NdotL, NdotV, LdotH, VdotH);
    float V = G / (4.0 * abs(NdotL) * abs(NdotV));
    float D = microfacetDistribution(alphaRoughness, NdotH);

    return V * D;
}
```

The geometricOcclusion (G) and microfacetDistribution (D) functions are defined in the SIGGRAPH 2013 course notes on PBR by Epic Games.

```
// The following equation(s) model the distribution of microfacet normals across the area being drawn (aka D())
// Implementation from "Average Irregularity Representation of a Roughened Surface for Ray Reflection" by T. S. Trowbridge, and K. P. Reitz
// Follows the distribution function recommended in the SIGGRAPH 2013 course notes from EPIC Games [1], Equation 3.
float microfacetDistribution(float alphaRoughness, float NdotH)
{
    float alphaRoughnessSq = alphaRoughness * alphaRoughness;
    float f = NdotH * NdotH * (alphaRoughnessSq - 1.0) + 1.0;

    float hs = heaviside(NdotH);
    return hs * alphaRoughnessSq / (c_Pi * f * f);
}

// This calculates the specular geometric attenuation (aka G())
// where rougher material will reflect less light back to the viewer.
// This implementation is based on [1] Equation 4, and we adopt their modifications to
// alphaRoughness as input as originally proposed in [2].
// Note that G has been confused with Vis.  Sample Viewer uses Vis = (G / (4 * NdotL * NdotV))
// If this ever gets switched to Vis, be sure to drop the 4 * NdotL * NdotV normalization term
float geometricOcclusion(float alphaRoughness, float NdotL, float NdotV, float LdotH, float VdotH)
{
    float alphaRoughnessSq = alphaRoughness * alphaRoughness;
    // Positive characteristic function: one if a > 0 and zero if a <= 0
    float hs_l = heaviside(LdotH);
    float hs_v = heaviside(VdotH);

     float attenuationL = 2.0 * abs(NdotL) / (abs(NdotL) + sqrt(alphaRoughnessSq + (1.0 - alphaRoughnessSq) * (NdotL * NdotL)));
     float attenuationV = 2.0 * abs(NdotV) / (abs(NdotV) + sqrt(alphaRoughnessSq + (1.0 - alphaRoughnessSq) * (NdotV * NdotV)));
     return hs_l * hs_v * attenuationL * attenuationV;
}
```

Finally, here are my implementations of `conductor_fresnel` and `fresnel_mix`:

```
vec3 conductor_fresnel(vec3 f0, vec3 bsdf, float VdotH)
{
    vec3 f90 = vec3(1.0);
    float invVoH = (clamp(1.0 - VdotH, 0.0, 1.0)); // 1.0 - abs(VdotH);
    float pow5 = invVoH * invVoH * invVoH * invVoH * invVoH;
    return bsdf * (f0 + (f90 - f0) * pow5);
}

vec3 fresnel_mix(float ior, vec3 base, vec3 layer, float VdotH)
{
    float f0 = ((1.0 - ior) / (1.0 + ior)) * ((1.0 - ior) / (1.0 + ior));

    float invVoH = 1.0 - abs(VdotH);
    float pow5 = invVoH * invVoH * invVoH * invVoH * invVoH;
    float fr = f0 + (1.0 - f0) * pow5;

    return mix(base, layer, fr);
}
```
