# Notes on Writing glTF PBR

These notes are my own personal summary on a rewrite of my PBR shaders to be a bit more maintainable, readable, and accurate. This rewrite was not possible without the following resources:

 * Appendix B of the glTF 2.0 Specification
 * glTF Sample Viewer Source
 * Enterprise PBR (cross referencing the above two references)
 * [Energy conserving IBL](https://bruop.github.io/ibl/#single_scattering_results)

## Core glTF 2.0 PBR (Metallic-Roughness)

The key starting point for this glTF PBR rewrite was to begin in Appendix B of the [glTF specification](https://registry.khronos.org/glTF/specs/2.0/glTF-2.0.html#appendix-b-brdf-implementation), implementing core glTF 2.0 closest to the intent of the Khronos PBR working group. Understand that this particular implementation aims to convey spirit of the material rather than promoting optimizations and PBR witchraft. Perfect for my use.

Rather than recite Appendix B back here, I'll simply fill in the blanks. Here are the nuts and bolts of my approach, any bugs are my own. We'll start with how punctual lighting (spot, point, direction lights) is handled on the PBR material.

I'll start at a high level implementation modeled after Appendix B, and drill down into each of these building blocks individually. First we'll need to calculate the specular and diffuse surface reflections, which will then be the core inputs into dielectric and metallic brdfs. The final color is selected based on the metallic value defined on the material.

```
    float specularBrdf = specular_brdf(alphaRoughness, NdotL, NdotV, NdotH, LdotH, VdotH);
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

float specular_brdf(float alphaRoughness, float NdotL, float NdotV, float NdotH, float LdotH, float VdotH)
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

## Adding Image-Based Lighting
The technique for Image-Based Lighting (IBL) is largely based off of the [SIGGRAPH 2013 course notes](https://blog.selfshadow.com/publications/s2013-shading-course/karis/s2013_pbs_epic_notes_v2.pdf) by Brian Karis at Epic Games. An approachable supplemental overview can also be found at [learnopengl.com](https://learnopengl.com/PBR/IBL/Diffuse-irradiance). There are several different alternative approaches to IBL available depending on your target hardware and runtime criteria. Emmett Lalish from Google developed a novel approach to IBL in three.js that is fast enough to be computed on the fly at runtime. Cesium.js also has an implementation that encodes the environment maps into an octahedron, which can be unrolled into a 2D texture (this technique is known as oct-encoding). Again, I just stuck with what I had originally written years ago and is well documented online. The approach I used is energy conserving thanks to the excellent summary presented by [Bruno Opsenica](https://bruop.github.io/ibl/#single_scattering_results).

First I'll present how surface light contributes are calculated for both the diffuse BRDF and the specular BRDF.

```
// IBL Irradiance represents the average lighting from any direction.
vec3 ibl_irradiance(vec3 diffuseColor, vec3 n, float roughness, float NdotV, vec3 F0, vec2 brdfLUT)
{
    n =  mat3(pbrInputs.environmentMapTransform) * n;
    // You're reading this right, the prefiltered diffuse IBL component 
    // is jammed into this miplevel intentionally. Save on samplers.
    vec3 diffuseLight = textureLod(u_SpecularEnvSampler, n, c_MaxLod).rgb;

    // The following models energy conservation for IBL
    // see https://bruop.github.io/ibl/
    // ss = single scattering, ms = multiple scattering
    vec3 Fr = max(vec3(1.0 - roughness), F0) - F0;
    vec3 k_S = F0 + Fr * pow(1.0 - NdotV, 5.0);
    vec3 FssEss = k_S * brdfLUT.x + brdfLUT.y;

    // Multiple scattering, from Fdez-Aguera
    float Ems = (1.0 - (brdfLUT.x + brdfLUT.y));
    vec3 F_avg = F0 + (1.0 - F0) / 21.0;
    vec3 FmsEms = Ems * FssEss * F_avg / (1.0 - F_avg * Ems);
    vec3 k_D = diffuseColor * (1.0 - FssEss + FmsEms);

    return (FmsEms + k_D) * diffuseLight;
}

vec3 ibl_specular(float roughness, vec3 n, vec3 v, vec3 F0, vec2 brdfLUT)
{
    const float roughnessOneLOD = c_MaxLod - 1.0; // c_MaxLod contains the diffuseIBL
    float lod = roughnessOneLOD * roughness * (2.0 - roughness);
    vec3 reflection = normalize(reflect(-v, n));
    reflection = mat3(pbrInputs.environmentMapTransform) * reflection;
    reflection = normalize(reflection);
    vec3 specularLight = textureLod(u_SpecularEnvSampler, reflection, lod).rgb;

    // The following models energy conservation for IBL
    // see https://bruop.github.io/ibl/
    // ss = single scattering, ms = multiple scattering
    vec3 Fr = max(vec3(1.0 - roughness), F0) - F0;
    vec3 k_S = F0 + Fr * pow(1.0 - dot(n,v), 5.0);
    vec3 FssEss = k_S * brdfLUT.x + brdfLUT.y;

    return specularLight * FssEss;
}
```

Below I show how these terms are combined with the results of punctual lighting covered earlier in this document.

```
    vec2 uv = clamp(vec2(abs(NdotV), 1.0 - perceptualRoughness), vec2(0.0, 0.0), vec2(1.0, 1.0));
    vec2 brdfLUT = texture(u_brdfLUT, uv).rg;

    float specularBrdf = specular_brdf(alphaRoughness, NdotL, NdotV, NdotH, LdotH, VdotH);
    vec3 f_specular = NoL * vec3(specularBrdf);
    vec3 f_diffuse = NoL * diffuse_brdf(baseColor.rgb);

    vec3 iblIrradiance = ibl_irradiance(
        diffuseColor, 
        n, 
        perceptualRoughness, 
        NoV, 
        reflectanceF0, 
        brdfLUT);

    vec3 iblSpecularColor = ibl_specular(
        perceptualRoughness, 
        n, 
        v, 
        reflectanceF0, 
        brdfLUT);

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

    dielectric_brdf += iblIrradiance + iblSpecularColor;
    metal_brdf += iblIrradiance + iblSpecularColor;
    vec3 color = mix(dielectric_brdf, metal_brdf, metallic);
```

As you can see, the approach I took was to combine irradiance (average light from all directions) and specular IBL contributions with the dielectric and metallic BRDFs at the end.

## Extending Metallic-Roughness for Transmission and Refractive Volumes

