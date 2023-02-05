# Notes on Writing glTF PBR

These notes are my own personal summary on a rewrite of my PBR shaders to be a bit more maintainable, readable, and accurate. My approach was to follow the definition of PBR as defined in Appendix B of the glTF 2.0 specification as closely as possible. Sample code is provided below in GLSL.

TODO:
 * Summary of creating an offscreen scene for transmission materials
 * Extending thin-walled transmission to handle refractive volumes
 * Include diagrams from Appendix B and Extensions
 * Include Renderings of the materials.

## Core glTF 2.0 PBR (Metallic-Roughness)

The key starting point for this glTF PBR rewrite was to begin in Appendix B of the [glTF specification](https://registry.khronos.org/glTF/specs/2.0/glTF-2.0.html#appendix-b-brdf-implementation), implementing core glTF 2.0 closest to the intent of the Khronos PBR working group. Understand that this particular implementation aims to convey spirit of the material rather than promoting optimizations and PBR witchraft. Perfect for my use.

Rather than recite Appendix B back here, I'll simply fill in the blanks. Here are the nuts and bolts of my approach, any bugs are my own. We'll start with how punctual lighting (spot, point, direction lights) is handled on the PBR material.

I'll start at a high level implementation modeled after Appendix B, and drill down into each of these building blocks individually. First we'll need to calculate the specular and diffuse surface reflections, which will then be the core inputs into dielectric and metallic brdfs. The final color is selected based on the metallic value defined on the material.

```GLSL
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
```GLSL
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

```GLSL
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

```GLSL
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

```GLSL
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

```GLSL
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

## Extending Metallic-Roughness for Transmission

For modeling glass like surfaces, we need to extend the core metallic roughness model to include transmission. The [KHR_materials_transmission](https://github.com/KhronosGroup/glTF/tree/main/extensions/2.0/Khronos/KHR_materials_transmission) extension allows for light to transmit through the surface. That might sound like alpha blending to some, but these two concepts are very different. In short, alpha blending is used to model the precense or absence of the material (think screen door). 50% alpha implies that only 50% of the material is actually there in the fragment. If that surface was shiny, you would only render 50% of that shiny surface and lose half of the highlight. For transmission, the entire material is present but we now signal to the renderer that light travels through the surface and we must therefore combine surface highlights, tint, etc on this material with the scene that exists behind this fragment. This is covered in our [Khronos webinar](https://www.khronos.org/events/advanced-pbr-material-parameters-in-gltf) and in the specification.

### Preparing the offscreen opaque scene for Transmission

TODO. See the second half of our [Khronos webinar](https://www.khronos.org/events/advanced-pbr-material-parameters-in-gltf) for the overview.

### Shader updates for transmission.

We'll start by defining our transmission response with punctual lighting. Transmission extends the core glTF model by inserting a new node, the `specular_btdf` (Bi-directional Transmission Function). At its core, the glTF model will provide a `transmissionFactor` that we will use to blend between the `diffuse_brdf` and the `specular_btdf`.

```GLSL
float heaviside(float v)
{
    // return 1.0 if v is > 0, 0.0 otherwise
    return clamp(sign(v), 0.0, 1.0);
}

float specular_btdf(float alphaRoughness, vec3 n, vec3 l, vec3 v, float ior)
{
    l = l + 2.0 * n * (dot(-l, n)); // mirror light reflection vector on surface
    vec3 Ht = normalize(v + l);

    float D_heaviside = heaviside(clampdot(n,Ht));
    float G_heaviside = heaviside(clampdot(Ht,v) / clampdot2(n,v));
    
    float Dt = D_heaviside * microfacetDistribution(alphaRoughness, clampdot(n,Ht));
    float Gt = G_heaviside * geometricOcclusion(alphaRoughness, clampdot(n,l), clampdot(n,v), clampdot(l,Ht), clampdot(v,Ht));

    float Vt = Gt / (4.0 * abs(dot(n,l)) * abs(dot(n,v))); 

    return clampdot(n, l) * Vt * Dt;
}
```

The specular_btdf should look very familiar, this is essentially the specular_brdf function, but we've modified it to be calculated with the light vector mirrored on the surface, and flipped our hemisphere of acceptable values (via the `heaviside`). The `microfacetDistribution` and `geometricOcclusion` are the same functions defined above. You may be tempted to simplify things and send different inputs into the specular_brdf function, but that unused `float ior` input is foreshadowing for changes that we will make to extend this effect to include refraction. 

To factor the `specular_btdf` into our material, we will need to mix the resulting value with the `diffuse_brdf`.

```GLSL
    float transmissionFactor = pbrInputs.transmissionFactor;
    transmissionFactor *= texture(u_TransmissionSampler, texCoords).r;
    vec3 f_transmission = specular_btdf(alphaRoughness, n, l, v, pbrInputs.ior) * baseColor.rgb;
    f_diffuse = mix(f_diffuse, f_transmission, transmissionFactor);
```


No surprises here. In my implementation, I will bind a 1x1 white texture to the u_TransmissionSampler descriptor set if there is no transmission texture available. This is merely to reduce the number of shader permutations required. For transmission surfaces, we also require a slightly modified `fresnel_mix` to use the modified half-vector.

```GLSL
    vec3 fresnel_mix(float ior, vec3 base, vec3 layer, vec3 n, vec3 l, vec3 v)
    {
        l = l + 2.0 * n * (dot(-l, n)); // mirror light reflection vector on surface
        vec3 Ht = normalize(v + l);

        float f0 = ((1.0 - ior) / (1.0 + ior)) * ((1.0 - ior) / (1.0 + ior));

        float invVoH = 1.0 - abs(clampdot(v, Ht));
        float pow5 = invVoH * invVoH * invVoH * invVoH * invVoH;
        float fr = f0 + (1.0 - f0) * pow5;

        return mix(base, layer, fr);
    }
```

At a minimum, rasterizers like this implementation must "reveal" the opaque scene through a transmission surface. More advanced renders may take this requirement further by rendering layers upon layers of transmission materials to allow for stacking of the effect, though it is really hard to get this right. It turns out that modeling a mug of beer is really hard!. 

Below is how I sample into the offscreen scene of opaque objects. This scene includes the Skybox and all solid objects, rendered to a 1024x1024 UNORM texture. Emphasis added to the UNORM, as this buffer isn't here to represent colors, but linear light values from this scene. As this scene is a mipmapped texture, each level in the mipmap tree contains a lower resolution and blurrier scene from the level above. As roughness increases, we want to sample into these lower levels. Note that this LOD selection is also influenced by the IOR of the material to retain two key properties; (1) as IOR approaches 1.0 sensitivity to roughness matters less (1.0 is air and we will always sample miplevel 0), and (2) a default IOR of 1.5 represents our identity value for IOR influence (no skew).

```GLSL
vec3 ibl_transmission(float roughness, float ior, vec3 absorptionColor, vec3 v, vec3 n)
{
    const float maxLod = 6; // scene is now always 1024x1024 with mipLevels = 10
    const float roughnessOneLOD = maxLod - 1.0;
    // maps 1.0 ior to LOD 0, and 1.5 ior (default) to put IOR influence at 1.0
    const float iorLODInfluence = (ior - 1.0) * 2.0;
    float lod = iorLODInfluence * roughnessOneLOD * roughness * (2.0 - roughness);
    vec2 uv = gl_FragCoord.xy / scene.viewport.zw;

    vec3 transmittedLight = textureLod(u_TransmissionScene, uv, lod).rgb;
    vec3 transmittanceColor = transmittedLight * absorptionColor;
    return transmittanceColor;
}
```

I understand that calling this function "ibl" is a stretch, as it has nothing to do with our IBL cubemap. We are sampling an image for lighting though, so it stays. You'll see below that it just fits in nicely with our other IBL sampling. Note that as covered above for the punctual lighting case, iblTransmission is not a component of the metal_brdf.

```GLSL
    vec3 iblTransmission = transmissionFactor * ibl_transmission(
        perceptualRoughness, 
        pbrInputs.ior, 
        baseColor.rgb, 
        v,
        n);

    vec3 dielectricIblColor = mix(iblIrradiance, iblTransmission, transmissionFactor);
    dielectric_brdf += dielectricIblColor + iblSpecularColor;
    metal_brdf += iblIrradiance + iblSpecularColor;
    vec3 color = mix(dielectric_brdf, metal_brdf, metallic);
```

## Extending Transmission to include Refractive Volumes
