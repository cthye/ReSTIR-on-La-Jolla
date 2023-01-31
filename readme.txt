## Diffuse
1. Compare the two BRDFs with a Lambertian BRDF (try to render images with all three BRDFs):
what differences do you see? Why?
- For rough material, subsurface BRDF & base diffuse BRDF are both brighter at the edge and stronger grazing retroreflection compared with Lambertian BRDF.
- They both has higher F_D/F_SS term at grazing angle for diffuse materials, which produce a retroreflection peak at grazing angle.

2. Play with the roughness parameter: how does it affect the appearance?
- For smooth material (roughness is small), the edge is darken. 
- The f is attenuated by the lower roughness at grazing angle.

3. Compare the base diffuse BRDF (fbaseDiffuse) with the subsurface BRDF (fsubsurface) by playing with the subsurface parameter. What differences do you see? Why? In what lighting condition does the base diffuse BRDF differ the most from the subsurface BRDF? (Play with the light position in simple_sphere.xml for your experimentation)
- It becomes brighter if subsurface parameter is larger.
- In subsurface model has a higher fresnel peak and make the model look more transparent and brighter.
- When the light is in the opposite direction of view, the difference is more obvious.

4. (Optional,bonus3%) Another popular option for modeling diffuse surfaces is the Oren-Nayar BRDF[9], which is used in the Autodesk’s Standard Surface BSDF. What is the difference between the Oren-Nayar BRDF and the Disney diffuse BRDF? What are the pros and cons? What is your preference?
- The Oren-Nayar fails to get a fall at grazing retroreflection for smooth materials compared with DisneyBSDF.
- Pros: Oren-Nayer behaves well for rough materials (still get a grazing retro-reflective peak) with simpler computation.
- Cons: Not behaves well for smooth materials.
- I'd prefer the DisneyBSDF for more realistic looking results and the computation overhead is not that big.

## Metal 
1. Compare DisneyMetal with the roughplastic material (try to render images with both BRDFs). What
differences do you see? Why?
- DisneyMetal is anisotropic, having a stronger concentration and longer tail highlihgt.
- DisenyMetal has a NDF with longer tails as well as considering anisotropic parameter (a ellipsoid distribution of normals).

2. Change the roughness parameters. Apart from how specular the surface it, do you observe any other differences?
- It looks more like isotropic as roughness grows.

3. A popular alternative over the Trowbridge-Reitz normal distribution function is the Beckmann distri- bution (a Gaussian distribution on the slopes of the normals). What are the differences
between Trowbridge-Reitz and Beckmann? Why did Disney folks choose to use Trowbridge-Reitz instead of Beckmann? (You might want to read the awesome article Slope Space in BRDF Theory from Nathan Reed.)
- GGX NDF is derived from multivariate student-t distribution with v = 2. When v -> unlimited, it becomes the Beckmann Guassian.
- The difference lies on the tail and peak: GGX has a longer tail and falls more slowerly than Beckmann.
- GGX looks more realistic.

4. (Optional, bonus 3%) What are the pros and cons of the Schlick approximation compared to the actual Fresnel equation? What is your preference? (You may want to read/watch the Naty Hoffman presentation in the footnote above.)
- Pros: computation is simpler; in the context of RGB, its approximation is more principled using quantities both physically and perceptually meaningful; provides artistic edge lighting control
- Cons: not physically accurate in spectral specturm.

## Clearcoat 
1. Compare DisneyClearcoat with DisneyMetal using similar roughness. What differences do you see? Where do the differences come from?
- Clearcoat is isotropic: only single isotropic roughness is considered
- Clearcoat is more blur: its DNF has longer tail

2. For coating, Autodesk Standard Surface uses a standard Trowbridge-Reitz microfacet distribution, instead of the modified normal distribution function proposed by Burley. What are the pros and cons? What is your preference?
- Pros: can be anisotropic; has geometric meaning; has a longer tail distribution
- Cons: more computation
- I'd prefer standard Trowbridge-Reitz microfacet distribution because hard-coded ior is confusing and may not be accurate all the time.

3. Why do Burley limit the clearcoat BRDF to be isotropic?
- The metal lobe is already anisotropic and metallic. Therefore, clearcoat lobe is always isotropic and non-metallic. 
- Also, being isotropic simpifies the computation.

## Glass 
1. Why do we take a square root of baseColor in the refractive case?
- to account for the two refractions when going through the material 

2. Play with the index of refraction parameter η (the physically plausible range is [1,2]). How does it affect appearance?
- It appears a stronger and wider peak when η is larger.
- Also the reflection is stronger when η is larger and the outer layer looks transparent when η is near to 1.

3. If a refractive object is not a closed object and we assign it to be a glass BSDF, will everything still work properly? Why? Propose a solution if it will not work properly (hint: you may want to read the thin-surface BSDF in Burley’s note [2].).
- The light transmitted into object will escape from the open boundary (and darken the object?). Using thin-surface BSDF can simulate both the entrance and exit scattering events at a single shading point and get a appearance for a refractive object as if it's closed. (???

4. (Optional, 3%) Replace the dielectric Fresnel equation with a Schlick approximation (see Burley’s course notes [2] on the fix to the Schlick approximation to make it work for η < 1). Do you observe any differences when η = 1.5? What about η = 1.01?
- for eta = 1.5, not obvious difference (though Schlick approximation's result seems darker than Fresnel equation's)
- for eta = 1.01, the Schlick approximation's result is brighter than the other's.

## Sheen
1. Render the simple_sphere scene with the sheen BRDF. What do you see? Why? What happens if you
change the position of the light source?
- If the light comes from the view direction, it shows nothing...
- If the light direction is the vertical/opposite direction of view direction, it shows tinted specular at grazing angle because more light is transmitted through the object.

2. Play with the parameter sheenTint, how does it affect the appearance? Why?
- It looks more colorful as sheenTint is larger. This parameter decides the achromatism (the weight of base color). 

3. In Autodesk Standard Surface, the sheen is modeled by a microfacet BRDF [4]. What are the pros and cons between the Autodesk approach and the Disney approach? What is your preference?
- Pros: Microfacet BRDF based sheen looks more physically realistic and softer; stimulates backward scattering.
- Cons: More complex
- I'd prefer the autodesk's since it's more physically realisitc and accurate

## Put together
1. What are the differences between the specular and metallic parameters? How do they affect the
appearance?
- Specular describes how the dielectric layer (the one upon the diffuse layer in diffuse lobe) reflects while metallic describes how it looks like a metal (at least their semestic meaning are different)
- If the specular parameter is larger, the object has a wider specular reflection (looks like a rough plastic)
- If the metallic parameter is larger, the object has sharper and stronger specular

2. What are the differences between the roughness and clearcoat_gloss parameters? How do they affect the appearance?
- Clear coat gloss describe how the object looks reflective with a given roughness (despite how rough it's, it always shows reflective like a dielectric material) 
- Roughness shows how specular or diffuse is the object.
....???

3. Play with the specularTint parameter. How does it affect the appearance?
- It tints the incident specular towards the base color. If the specularTint is larger, the specular is more colorful (like the base color)

4. The roughness parameter affects many components of the BSDFs at once (e.g., both the diffuse and metal BRDF use the roughenss parameter). How do you feel about this? If you are an artist using the Disney BSDF, would you want to have a separate roughness parameter for each component?
- Hard to understand and control the visual affect. Yes.

## The final task
Please check with the scenes in side the "proj1" directory. Besides the test scene, I also rendered other scenes made in Blender.