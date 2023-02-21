## Task1
1. Change the absorption parameters to zero in scenes/volpath_test/volpath_test1.xml. What do you see? Why?
The result shows its orginal color of the emitted surface because no absorption attenuates the color.
2. In the homework, we assume the volume being not emissive. If you were tasked to modify the pseudo code above to add volume emission, how would you do it? Briefly describe your approach.
Using Monte Carlo to integral over the path and add up all the emission from particles. Multiple the particle_Le with the transmittance and integrate them over the (0 - t_hit).
#....some code before
camera_ray = sample_primary(camera, screen_pos, rng)
isect = intersect(scene, camera_ray)
u = next(rng) # u \in [0, 1]
t = -log(1 - u) / sigma_t
if t < isect.t_hit:
    trans_pdf = exp(-sigma_t * t) * sigma_t
    transmittance = exp(-sigma_t * t) * particle_Le
    p = camera_ray.org + t * camera_ray.dir
    # change here
    current_path_throughput *= transmittance / trans_pdf * particle_Le
else:
    # hit a surface, account for surface emission
    trans_pdf = exp(-sigma_t * isect.t_hit)
    transmittance = exp(-sigma_t * isect.t_hit)
    p = camera_ray.org + t_hit * camera_ray.dir
    current_path_throughput *= transmittance / trans_pdf 

## Task2
1. In the derivation above, how did we get from p(t) ∝ exp(−σtt) to p(t) = σt exp(σtt)?
- Because of the normalization of the transmittance exp(-sigma_t * t) to make the pdf. 
- Integrate from 0 to infinite over exp(-sigma_t * t) is 1/sigma_t. 
- exp(-sigma_t * t) / (1 / sigma_t) = sigma_t * exp(-sigma_t * t). 
2. How was Equation (13) incorporated into the pseudo code above? Why is it done this way? 
- The equation 13 is incorporated in the else branch when t is outside of [0, t_hit].
- In this situation, it hits the background radiance. The pdf of background radiance is exp(-sigma_t * t_hit).
Thus, the radiance returned is transmittance / pdf * Le. And for the reason why we do this, it's just because we use Monte Carlos 
to integrate the equation and the equation13 is exactly the last "light source term".
3. Play with the parameters σs and σa, how do they affect the final image? Why? (we assume monochro- matic volumes, so don’t set the parameters to be different for each channel.)
- If sigma_s is larger, the result is brighter for more in-scattering; If sigma_a is larger, the result is darker for more absorption.
4. Change the phase function from isotropic (the default one) to Henyey-Greenstein by uncommenting the phase function specification in the scene scenes/volpath_test/volpath_test2.xml. Play with the parameter g (the valid range is (−1,1)). What does g mean? How does the g parameter change the appearance? Why?
- If g is positve, it shows forward scattering and the result is brighter, else it shows backward scattering and the result is darker. 
- Reason: when g is positive, and when the cosine of the angle between omega and omega_prime becomes larger (which means the -ray.dir and dir_light are in opposite direction and a larger angle between them),
the phase_function gives a higher values; if g is negative, when their cosine is smaller (which means the -ray.dir and dir_light are in the some direction and the angle between them is smaller), the phase_function would 
also gives a higher values. In this way, it shows the property of forward/backward scattering.

## Task3
1. Play with the parameters σs, σa of different volumes, and change max_depth, how do they affect the final image? How does increasing/decreasing σs and σa of medium1 and medium2 affect the appearance, respectively? Why? Do different σs and σa values affect how high you should set max_depth? 
- Camera + left-top circle --> Medium1; right-bottom circle --> Medium2
- For medium1: If increase the sigma_s, the result turns to be smoother. The difference between in two circles and background is smaller: 
the background and the left-top are obviously darker. It's because because both the camera and left-top circle medium scatters more light; 
If increase the sigma_a, the overall result becomes even darker because both the camera and left-top circle medium absorbs more light.
- For medium2: If increase the sigma_s, the right-bottom circle is less transparent since more light scatter out and less light goes through it.
If increase the sigma_a, the right-bottom circle is darker since it absorbs more light.
- If max_depth is higher, the overall result looks brighter
- Different sigma_s and sigma_a should not affect the max_depth...they are two semantically different paramters (their meanings are different).
But assume that we want to get the similar result for different sigma_s/sigma_s/max_depth, if we increase the sigma_s, more rays would be scattered out and need more bounces (increase the max_depth) 
and hit the light. Also, if we increase the sigma_a, the transmittance is larger and we also need to increase the max_depth.
2. Switch to the Henyey-Greenstein phase function again. How does changing the g parameter affect the appearance? Why? 
- It stills change the proportion of backward/forward-scattering. 
- If medium1's g is higher, the overall result looks brighter and versus darker
- If medium2's g is higher, the right-bottom circle looks darker (not obvious change for the background and left-top circle)
3. Propose a phase function yourself (don’t have to describe the exact mathematical form). How would you design the shape of the phase function? What parameter would you set to control it?
Properties of phase function:
- From mathematics aspect, the phase function is able to be importance sampled with PDF. Thus it need to be anti-derivative.
- Need to be normalized over sphere otherwise it would add or substract radiance to a scattering.
- Need to be reciprocity like BSDF. 
- Isotropic medium should has the same possiblity to scatter ray in any direction (contains 1/4PI)??
- Parameters: at lease it should contain the dir_in and dir_out...and another parameter to control backward/forward-scattering....
reference: https://graphics.pixar.com/library/ProductionVolumeRendering/paper.pdf

## Task4
1. When will next event estimation be more efficient than phase function sampling? In our test scenes, which one is more efficient? Why?
- When the light source is small. The right one is more suitable since its light source is smaller. 
- Also, the "air" is vacumn for the right one, which means the ray has lower chance to scatter and hit the light.
2. In scenes/volpath_test/volpath_test4_2.xml, we render a scene with an object composed of dense volume. How does it compare to rendering the object directly with a Lambertian material? Why are they alike or different? 
- They are alike. For dense volume and diffuse surface, ray will scatter/reflect randomly.
3. Jim Kajiya famously has predicted in 1991 that in 10 years, all rendering will be volume rendering. What do you think that makes him think so? Why hasn’t it happened yet?
- Because it's more physically based...but still too complex and time-consuming even with nowadays hardwares.
- What's more, volumetric rendering is not artist-friendly (Hard to define a material with it)

## Task5
1. Play with the index of refraction parameter of the dielectric interface in scenes/volpath_test/volpath_test5_2.xml. How does that affect appearance? Why? 
- If the ior is higher, the inner sphere will be brighter since more light is refracted inside and vice versa.
2. In the scene scenes/volpath_test/vol_cbox_teapot.xml, we model the glass teapot as a transparent glass with blue homogeneous medium inside. What is the difference in terms of appearance between this approach and just making the color of the glass blue without any medium inside?
- If just making the color of the glass blue without any medium inside (any also without the roughness parameter), the result will look totally transparent. But after adding the volume, we can only see part of the other side (not totally transparent).

## Task6 
1. For heterogeneous volumes, what kind of distribution of the volume density makes the null scattering efficient/inefficient? Can you think of a way to improve our current sampling scheme in the inefficient case? 
- If the volume is too heterogeneous (which means the extinction coefficient sigma_t is not tight to the majorant sigma_m), it would increase the number of 
rejected interactions (null-collision). To optimize this situation, we can localize the majorant and bound it the extinction coefficient locally.
reference: https://cs.dartmouth.edu/wjarosz/publications/novak14residual.pdf
2. How do we make the null-scattering work for emissive volumes? Briefly describe a solution. 
- Add contribution of the particle emission when hitting a real particle. The contribution may be emssion * curreny_path_thoughput.
3. Why is it important to have an unbiased solution for volume rendering? Would it be sensible to have something that is biased but faster? How would you do it?
- Without an unbiased solution, it needs to sample a lot of times to converge to a realistic result.
- If want to have something that is biased but faster, it need a dedicated sampler as well as strong computating resources.
