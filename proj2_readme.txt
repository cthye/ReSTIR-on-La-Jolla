## Task1
1. The result shows its orginal color of the emitted surface because no absorption attenuates the color.
2. Using Monte Carlo to integral over the path and add up all the emission from particles.(??)

## Task2
1. Because of the normalization of the transmittance exp(-sigma_t * t) to make the pdf. 
- Integrate from 0 to infinite over exp(-sigma_t * t) is 1/sigma_t. 
- exp(-sigma_t * t) / (1 / sigma_t) = sigma_t * exp(-sigma_t * t). 
2. The equation 13 is incorporated in the else branch when t is outside of [0, t_hit].
In this situation, it hits the background radiance. The pdf of background radiance is exp(-sigma_t * t_hit).
Thus, the radiance returned is transmittance / pdf * Le.
3. If sigma_s is larger, the result is brighter for more in-scattering; If sigma_a is larger, the result is darker for more absorption.
4. If g is positve, it shows forward scattering and the result is brighter, else it shows backward scattering and the result is darker. (why?)

## Task3
1. Camera + left-top circle --> Medium1; right-bottom circle --> Medium2
- For medium1: If increase the sigma_s, the result turns to be smoother. The difference between in two circles and background is smaller: 
the background and the left-top are obviously darker. It's because because both the camera and left-top circle medium scatters more light; 
If increase the sigma_a, the overall result becomes even darker because both the camera and left-top circle medium absorbs more light.
- For medium2: If increase the sigma_s, the right-bottom circle is less transparent since more light scatter out and less light goes through it.
If increase the sigma_a, the right-bottom circle is darker since it absorbs more light.
- If max_depth is higher, the overall result looks brighter
- Different sigma_s and sigma_a should not affect the max_depth...(?)
2. It stills change the proportion of backward/forward-scattering. 
If medium1's g is higher, the overall result looks brighter and versus darker
If medium2's g is higher, the right-bottom circle looks darker (not obvious change for the background and left-top circle)
3. Properties of phase function:
- From mathematics aspect, the phase function is able to be importance sampled with PDF. Thus it need to be anti-derivative.
- Need to be normalized over sphere otherwise it would add or substract radiance to a scattering.
- Need to be reciprocity like BSDF. 
- Isotropic medium should has the same possiblity to scatter ray in any direction (contains 1/4PI)??
- Parameters: at lease it should contain the dir_in and dir_out...and another parameter to control backward/forward-scattering....
reference: https://graphics.pixar.com/library/ProductionVolumeRendering/paper.pdf

## Task4
1. When the light source is small. The right one is more suitable since its light source is smaller. 
Also, the "air" is vacumn for the right one, which means the ray has lower chance to scatter and hit the light.(?)
2. They are alike. For dense volume and diffuse surface, ray will scatter/reflect randomly.
3. More physically based...but still too complex and time-consuming even with nowadays hardwares.

## Task5
1. If the ior is higher, the inner sphere will be brighter since more light is refracted inside and vice versa.
2. If just making the color of the glass blue without any medium inside, the result will look too average (flat) and not so realistic without the scattering
behaviour inside the volume (????) 

## Task6 
1. If the volume is too heterogeneous (which means the extinction coefficient sigma_t is not tight to the majorant sigma_m), it would increase the number of 
rejected interactions (null-collision). To optimize this situation, we can localize the majorant and bound it the extinction coefficient locally.
reference: https://cs.dartmouth.edu/wjarosz/publications/novak14residual.pdf
2. Add contribution of the particle emission when hitting a particle. The contribution may be emssion * curreny_path_thoughput.(??)
3. Without an unbiased solution, it needs to sample a lot of times to converge to a realistic result.
If want to have something that is biased but faster, it need a dedicated sampler as well as strong computating resources.
