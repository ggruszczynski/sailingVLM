import numba
import numpy as np
import time

from numba.typed import List

from numpy.linalg import norm


# collocation points
collocation_points = np.array([[  0.1497944302, -90.          ,  -0.0078503934],
       [  0.1497944302, -70.          ,  -0.0078503934],
       [  0.1497944302, -50.          ,  -0.0078503934],
       [  0.1497944302, -30.          ,  -0.0078503934],
       [  0.1497944302, -10.          ,  -0.0078503934],
       [  0.1497944302,  10.          ,  -0.0078503934],
       [  0.1497944302,  30.          ,  -0.0078503934],
       [  0.1497944302,  50.          ,  -0.0078503934],
       [  0.1497944302,  70.          ,  -0.0078503934],
       [  0.1497944302,  90.          ,  -0.0078503934],
       [  0.3495203372, -90.          ,  -0.0183175847],
       [  0.3495203372, -70.          ,  -0.0183175847],
       [  0.3495203372, -50.          ,  -0.0183175847],
       [  0.3495203372, -30.          ,  -0.0183175847],
       [  0.3495203372, -10.          ,  -0.0183175847],
       [  0.3495203372,  10.          ,  -0.0183175847],
       [  0.3495203372,  30.          ,  -0.0183175847],
       [  0.3495203372,  50.          ,  -0.0183175847],
       [  0.3495203372,  70.          ,  -0.0183175847],
       [  0.3495203372,  90.          ,  -0.0183175847],
       [  0.5492462441, -90.          ,  -0.0287847759],
       [  0.5492462441, -70.          ,  -0.0287847759],
       [  0.5492462441, -50.          ,  -0.0287847759],
       [  0.5492462441, -30.          ,  -0.0287847759],
       [  0.5492462441, -10.          ,  -0.0287847759],
       [  0.5492462441,  10.          ,  -0.0287847759],
       [  0.5492462441,  30.          ,  -0.0287847759],
       [  0.5492462441,  50.          ,  -0.0287847759],
       [  0.5492462441,  70.          ,  -0.0287847759],
       [  0.5492462441,  90.          ,  -0.0287847759],
       [  0.7489721511, -90.          ,  -0.0392519672],
       [  0.7489721511, -70.          ,  -0.0392519672],
       [  0.7489721511, -50.          ,  -0.0392519672],
       [  0.7489721511, -30.          ,  -0.0392519672],
       [  0.7489721511, -10.          ,  -0.0392519672],
       [  0.7489721511,  10.          ,  -0.0392519672],
       [  0.7489721511,  30.          ,  -0.0392519672],
       [  0.7489721511,  50.          ,  -0.0392519672],
       [  0.7489721511,  70.          ,  -0.0392519672],
       [  0.7489721511,  90.          ,  -0.0392519672],
       [  0.948698058 , -90.          ,  -0.0497191584],
       [  0.948698058 , -70.          ,  -0.0497191584],
       [  0.948698058 , -50.          ,  -0.0497191584],
       [  0.948698058 , -30.          ,  -0.0497191584],
       [  0.948698058 , -10.          ,  -0.0497191584],
       [  0.948698058 ,  10.          ,  -0.0497191584],
       [  0.948698058 ,  30.          ,  -0.0497191584],
       [  0.948698058 ,  50.          ,  -0.0497191584],
       [  0.948698058 ,  70.          ,  -0.0497191584],
       [  0.948698058 ,  90.          ,  -0.0497191584]])

# rings
rings = np.array([[[   0.2496573837, -100.          ,   -0.0130839891],
        [   0.0499314767, -100.          ,   -0.0026167978],
        [   0.0499314767,  -80.          ,   -0.0026167978],
        [   0.2496573837,  -80.          ,   -0.0130839891]],

       [[   0.2496573837,  -80.          ,   -0.0130839891],
        [   0.0499314767,  -80.          ,   -0.0026167978],
        [   0.0499314767,  -60.          ,   -0.0026167978],
        [   0.2496573837,  -60.          ,   -0.0130839891]],

       [[   0.2496573837,  -60.          ,   -0.0130839891],
        [   0.0499314767,  -60.          ,   -0.0026167978],
        [   0.0499314767,  -40.          ,   -0.0026167978],
        [   0.2496573837,  -40.          ,   -0.0130839891]],

       [[   0.2496573837,  -40.          ,   -0.0130839891],
        [   0.0499314767,  -40.          ,   -0.0026167978],
        [   0.0499314767,  -20.          ,   -0.0026167978],
        [   0.2496573837,  -20.          ,   -0.0130839891]],

       [[   0.2496573837,  -20.          ,   -0.0130839891],
        [   0.0499314767,  -20.          ,   -0.0026167978],
        [   0.0499314767,    0.          ,   -0.0026167978],
        [   0.2496573837,    0.          ,   -0.0130839891]],

       [[   0.2496573837,    0.          ,   -0.0130839891],
        [   0.0499314767,    0.          ,   -0.0026167978],
        [   0.0499314767,   20.          ,   -0.0026167978],
        [   0.2496573837,   20.          ,   -0.0130839891]],

       [[   0.2496573837,   20.          ,   -0.0130839891],
        [   0.0499314767,   20.          ,   -0.0026167978],
        [   0.0499314767,   40.          ,   -0.0026167978],
        [   0.2496573837,   40.          ,   -0.0130839891]],

       [[   0.2496573837,   40.          ,   -0.0130839891],
        [   0.0499314767,   40.          ,   -0.0026167978],
        [   0.0499314767,   60.          ,   -0.0026167978],
        [   0.2496573837,   60.          ,   -0.0130839891]],

       [[   0.2496573837,   60.          ,   -0.0130839891],
        [   0.0499314767,   60.          ,   -0.0026167978],
        [   0.0499314767,   80.          ,   -0.0026167978],
        [   0.2496573837,   80.          ,   -0.0130839891]],

       [[   0.2496573837,   80.          ,   -0.0130839891],
        [   0.0499314767,   80.          ,   -0.0026167978],
        [   0.0499314767,  100.          ,   -0.0026167978],
        [   0.2496573837,  100.          ,   -0.0130839891]],

       [[   0.4493832906, -100.          ,   -0.0235511803],
        [   0.2496573837, -100.          ,   -0.0130839891],
        [   0.2496573837,  -80.          ,   -0.0130839891],
        [   0.4493832906,  -80.          ,   -0.0235511803]],

       [[   0.4493832906,  -80.          ,   -0.0235511803],
        [   0.2496573837,  -80.          ,   -0.0130839891],
        [   0.2496573837,  -60.          ,   -0.0130839891],
        [   0.4493832906,  -60.          ,   -0.0235511803]],

       [[   0.4493832906,  -60.          ,   -0.0235511803],
        [   0.2496573837,  -60.          ,   -0.0130839891],
        [   0.2496573837,  -40.          ,   -0.0130839891],
        [   0.4493832906,  -40.          ,   -0.0235511803]],

       [[   0.4493832906,  -40.          ,   -0.0235511803],
        [   0.2496573837,  -40.          ,   -0.0130839891],
        [   0.2496573837,  -20.          ,   -0.0130839891],
        [   0.4493832906,  -20.          ,   -0.0235511803]],

       [[   0.4493832906,  -20.          ,   -0.0235511803],
        [   0.2496573837,  -20.          ,   -0.0130839891],
        [   0.2496573837,    0.          ,   -0.0130839891],
        [   0.4493832906,    0.          ,   -0.0235511803]],

       [[   0.4493832906,    0.          ,   -0.0235511803],
        [   0.2496573837,    0.          ,   -0.0130839891],
        [   0.2496573837,   20.          ,   -0.0130839891],
        [   0.4493832906,   20.          ,   -0.0235511803]],

       [[   0.4493832906,   20.          ,   -0.0235511803],
        [   0.2496573837,   20.          ,   -0.0130839891],
        [   0.2496573837,   40.          ,   -0.0130839891],
        [   0.4493832906,   40.          ,   -0.0235511803]],

       [[   0.4493832906,   40.          ,   -0.0235511803],
        [   0.2496573837,   40.          ,   -0.0130839891],
        [   0.2496573837,   60.          ,   -0.0130839891],
        [   0.4493832906,   60.          ,   -0.0235511803]],

       [[   0.4493832906,   60.          ,   -0.0235511803],
        [   0.2496573837,   60.          ,   -0.0130839891],
        [   0.2496573837,   80.          ,   -0.0130839891],
        [   0.4493832906,   80.          ,   -0.0235511803]],

       [[   0.4493832906,   80.          ,   -0.0235511803],
        [   0.2496573837,   80.          ,   -0.0130839891],
        [   0.2496573837,  100.          ,   -0.0130839891],
        [   0.4493832906,  100.          ,   -0.0235511803]],

       [[   0.6491091976, -100.          ,   -0.0340183716],
        [   0.4493832906, -100.          ,   -0.0235511803],
        [   0.4493832906,  -80.          ,   -0.0235511803],
        [   0.6491091976,  -80.          ,   -0.0340183716]],

       [[   0.6491091976,  -80.          ,   -0.0340183716],
        [   0.4493832906,  -80.          ,   -0.0235511803],
        [   0.4493832906,  -60.          ,   -0.0235511803],
        [   0.6491091976,  -60.          ,   -0.0340183716]],

       [[   0.6491091976,  -60.          ,   -0.0340183716],
        [   0.4493832906,  -60.          ,   -0.0235511803],
        [   0.4493832906,  -40.          ,   -0.0235511803],
        [   0.6491091976,  -40.          ,   -0.0340183716]],

       [[   0.6491091976,  -40.          ,   -0.0340183716],
        [   0.4493832906,  -40.          ,   -0.0235511803],
        [   0.4493832906,  -20.          ,   -0.0235511803],
        [   0.6491091976,  -20.          ,   -0.0340183716]],

       [[   0.6491091976,  -20.          ,   -0.0340183716],
        [   0.4493832906,  -20.          ,   -0.0235511803],
        [   0.4493832906,    0.          ,   -0.0235511803],
        [   0.6491091976,    0.          ,   -0.0340183716]],

       [[   0.6491091976,    0.          ,   -0.0340183716],
        [   0.4493832906,    0.          ,   -0.0235511803],
        [   0.4493832906,   20.          ,   -0.0235511803],
        [   0.6491091976,   20.          ,   -0.0340183716]],

       [[   0.6491091976,   20.          ,   -0.0340183716],
        [   0.4493832906,   20.          ,   -0.0235511803],
        [   0.4493832906,   40.          ,   -0.0235511803],
        [   0.6491091976,   40.          ,   -0.0340183716]],

       [[   0.6491091976,   40.          ,   -0.0340183716],
        [   0.4493832906,   40.          ,   -0.0235511803],
        [   0.4493832906,   60.          ,   -0.0235511803],
        [   0.6491091976,   60.          ,   -0.0340183716]],

       [[   0.6491091976,   60.          ,   -0.0340183716],
        [   0.4493832906,   60.          ,   -0.0235511803],
        [   0.4493832906,   80.          ,   -0.0235511803],
        [   0.6491091976,   80.          ,   -0.0340183716]],

       [[   0.6491091976,   80.          ,   -0.0340183716],
        [   0.4493832906,   80.          ,   -0.0235511803],
        [   0.4493832906,  100.          ,   -0.0235511803],
        [   0.6491091976,  100.          ,   -0.0340183716]],

       [[   0.8488351045, -100.          ,   -0.0444855628],
        [   0.6491091976, -100.          ,   -0.0340183716],
        [   0.6491091976,  -80.          ,   -0.0340183716],
        [   0.8488351045,  -80.          ,   -0.0444855628]],

       [[   0.8488351045,  -80.          ,   -0.0444855628],
        [   0.6491091976,  -80.          ,   -0.0340183716],
        [   0.6491091976,  -60.          ,   -0.0340183716],
        [   0.8488351045,  -60.          ,   -0.0444855628]],

       [[   0.8488351045,  -60.          ,   -0.0444855628],
        [   0.6491091976,  -60.          ,   -0.0340183716],
        [   0.6491091976,  -40.          ,   -0.0340183716],
        [   0.8488351045,  -40.          ,   -0.0444855628]],

       [[   0.8488351045,  -40.          ,   -0.0444855628],
        [   0.6491091976,  -40.          ,   -0.0340183716],
        [   0.6491091976,  -20.          ,   -0.0340183716],
        [   0.8488351045,  -20.          ,   -0.0444855628]],

       [[   0.8488351045,  -20.          ,   -0.0444855628],
        [   0.6491091976,  -20.          ,   -0.0340183716],
        [   0.6491091976,    0.          ,   -0.0340183716],
        [   0.8488351045,    0.          ,   -0.0444855628]],

       [[   0.8488351045,    0.          ,   -0.0444855628],
        [   0.6491091976,    0.          ,   -0.0340183716],
        [   0.6491091976,   20.          ,   -0.0340183716],
        [   0.8488351045,   20.          ,   -0.0444855628]],

       [[   0.8488351045,   20.          ,   -0.0444855628],
        [   0.6491091976,   20.          ,   -0.0340183716],
        [   0.6491091976,   40.          ,   -0.0340183716],
        [   0.8488351045,   40.          ,   -0.0444855628]],

       [[   0.8488351045,   40.          ,   -0.0444855628],
        [   0.6491091976,   40.          ,   -0.0340183716],
        [   0.6491091976,   60.          ,   -0.0340183716],
        [   0.8488351045,   60.          ,   -0.0444855628]],

       [[   0.8488351045,   60.          ,   -0.0444855628],
        [   0.6491091976,   60.          ,   -0.0340183716],
        [   0.6491091976,   80.          ,   -0.0340183716],
        [   0.8488351045,   80.          ,   -0.0444855628]],

       [[   0.8488351045,   80.          ,   -0.0444855628],
        [   0.6491091976,   80.          ,   -0.0340183716],
        [   0.6491091976,  100.          ,   -0.0340183716],
        [   0.8488351045,  100.          ,   -0.0444855628]],

       [[   1.0485610115, -100.          ,   -0.0549527541],
        [   0.8488351045, -100.          ,   -0.0444855628],
        [   0.8488351045,  -80.          ,   -0.0444855628],
        [   1.0485610115,  -80.          ,   -0.0549527541]],

       [[   1.0485610115,  -80.          ,   -0.0549527541],
        [   0.8488351045,  -80.          ,   -0.0444855628],
        [   0.8488351045,  -60.          ,   -0.0444855628],
        [   1.0485610115,  -60.          ,   -0.0549527541]],

       [[   1.0485610115,  -60.          ,   -0.0549527541],
        [   0.8488351045,  -60.          ,   -0.0444855628],
        [   0.8488351045,  -40.          ,   -0.0444855628],
        [   1.0485610115,  -40.          ,   -0.0549527541]],

       [[   1.0485610115,  -40.          ,   -0.0549527541],
        [   0.8488351045,  -40.          ,   -0.0444855628],
        [   0.8488351045,  -20.          ,   -0.0444855628],
        [   1.0485610115,  -20.          ,   -0.0549527541]],

       [[   1.0485610115,  -20.          ,   -0.0549527541],
        [   0.8488351045,  -20.          ,   -0.0444855628],
        [   0.8488351045,    0.          ,   -0.0444855628],
        [   1.0485610115,    0.          ,   -0.0549527541]],

       [[   1.0485610115,    0.          ,   -0.0549527541],
        [   0.8488351045,    0.          ,   -0.0444855628],
        [   0.8488351045,   20.          ,   -0.0444855628],
        [   1.0485610115,   20.          ,   -0.0549527541]],

       [[   1.0485610115,   20.          ,   -0.0549527541],
        [   0.8488351045,   20.          ,   -0.0444855628],
        [   0.8488351045,   40.          ,   -0.0444855628],
        [   1.0485610115,   40.          ,   -0.0549527541]],

       [[   1.0485610115,   40.          ,   -0.0549527541],
        [   0.8488351045,   40.          ,   -0.0444855628],
        [   0.8488351045,   60.          ,   -0.0444855628],
        [   1.0485610115,   60.          ,   -0.0549527541]],

       [[   1.0485610115,   60.          ,   -0.0549527541],
        [   0.8488351045,   60.          ,   -0.0444855628],
        [   0.8488351045,   80.          ,   -0.0444855628],
        [   1.0485610115,   80.          ,   -0.0549527541]],

       [[   1.0485610115,   80.          ,   -0.0549527541],
        [   0.8488351045,   80.          ,   -0.0444855628],
        [   0.8488351045,  100.          ,   -0.0444855628],
        [   1.0485610115,  100.          ,   -0.0549527541]]])


##### vlm

#  "It's strongly recommended that you don't specify type 
# signatures unless you really need to, Numba is pretty good at just working them out." ~ numba collaborator
@numba.jit(nopython=True)
def is_in_vortex_core(vector_list):
    #todo: polepszyc to
    for vec in vector_list:
        if norm(vec) < 1e-9:
            return True
    return False
        
#@numba.jit(nopython=True) -> slower than version below
@numba.jit(numba.float64[::1](numba.float64[::1], numba.float64[::1], numba.float64[::1], numba.optional(numba.float64)), nopython=True, debug = True) 
def vortex_line(p: np.array, p1: np.array, p2: np.array, gamma: float = 1.0) -> np.array:
#def vortex_line(p, p1,  p2,  gamma = 1.0):
    # strona 254

    r0 = np.asarray(p2 - p1)
    #print(numba.typeof(r0))
    r1 = np.asarray(p - p1)
    r2 = np.asarray(p - p2)
    
    r1_cross_r2 = np.cross(r1, r2)
    
    q_ind = np.array([0.0, 0.0, 0.0], dtype=np.float64)
    # in nonpython mode must be list reflection to convert list to non python type
    # nested python oject can be badly converted -> recommend to use numba.typed.List
    b = is_in_vortex_core(List([r1, r2, r1_cross_r2]))
    if b:
        return np.asarray([0.0, 0.0, 0.0], dtype=np.float64)
    else:
        q_ind = r1_cross_r2 / np.square(np.linalg.norm(r1_cross_r2))
        q_ind *= np.dot(r0, (r1 / np.linalg.norm(r1) - r2 / np.linalg.norm(r2)))
        q_ind *= gamma / (4 * np.pi)

    return q_ind

##### end vilm functions


# DO NOT REPORT THIS... COMPILATION TIME IS INCLUDED IN THE EXECUTION TIME!
start = time.time()
vortex_line(collocation_points[0], rings[0,0], rings[0,1], gamma= 1.0)
end = time.time()
print("Elapsed (with compilation) = %s" % (end - start))

# NOW THE FUNCTION IS COMPILED, RE-TIME IT EXECUTING FROM CACHE
start = time.time()
vortex_line(collocation_points[0], rings[0,0], rings[0,1], gamma= 1.0)
end = time.time()
print("Elapsed (after compilation) = %s" % (end - start))


# > python sailingVLM/NewApproach/jit_testy.py 
# Elapsed (with compilation) = 0.00041294097900390625
# Elapsed (after compilation) = 1.6450881958007812e-05