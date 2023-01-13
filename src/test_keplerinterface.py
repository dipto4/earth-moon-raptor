import numpy as np
from ctypes import *
import os

script_dir = os.path.abspath(os.path.dirname(__file__))
lib_path = os.path.join(script_dir, "keplerPy.so")

libkepler = cdll.LoadLibrary(lib_path)

libkepler.integrate.argtypes = [POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double),
                                POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double),
                                c_double, POINTER(c_int)]


xs = [2.2781400516148997e-06, -2.0695749954758408e-06, 0.0]
xe = [-7.4921880741332536e-01,6.8064733945096878e-01,0.0]
xm = [-7.5182659042728350e-01,6.8137461212323758e-01,0.0000000000000000e+00]
xc = [1.2946637977571502e+01,4.2149185223069804e-02,1.4926795913517259e+01]

vxs = [1.2849921505053460e-05,1.3830297049357371e-05,0.0]
vxe = [-4.2254649390474457e+00,-4.5462155518353926e+00,0.0]
vxm = [-4.2835868243084372e+00,-4.7427451055342340e+00,0.0000000000000000e+00]
vxc = [-1.3028982911733853e+00,3.4814493538553193e-01,-1.4904828663885017e+00]
#xs = [0.0,0.0,0.0]
#xe = [0.9833000000000001, 0.0, 0.0]
#xm = [0.9857306852303347, 0.0, 0.0]
#xc = [0.999, 0.0, 0.0]

#vxs = [0.0,0.0,0.0]
#vxe = [0.0, 6.388894409407078, 0.0]
#vxm = [0.0, 6.61713329602396, 0.0]
#vxc = [0.0, 6.286210532923412, 0.0]

t_end = 20.0

xs_c = (c_double * 3)(*xs)
xe_c = (c_double * 3)(*xe)
xm_c = (c_double * 3)(*xm)
xc_c = (c_double * 3)(*xc)

vxs_c = (c_double * 3)(*vxs)
vxe_c = (c_double * 3)(*vxe)
vxm_c = (c_double * 3)(*vxm)
vxc_c = (c_double * 3)(*vxc)

coll_flag = c_int(0)

libkepler.integrate(xs_c,xe_c,xm_c,xc_c,vxs_c,vxe_c,vxm_c,vxc_c,c_double(t_end),byref(coll_flag))

print(coll_flag.value)

#for i in xc_c:
#    print(i)

#for i in vxc_c:
#    print(i)
