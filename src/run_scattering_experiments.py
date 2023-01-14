import rebound
from astropy import units as u
import numpy as np
import argparse
import multiprocessing
import time
import numba
import h5py
from ctypes import *
import os
#from tqdm import tqdm
#from joblib import Parallel, delayed
#from parallelbar import progress_imap
#multiprocessing.set_start_method("spawn")
import threading
from functools import wraps
import signal



script_dir = os.path.abspath(os.path.dirname(__file__))
lib_path = os.path.join(script_dir, "keplerPy.so")

libkepler = cdll.LoadLibrary(lib_path)
libkepler.integrate.argtypes = [POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double),
                                POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double),
                                c_double, POINTER(c_int)]

script_dir = os.path.abspath(os.path.dirname(__file__))
lib_path = os.path.join(script_dir, "twobodyPy.so")

libtwobody = cdll.LoadLibrary(lib_path)
libtwobody.integrate2body.argtypes = [POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double),
                                c_double, POINTER(c_double), POINTER(c_double)]



def stop_function():
    os.kill(os.getpid(), signal.SIGINT)


def stopit_after_timeout(s, raise_exception=True):
    def actual_decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            timer = threading.Timer(s, stop_function)
            try:
                timer.start()
                result = func(*args, **kwargs)
            except KeyboardInterrupt:
                msg = f'function {func.__name__} took longer than {s} s.'
                if raise_exception:
                    raise TimeoutError(msg)
                result = msg
            finally:
                timer.cancel()
            return result

        return wrapper

    return actual_decorator


def concatenate_res2(res_data,nsample):
    xe = np.zeros(nsample)
    ye = np.zeros(nsample)
    ze = np.zeros(nsample)
    vxe = np.zeros(nsample)
    vye = np.zeros(nsample)
    vze = np.zeros(nsample)

    xm = np.zeros(nsample)
    ym = np.zeros(nsample)
    zm = np.zeros(nsample)
    vxm = np.zeros(nsample)
    vym = np.zeros(nsample)
    vzm = np.zeros(nsample)
    for i in range(0,len(res_data)):
        xe[i] = res_data[i][0]
        ye[i] = res_data[i][1]
        ze[i] = res_data[i][2]
        vxe[i] = res_data[i][3]
        vye[i] = res_data[i][4]
        vze[i] = res_data[i][5]
        xm[i] = res_data[i][6]
        ym[i] = res_data[i][7]
        zm[i] = res_data[i][8]
        vxm[i] = res_data[i][9]
        vym[i] = res_data[i][10]
        vzm[i] = res_data[i][11]


    return xe,ye,ze,vxe,vye,vze,xm,ym,zm,vxm,vym,vzm



def concatenate_res(res_data,nsample):
    xe = np.zeros(nsample)
    ye = np.zeros(nsample)
    ze = np.zeros(nsample)
    vxe = np.zeros(nsample)
    vye = np.zeros(nsample)
    vze = np.zeros(nsample)

    xm = np.zeros(nsample)
    ym = np.zeros(nsample)
    zm = np.zeros(nsample)
    vxm = np.zeros(nsample)
    vym = np.zeros(nsample)
    vzm = np.zeros(nsample)

    for i in range(0,len(res_data)):
        xe[i*1000:(i+1)*1000] = res_data[i][0]
        ye[i*1000:(i+1)*1000] = res_data[i][1]
        ze[i*1000:(i+1)*1000] = res_data[i][2]
        vxe[i*1000:(i+1)*1000] = res_data[i][3]
        vye[i*1000:(i+1)*1000] = res_data[i][4]
        vze[i*1000:(i+1)*1000] = res_data[i][5]
        xm[i*1000:(i+1)*1000] = res_data[i][6]
        ym[i*1000:(i+1)*1000] = res_data[i][7]
        zm[i*1000:(i+1)*1000] = res_data[i][8]
        vxm[i*1000:(i+1)*1000] = res_data[i][9]
        vym[i*1000:(i+1)*1000] = res_data[i][10]
        vzm[i*1000:(i+1)*1000] = res_data[i][11]


    return xe,ye,ze,vxe,vye,vze,xm,ym,zm,vxm,vym,vzm

def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-b','--batchsize',dest='batchsize',action='store',default=100)
    parser.add_argument('-n','--nbatch',dest='nbatch',action='store',default=10000)
    parser.add_argument('-v','--vasym',dest='v_asym',action='store',default=1.0)
    args = parser.parse_args()
    return int(args.batchsize), int(args.nbatch), float(args.v_asym)

def generate_em_posvel(nsample):
    sim = rebound.Simulation()
    sim.units = ('yr', 'AU', 'Msun')

    #sim.add(m=1.0)
    me = 5.9724e24 * u.kg
    me = me.to(u.Msun).value
    mm = 7.34767309e22 * u.kg
    mm = mm.to(u.Msun).value

    am = 384748 * u.km
    am = am.to(u.AU).value
    xe = []
    vxe = []
    xm = []
    vxm = []
    rng = np.random.default_rng()

    for i in range(0,nsample):
        f1 = rng.random() * 2 * np.pi
        f2 = rng.random() * 2 * np.pi
        #f2 = np.random.random_sample()*2*np.pi
        sim.add(m=1.0)
        sim.add(primary=sim.particles[0],m=me,a=1,e=0.0167,f=f1)
        sim.add(primary=sim.particles[1],m=mm,a=am,e=0.0549,f=f2)


        earth = sim.particles[1]
        moon = sim.particles[2]

        xe.append([earth.x,earth.y,earth.z])
        vxe.append([earth.vx,earth.vy,earth.vz])

        xm.append([moon.x,moon.y,moon.z])
        vxm.append([moon.vx,moon.vy,moon.vz])

        del sim.particles

    sim = None

    return xe,xm,vxe,vxm

#@numba.jit(nopython=True)
def generate_vuv_and_r(nsample):
    rng = np.random.default_rng()
    r = 1e9 #in AU
    phi = 2*np.pi * rng.random(nsample)
    theta = np.pi * rng.random(nsample)
    x = r * np.cos(phi) * np.sin(theta)
    y = r * np.sin(phi) * np.sin(theta)
    z = r * np.cos(theta)

    #vuv = np.zeros(nsample)
    #xyz = np.zeros(nsample)
    xyz = np.array((x,y,z))
    xyz = xyz.T
    vuv = -xyz/r
    #for i in range(0,nsample):
    #    vuv[i] = np.array([-x[i]/r,-y[i]/r,-z[i]/r])
        #vuv.append(np.array([-x[i]/r,-y[i]/r,-z[i]/r]))
    #    xyz[i] = np.array([x[i],y[i],z[i]])
        #xyz.append(np.array([x[i],y[i],z[i]]))

    #vuv = np.array(vuv)
   #xyz = np.array(xyz)
    return vuv,xyz

#@numba.njit(parallel=True)
def get_tangent_plane(nsample,dist,vuv,xyz):
    rng = np.random.default_rng()
    vall = np.zeros((nsample,3))
    for i in range(0,nsample):
        x = -vuv[i][0]
        y = -vuv[i][1]
        z = -vuv[i][2]
        d = np.array([x,y,z])
        n = [2*x,2*y,2*z]
        j = np.random.randint(low=0,high=3)
        coords = np.delete([0,1,2],j)
        tx = -5 + 10 * rng.random()
        ty = -5 + 10 * rng.random()
        tz = (-(tx-d[coords[0]])*n[coords[0]] + -(ty-d[coords[1]])*n[coords[1]])/n[j] + d[j]
        t = np.zeros(3)
        t[coords[0]]=tx-d[coords[0]]
        t[coords[1]]=ty-d[coords[1]]
        t[j]=tz-d[j]
        modt = np.sqrt(t[0]**2 + t[1]**2 + t[2]**2)
        vhat = np.array([t[0]/modt, t[1]/modt, t[2]/modt])

        v = vhat * dist[i]
        v[0] += xyz[i][0]
        v[1] += xyz[i][1]
        v[2] += xyz[i][2]

        vall[i][0]=v[0]
        vall[i][1]=v[1]
        vall[i][2]=v[2]
        #vall.append(v)
        #print(i)

    return vall


#@numba.njit(parallel=True)
def velocity_at_position(r,rb,v_inf,G,me,mm,BATCHSIZE):
    vfinal = np.zeros(BATCHSIZE)
    for i in range(0,BATCHSIZE):
        e_inf = 0.5*v_inf[i]**2
        db = np.sqrt((r[i][0]-rb[i][0])**2 + (r[i][1]-rb[i][1])**2 + (r[i][2]-rb[i][2])**2)
        e_bar = -G*(me+mm+1) / db
    #rsun = r + rb
    #dsun = np.sqrt(rsun[0]**2 + rsun[1]**2 + rsun[2]**2)
    #e_sun = -G / dsun
        vfinal[i] = np.sqrt(2 * (e_inf-e_bar))
    return vfinal



def generate_all_components(nsample,v_a,xe,ye,ze,vxe,vye,vze,xm,ym,zm,vxm,vym,vzm):
    sim = rebound.Simulation()
    sim.units = ('yr', 'AU', 'Msun')

    me = 5.9724e24 * u.kg
    me = me.to(u.Msun).value
    mm = 7.34767309e22 * u.kg
    mm = mm.to(u.Msun).value

    v_asym = v_a #km/s
    v_asym_kms = v_asym * (u.km/u.s)#km/s
    print(v_asym_kms)
    v_asym_auyr = v_asym_kms.to(u.AU/u.yr).value
    s = np.zeros(nsample) + v_asym

    phiv = 2*np.pi * np.random.random_sample(nsample)
    thetav = np.pi * np.random.random_sample(nsample)

    vx = s * np.cos(phiv) * np.sin(thetav)
    vy = s * np.sin(phiv) * np.sin(thetav)
    vz = s * np.cos(thetav)

    vx = vx * (u.km / u.s)
    vy = vy * (u.km / u.s)
    vz = vz * (u.km / u.s)

    vx = vx.to(u.AU/u.yr).value
    vy = vy.to(u.AU/u.yr).value
    vz = vz.to(u.AU/u.yr).value

    v_inf = np.sqrt(vx**2 + vy**2 + vz**2)

    vuv, xyz = generate_vuv_and_r(nsample)
    upper_b = np.sqrt(1 + 2 * sim.G*(me+mm+1)/(v_asym_auyr**2))
    print(upper_b)
    b = np.random.random_sample(nsample) * (upper_b) #* (20*Re-Re) + Re

    xyz_tang = get_tangent_plane(nsample,b,vuv,xyz)


    #get velocity of barycenter
    vxb = (me*vxe + mm*vxm)/(me+mm+1)
    vyb = (me*vye + mm*vym)/(me+mm+1)
    vzb = (me*vze + mm*vzm)/(me+mm+1)

    #get position of barycenter
    xb = (me*xe + mm*xm)/(me+mm+1)
    yb = (me*ye + mm*ym)/(me+mm+1)
    zb = (me*ze + mm*zm)/(me+mm+1)

    #use conservation of energy to calculate the velocity
    rb = np.array((xb,yb,zb))
    rb = rb.T
    vfinal = velocity_at_position(xyz_tang,rb,v_inf,sim.G,me,mm)


    vel_parts=[]
    for i in range(0,nsample):
        vel_parts.append(vuv[i]*vfinal[i])

    vel_parts = np.array(vel_parts)
    xyz_tang[:,0] += xb
    xyz_tang[:,1] += yb
    xyz_tang[:,2] += zb

    vel_parts[:,0] += vxb
    vel_parts[:,1] += vyb
    vel_parts[:,2] += vzb

    return upper_b, b, xyz_tang, vel_parts

def run_kepler_mod(xb,xc,vxb,vxc,t_end):

    xci = [0.0,0.0,0.0]
    vxci = [0.0,0.0,0.0]
    xb_c = (c_double * 3)(*xb)
    xc_c = (c_double * 3)(*xc)
    xci_c = (c_double *3)(*xci)
    vxb_c = (c_double * 3)(*vxb)
    vxc_c = (c_double * 3)(*vxc)
    vxci_c = (c_double *3)(*vxci)

    libtwobody.integrate2body(xb_c, xc_c, vxb_c, vxc_c, c_double(t_end), xci_c, vxci_c)

    comet_params = np.array([xci_c[0],xci_c[1],xci_c[2],vxci_c[0],vxci_c[1],vxci_c[2]])
    return comet_params

def run_forward_kepler(xe,xm,xc,vxe,vxm,vxc,t_end):
    xs = [0.0,0.0,0.0]
    vxs = [0.0,0.0,0.0]

    xs_c = (c_double * 3)(*xs)

    xe_c = (c_double * 3)(*xe)
    xm_c = (c_double * 3)(*xm)
    xc_c = (c_double * 3)(*xc)

    vxs_c = (c_double * 3)(*vxs)
    vxe_c = (c_double * 3)(*vxe)
    vxm_c = (c_double * 3)(*vxm)
    vxc_c = (c_double * 3)(*vxc)

    collb_flag = c_int(0)

    #sim.move_to_com()
    #t_end = 2*orbit.T
    libkepler.integrate(xs_c,xe_c,xm_c,xc_c,vxs_c,vxe_c,vxm_c,vxc_c,c_double(t_end),byref(collb_flag))


    #pack the data properly and return it


    #if(collb_flag.value == 4):
        #print("aye!")
    #    isunbound = False

    #DEBUG START
    #if(collb_flag.value == 4):
    #    sim = rebound.Simulation()
    #    sim = rebound.Simulation()
    #    sim.units = ('yr', 'AU', 'Msun')

    #    me = 5.9724e24 * u.kg
    #    me = me.to(u.Msun).value
    #    mm = 7.34767309e22 * u.kg
    #    mm = mm.to(u.Msun).value

    #    sim.add(m=1.0)
    #    sim.add(m=me,x=xe[0],y=xe[1],z=xe[2],vx=vxe[0],vy=vxe[1],vz=vxe[2])
    #    sim.add(m=mm,x=xm[0],y=xm[1],z=xm[2],vx=vxm[0],vy=vxm[1],vz=vxm[2])
    #    sim.add(m=0.0,x=xc[0],y=xc[1],z=xc[2],vx=vxc[0],vy=vxc[1],vz=vxc[2])
    #    sim.move_to_com()
    #    sim.integrate(t_end)
    #    comet = sim.particles[3]
    #    print('-------------------')
    #    print((comet.x-xc_c[0])/comet.x,(comet.y-xc_c[1])/comet.y,(comet.z-xc_c[2])/comet.z)
    #    print((comet.vx-vxc_c[0])/comet.vx,(comet.vy-vxc_c[1])/comet.vy,(comet.vz-vxc_c[2])/comet.vz)
    #    orbit = sim.particles[3].calculate_orbit(primary=sim.particles[0])
    #    print(orbit.a)
    #    print('-------------------')
    #DEBUG END

    xsi = np.array([xs_c[0],xs_c[1],xs_c[2]])
    xei = np.array([xe_c[0],xe_c[1],xe_c[2]])
    xmi = np.array([xm_c[0],xm_c[1],xm_c[2]])
    xpi = np.array([xc_c[0],xc_c[1],xc_c[2]])

    vsi = np.array([vxs_c[0],vxs_c[1],vxs_c[2]])
    vei = np.array([vxe_c[0],vxe_c[1],vxe_c[2]])
    vmi = np.array([vxm_c[0],vxm_c[1],vxm_c[2]])
    vpi = np.array([vxc_c[0],vxc_c[1],vxc_c[2]])

    return (xsi,xei,xmi,xpi,vsi,vei,vmi,vpi,collb_flag.value)
    #return (isunbound,xsi,xei,xmi,xpi,vsi,vei,vmi,vpi)


@stopit_after_timeout(10, raise_exception=True)
def run_scattering_experiment(args):
    #set BATCHSIZE to be 100
    #generate BATCHSIZE random earth-moon posvels first
    v_a = args[0]
    BATCHSIZE=int(args[1])

    xe,xm,vxe,vxm = generate_em_posvel(BATCHSIZE)

    #generate the positions of the comets next
    sim = rebound.Simulation()
    sim.units = ('yr', 'AU', 'Msun')

    me = 5.9724e24 * u.kg
    me = me.to(u.Msun).value
    mm = 7.34767309e22 * u.kg
    mm = mm.to(u.Msun).value

    v_asym = v_a #km/s
    v_asym_kms = v_asym * (u.km/u.s)#km/s
    #print(v_asym_kms)
    v_asym_auyr = v_asym_kms.to(u.AU/u.yr).value
    s = np.zeros(BATCHSIZE) + v_asym

    #use rng for better randomness
    rng = np.random.default_rng()


    phiv = 2*np.pi * rng.random(BATCHSIZE)
    thetav = np.pi * rng.random(BATCHSIZE)

    vx = s * np.cos(phiv) * np.sin(thetav)
    vy = s * np.sin(phiv) * np.sin(thetav)
    vz = s * np.cos(thetav)

    vx = vx * (u.km / u.s)
    vy = vy * (u.km / u.s)
    vz = vz * (u.km / u.s)

    vx = vx.to(u.AU/u.yr).value
    vy = vy.to(u.AU/u.yr).value
    vz = vz.to(u.AU/u.yr).value

    v_inf = np.sqrt(vx**2 + vy**2 + vz**2)

    vuv, xyz = generate_vuv_and_r(BATCHSIZE)
    upper_b = np.sqrt(1 + 2 * sim.G*(me+mm+1)/(v_asym_auyr**2))
    #print(upper_b)
    b = rng.random(BATCHSIZE) * (upper_b) #* (20*Re-Re) + Re

    #DEBUG STUFF START
    #a = sim.G * (1+me+mm) / v_inf**2
    #e2 = 1 + b**2 / a**2
    #ec = np.sqrt(e2)
    #rp = a * (ec - 1)
    #print(rp)
    #DEBUG STUFF ENDS


    xyz_tang = get_tangent_plane(BATCHSIZE,b,vuv,xyz)

    #get velocity of barycenter
    vxb = []
    vyb = []
    vzb = []
    xb = []
    yb = []
    zb = []

    for i in range(0,BATCHSIZE):
        vxb.append((me*vxe[i][0] + mm*vxm[i][0])/(me+mm+1))
        vyb.append((me*vxe[i][1] + mm*vxm[i][1])/(me+mm+1))
        vzb.append((me*vxe[i][2] + mm*vxm[i][2])/(me+mm+1))

    #get position of barycenter
        xb.append((me*xe[i][0] + mm*xm[i][0])/(me+mm+1))
        yb.append((me*xe[i][1] + mm*xm[i][1])/(me+mm+1))
        zb.append((me*xe[i][2] + mm*xm[i][2])/(me+mm+1))

    vxb = np.array(vxb)
    vyb = np.array(vyb)
    vzb = np.array(vzb)
    xb = np.array(xb)
    yb = np.array(yb)
    zb = np.array(zb)

    #use conservation of energy to calculate the velocity
    rb = np.array((xb,yb,zb))
    rb = rb.T
    vfinal = velocity_at_position(xyz_tang,rb,v_inf,sim.G,me,mm,BATCHSIZE)

    vel_parts=[]
    for i in range(0,BATCHSIZE):
        vel_parts.append(vuv[i]*vfinal[i])

    vel_parts = np.array(vel_parts)
    xyz_tang[:,0] += xb
    xyz_tang[:,1] += yb
    xyz_tang[:,2] += zb

    vel_parts[:,0] += vxb
    vel_parts[:,1] += vyb
    vel_parts[:,2] += vzb

    #now run the simulation in two steps
    #1. advance kepler solution to 20 AU
    #2. scattering experiment until bound object found

    #find time of periastron passage

    #DEBUG STUFF START

    #rp_debug = []
    #rp_debug2 = []
#     tp = []

    #for i in range(0,BATCHSIZE):
        #sim = rebound.Simulation()
        #sim.units = ('yr', 'AU', 'Msun')


        #sim.add(m=(1+me+mm),x=xb[i],y=yb[i],z=zb[i],vx=vxb[i],vy=vyb[i],vz=vzb[i])
        #sim.add(m=0.0,x=xyz_tang[i][0],y=xyz_tang[i][1],z=xyz_tang[i][2],vx=vel_parts[i][0],vy=vel_parts[i][1],vz=vel_parts[i][2])
        #sim.move_to_com()
# #         comet = sim.particles[1]
        #orbit = sim.particles[1].calculate_orbit(primary=sim.particles[0])
# #         #print(orbit.T)
# #         tp.append(orbit.T)
#         rp_debug.append(orbit.a*(1-orbit.e))
        #sim.integrate(orbit.T-20)
# #         orbit = sim.particles[1].calculate_orbit(primary=sim.particles[0])
# #         rp_debug2.append(orbit.a*(1-orbit.e))

        #del sim.particles
        #sim = None
# #     #print((np.array(rp_debug)-np.array(rp_debug2))/np.array(rp_debug))
# #     print(tp)
    #DEBUG STUFF END
    sim.add(m=(1+me+mm),x=xb[0],y=yb[0],z=zb[0],vx=vxb[0],vy=vyb[0],vz=vzb[0])
    sim.add(m=0.0,x=xyz_tang[0][0],y=xyz_tang[0][1],z=xyz_tang[0][2],vx=vel_parts[0][0],vy=vel_parts[0][1],vz=vel_parts[0][2])
    sim.move_to_com()
    orbit = sim.particles[1].calculate_orbit(primary=sim.particles[0])
    t_end = orbit.T+10 #buffer
    del sim.particles
    res_batch = []
    for i in range(0,BATCHSIZE):
        xxb = [xb[i],yb[i],zb[i]]
        vvxb = [vxb[i],vyb[i],vzb[i]]
        xxc = [xyz_tang[i][0],xyz_tang[i][1],xyz_tang[i][2]]
        vvxc = [vel_parts[i][0],vel_parts[i][1],vel_parts[i][2]]
        cp = run_kepler_mod(xxb,xxc,vvxb,vvxc,t_end)

        sim.add(m=1,x=0.0,y=0.0,z=0.0,vx=0.0,vy=0.0,vz=0.0)
        sim.add(m=0.0,x=cp[0],y=cp[1],z=cp[2],vx=cp[3],vy=cp[4],vz=cp[5])
        sim.move_to_com()
        orbit = sim.particles[1].calculate_orbit(primary=sim.particles[0])
        qi = (orbit.a*(1-orbit.e))
        #rp_debug.append(qi)
        t_end_s = 2*orbit.T
        del sim.particles
        xc = [cp[0],cp[1],cp[2]]
        vxc = [cp[3],cp[4],cp[5]]

        res_s = run_forward_kepler(xe[i],xm[i],xc,vxe[i],vxm[i],vxc,t_end_s)
        if(res_s[-1]>0):
            #xs,xe,xm,xc,vxs,vxe,vxm,vxc,initial_b,initial_rp,collb_flag
            res_batch.append([res_s[0],res_s[1],res_s[2],res_s[3],res_s[4],res_s[5],res_s[6],res_s[7],b[i],qi,res_s[-1]])
    #print("finished batch!")
    #print(((np.array(rp_debug)-rp)/rp))
    return res_batch
    #print((np.array(rp_debug)-np.array(rp_debug2))/np.array(rp_debug))
    #print(rp_debug)
if __name__ == "__main__":
    #multiprocessing.set_start_method('forkserver')
    #BATCHSIZE = 1000
    #v_a = 1.0
    #args = (v_a,BATCHSIZE)
    #NBATCHES = 2000
    BATCHSIZE,NBATCHES,v_a = get_arguments()
    args = (v_a, BATCHSIZE)
    args_rep = ((args, ) * NBATCHES)
    #run_scattering_experiment(args)
#     nsample, v_a = get_arguments()
    start_time = time.time()
    res_data = []
    count = 0
    count_capbound = 0

    print("RAPTOR v1.0. Made by Diptajyoti Mukherjee. Jan 13, 2023",flush=True)
    print("System: Sun-Earth-Moon",flush=True)
    print("BATCHSIZE: ",BATCHSIZE,flush=True)
    print("NBATCHES: ", NBATCHES,flush=True)
    print("Total experiments: ", BATCHSIZE*NBATCHES,flush=True)
    print("V_asym (km/s): ", v_a,flush=True)

    sim = rebound.Simulation()
    sim.units = ('yr', 'AU', 'Msun')

    me = 5.9724e24 * u.kg
    me = me.to(u.Msun).value
    mm = 7.34767309e22 * u.kg
    mm = mm.to(u.Msun).value

    v_asym = v_a #km/s
    v_asym_kms = v_asym * (u.km/u.s)#km/s
    #print(v_asym_kms)
    v_asym_auyr = v_asym_kms.to(u.AU/u.yr).value
    upper_b = np.sqrt(1 + 2 * sim.G*(me+mm+1)/(v_asym_auyr**2))
    print("Max b (AU): {:21.16f}".format(upper_b),flush=True)
    #print("generated earth, moon pos")
    sim = None

    UP = "\x1B[3A"
    CLR = "\x1B[0K"
    print("\n\n")

    #for saving the data
    xs = []
    ys = []
    zs = []
    vxs = []
    vys = []
    vzs = []

    xe = []
    ye = []
    ze = []
    vxe = []
    vye = []
    vze = []

    xm = []
    ym = []
    zm = []
    vxm = []
    vym = []
    vzm = []

    xc = []
    yc = []
    zc = []
    vxc = []
    vyc = []
    vzc = []

    q = []
    b = []
    collb_flag = []
    rrr=[]
    rrre=[]
    with multiprocessing.Pool() as pool:
        iterator = pool.imap_unordered(run_scattering_experiment, args_rep)
        while True:
            try:
                count+=1
                results = iterator.next()
                count_capbound += len(results)
                processed_percent = count/NBATCHES * 100
                print(f"{UP}processed batch: {processed_percent}%{CLR}\ncaptured/bound: {count_capbound}{CLR}\n",flush=True)

                if(len(results)>0):
                    for result in results:
                        xs.append(result[0][0])
                        ys.append(result[0][1])
                        zs.append(result[0][2])

                        xe.append(result[1][0])
                        ye.append(result[1][1])
                        ze.append(result[1][2])

                        xm.append(result[2][0])
                        ym.append(result[2][1])
                        zm.append(result[2][2])

                        xc.append(result[3][0])
                        yc.append(result[3][1])
                        zc.append(result[3][2])

                        vxs.append(result[4][0])
                        vys.append(result[4][1])
                        vzs.append(result[4][2])

                        vxe.append(result[5][0])
                        vye.append(result[5][1])
                        vze.append(result[5][2])

                        vxm.append(result[6][0])
                        vym.append(result[6][1])
                        vzm.append(result[6][2])

                        vxc.append(result[7][0])
                        vyc.append(result[7][1])
                        vzc.append(result[7][2])

                        b.append(result[8])
                        q.append(result[9])
                        collb_flag.append(result[10])

                #print(len(result))
                #rrr.append(result)

            except StopIteration:
                break
            except Exception as e:
                # do something
                rrre.append(e)
        pool.terminate()
        #for results in pool.imap_unordered(run_scattering_experiment, args_rep):
            #if(len(results)>0):
                #count_capbound += len(results)
                #for result in results:
#                     xs.append(result[0][0])
#                     ys.append(result[0][1])
#                     zs.append(result[0][2])

#                     xe.append(result[1][0])
#                     ye.append(result[1][1])
#                     ze.append(result[1][2])

#                     xm.append(result[2][0])
#                     ym.append(result[2][1])
#                     zm.append(result[2][2])

#                     xc.append(result[3][0])
#                     yc.append(result[3][1])
#                     zc.append(result[3][2])

#                     vxs.append(result[4][0])
#                     vys.append(result[4][1])
#                     vzs.append(result[4][2])

#                     vxe.append(result[5][0])
#                     vye.append(result[5][1])
#                     vze.append(result[5][2])

#                     vxm.append(result[6][0])
#                     vym.append(result[6][1])
#                     vzm.append(result[6][2])

#                     vxc.append(result[7][0])
#                     vyc.append(result[7][1])
#                     vzc.append(result[7][2])

#                     b.append(result[8])
#                     q.append(result[9])
#                     collb_flag.append(result[10])
            #processed_percent = count/NBATCHES * 100
            #print(f"{UP}processed batch: {processed_percent}%{CLR}\ncaptured/bound: {count_capbound}{CLR}\n",flush=True)
            #count+=1
            #for result in pool.imap(generate_em_posvel, range(0,nsample)):
            #res_data.append(result)
        #pool.terminate()
#     xe,ye,ze,vxe,vye,vze,xm,ym,zm,vxm,vym,vzm = concatenate_res(res_data,nsample)
#     #print(xe)
    print(count_capbound)
    print("\n\n")
    print(rrre)
    print("total time taken: ",time.time()-start_time)
#     upper_b, b, xyz,vel = generate_all_components(nsample,v_a,xe,ye,ze,vxe,vye,vze,xm,ym,zm,vxm,vym,vzm)
#     print("time taken to generate: ",time.time()-start_time)
    print("writing to file...")

    hf = h5py.File('captured_bound_data.h5', 'w')
    hf.create_dataset('xs',data=xs)
    hf.create_dataset('ys',data=ys)
    hf.create_dataset('zs',data=zs)
    hf.create_dataset('vxs',data=vxs)
    hf.create_dataset('vys',data=vys)
    hf.create_dataset('vzs',data=vzs)


    hf.create_dataset('xe',data=xe)
    hf.create_dataset('ye',data=ye)
    hf.create_dataset('ze',data=ze)
    hf.create_dataset('xm',data=xm)
    hf.create_dataset('ym',data=ym)
    hf.create_dataset('zm',data=zm)

    hf.create_dataset('vxe',data=vxe)
    hf.create_dataset('vye',data=vye)
    hf.create_dataset('vze',data=vze)
    hf.create_dataset('vxm',data=vxm)
    hf.create_dataset('vym',data=vym)
    hf.create_dataset('vzm',data=vzm)

    hf.create_dataset('xc',data=xc)
    hf.create_dataset('yc',data=yc)
    hf.create_dataset('zc',data=zc)
    hf.create_dataset('vxc',data=vxc)
    hf.create_dataset('vyc',data=vyc)
    hf.create_dataset('vzc',data=vzc)

    hf.create_dataset('b',data=b)
    hf.create_dataset('q',data=q)
    hf.create_dataset('collb_flag',data=collb_flag)

    info_dset = hf.create_dataset("info",dtype="f8")
    info_dset.attrs.create("v_asym",data=v_a)
    info_dset.attrs.create("upper_b",data=upper_b)
    info_dset.attrs.create("total_experiments",data=float(count*BATCHSIZE))
    hf.close()

    print("writing file completed")
