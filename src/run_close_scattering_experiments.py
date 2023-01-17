#modified version of the original RAPTOR code to handle close scattering events with the earth-moon-system
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


#load the pertinent libraries for the custom integration
script_dir = os.path.abspath(os.path.dirname(__file__))
lib_path = os.path.join(script_dir, "keplerPy.so")

libkepler = cdll.LoadLibrary(lib_path)
libkepler.integrate.argtypes = [POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double),
                                POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double),
                                c_double, POINTER(c_int)]

#script_dir = os.path.abspath(os.path.dirname(__file__))
#lib_path = os.path.join(script_dir, "twobodyPy.so")

#libtwobody = cdll.LoadLibrary(lib_path)
#libtwobody.integrate2body.argtypes = [POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double),
#                                c_double, POINTER(c_double), POINTER(c_double)]


#required for some bizarre multiprocessing bug that prevents the code from finishing. Crude solution but it works!
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

def velocity_at_peri(xe,xc,v_inf,G,me):
    #vfinal = np.zeros(BATCHSIZE)
    #for i in range(0,BATCHSIZE):
    e_inf = 0.5*v_inf**2
    xei = np.array(xe)
    #xmi = np.array(xm[i])
    xci = xei + xc #comet distance is w.r.t. earth
    ds = np.sqrt(np.sum(xci**2,axis=0))
    de = np.sqrt(np.sum((xci-xei)**2,axis=0))
    #dm = np.sqrt(np.sum((xc[i]-xmi)**2,axis=0))
    e_s = -G / ds
    e_e = -G * me / de
    #e_m = -G * mm / dm
    #rsun = r + rb
    #dsun = np.sqrt(rsun[0]**2 + rsun[1]**2 + rsun[2]**2)
    #e_sun = -G / dsun
    vfinal = np.sqrt(2 * (e_inf-e_s-e_e))

    return vfinal

def get_energy(G,me,mm,xsi,xei,xmi,xci,vci):
    v2 = np.sum(vci**2, axis=0)
    ds = np.sqrt(np.sum((xsi-xci)**2,axis=0))
    de = np.sqrt(np.sum((xsi-xci)**2,axis=0))
    dm = np.sqrt(np.sum((xsi-xci)**2,axis=0))
    ene = 0.5*v2 - G*(1.0/ds - me/de - mm / dm)
    return ene



def get_velocity(vperi,x,y,z):
    #vall = []
    rng = np.random.default_rng()
    #get (tx,ty) where t denotes tangent
    #for i in range(0,BATCHSIZE):
    d = np.array([x,y,z])
    n = [2*x,2*y,2*z]

    j = np.random.randint(low=0,high=3)
    coords = np.delete([0,1,2],j)
    tx = -0.01 + 0.02 * rng.random()
    ty = -0.01 + 0.02 * rng.random()
    #solving equation for plane to get tz
    tz = (-(tx-d[coords[0]])*n[coords[0]] + -(ty-d[coords[1]])*n[coords[1]])/n[j] + d[j]
    t = np.zeros(3)
    t[coords[0]]=tx-d[coords[0]]
    t[coords[1]]=ty-d[coords[1]]
    t[j]=tz-d[j]
    #t = np.array([(tx-x),(ty-y),(tz-z)])

    modt = np.sqrt(t[0]**2 + t[1]**2 + t[2]**2)
    vhat = np.array([t[0]/modt, t[1]/modt, t[2]/modt])

    #set a sample veloctiy
    #v = vperi
    v = vhat*vperi
    return v
        #vall.append(v)
    #return vall


def get_impact_parameter(rp,xe,vxe,xyz,vel):
    sim2 = rebound.Simulation()
    sim2.units = ('yr', 'AU', 'Msun')


    me = 5.9724e24 * u.kg
    me = me.to(u.Msun).value
    mm = 7.34767309e22 * u.kg
    mm = mm.to(u.Msun).value
    sim2.add(m=me,x=xe[0],y=xe[1],z=xe[2],vx=vxe[0],vy=vxe[1],vz=vxe[2])
    sim2.add(m=0.0,x=xyz[0],y=xyz[1],z=xyz[2],vx=vel[0],vy=vel[1],vz=vel[2])
    orbit = sim2.particles[1].calculate_orbit(primary=sim2.particles[0])
    a = orbit.a
    del sim2.particles

    b = np.sqrt((rp - a)**2 - a**2)
    return b


def generate_unbound_orbits(BATCHSIZE,v_inf,G,xe,vxe,xm,vxm):
    rng = np.random.default_rng()
    Re = 6378.1 * u.km
    Re = Re.to(u.AU).value
    me = 5.9724e24 * u.kg
    me = me.to(u.Msun).value
    mm = 7.34767309e22 * u.kg
    mm = mm.to(u.Msun).value
    mb = me+mm
    am = 384748 * u.km
    am = am.to(u.AU).value

    count_unbound = 0
    i = 0
    xs_all = []
    xe_all = []
    xm_all = []
    vxs_all = []
    vxe_all = []
    vxm_all = []
    xyz_all = []
    vel_all = []
    rp_all = []
    b_all = []
    while count_unbound < BATCHSIZE:
        xb = (me*xe[i][0] + mm*xm[i][0]) / mb
        yb = (me*xe[i][1] + mm*xm[i][1]) / mb
        zb = (me*xe[i][2] + mm*xm[i][2]) / mb
        
        vxb = (me*vxe[i][0] + mm*vxm[i][0]) / mb
        vyb = (me*vxe[i][1] + mm*vxm[i][1]) / mb
        vzb = (me*vxe[i][2] + mm*vxm[i][2]) / mb
        
        xbi = np.array([xb,yb,zb])
        vxbi = np.array([vxb,vyb,vzb])
        xei = np.array(xe[i])
        xmi = np.array(xm[i])
        vxei = np.array(vxe[i])
        vxmi = np.array(vxm[i])
        rpos_e = xei-xbi
        rpos_m = xmi-xbi
        rvel_e = vxei-vxbi
        rvel_m = vxmi-vxbi
        dpos_e = np.sqrt(np.sum(rpos_e**2,axis=0))
        dpos_m = np.sqrt(np.sum(rpos_m**2,axis=0))
        am_mod = dpos_m
        Re_mod = Re-dpos_e
        rp = rng.random() * (am_mod-Re_mod) + Re_mod
        phi = (2*np.pi)*rng.random()
        theta = (np.pi)*rng.random()

        x = rp * np.cos(phi) * np.sin(theta)
        y = rp * np.sin(phi) * np.sin(theta)
        z = rp * np.cos(theta)
        xyz = np.array([x,y,z])
        vperi = velocity_at_peri(xe[i],xyz,v_inf,G,me)
        vel_part=get_velocity(vperi,x,y,z)
        
        
        xyz[0]+=xbi[0]
        xyz[1]+=xbi[1]
        xyz[2]+=xbi[2]

        vel_part[0]+=vxbi[0]
        vel_part[1]+=vxbi[1]
        vel_part[2]+=vxbi[2]
        xsi = np.array([0.0,0.0,0.0])
        #xmi = np.array([0.0,0.0,0.0])
        
        
        ene=get_energy(G,mb,0.0,xsi,xbi,xmi,xyz,vel_part)
        if(ene > 0):
            #xyz_all.append(xyz)
            #vel_all.append(vel_part)
            #rp_all.append(rp)

            simi = rebound.Simulation()
            simi.units = ('yr', 'AU', 'Msun')
            simi.add(m=1.0)
            simi.add(m=me, x=xbi[0],y=xbi[1],z=xbi[2],vx=vxbi[0],vy=vxbi[1],vz=vxbi[2])
            #simi.add(m=mm, x=xmi[0],y=xmi[1],z=xmi[2],vx=vxmi[0],vy=vxmi[1],vz=vxmi[2])
            simi.add(m=0.0,x=xyz[0],y=xyz[1],z=xyz[2],vx=vel_part[0],vy=vel_part[1],vz=vel_part[2])
            simi.move_to_com()
            simi.integrate(-2) #### IMPORTANT LINE!! HARDCODED ####
            sun = simi.particles[0]
            #earth = simi.particles[1]
            #moon = simi.particles[2]
            #comet = simi.particles[3]
            bary = simi.particles[1]
            comet=simi.particles[2]
            xbii = np.array([bary.x,bary.y,bary.z])
            xeii = np.array([bary.x+rpos_e[0],bary.y+rpos_e[1],bary.z+rpos_e[2]])
            xmii = np.array([bary.x+rpos_m[0],bary.y+rpos_m[1],bary.z+rpos_m[2]])
            veii = np.array([bary.vx+rvel_e[0],bary.vy+rvel_e[1],bary.vz+rvel_e[2]])
            vmii = np.array([bary.vx+rvel_m[0],bary.vy+rvel_m[1],bary.vz+rvel_m[2]])
            
            xsii = np.array([sun.x,sun.y,sun.z])
            #xeii = np.array([earth.x,earth.y,earth.z])
            #xmii = np.array([moon.x,moon.y,moon.z])
            xcii = np.array([comet.x,comet.y,comet.z])
            vsii = np.array([sun.vx,sun.vy,sun.vz])
            #veii = np.array([earth.vx,earth.vy,earth.vz])
            #vmii = np.array([moon.vx,moon.vy,moon.vz])
            vcii = np.array([comet.vx,comet.vy,comet.vz])
            enei = get_energy(G,me,mm,xsii,xeii,xmii,xcii,vcii)
            if(enei > 0):
                #simi.add(primary=simi.particles[1],m=mm,a=am,e=0.0549,f=rng.random() * 2 * np.pi)
                #moon = simi.particles[3]
                #xm_all.append(np.array([moon.x,moon.y,moon.z]))
                #vxm_all.append(np.array([moon.vx,moon.vy,moon.vz]))
                xm_all.append(xmii)
                vxm_all.append(vmii)
                xs_all.append(xsii)
                xe_all.append(xeii)
                vxs_all.append(vsii)
                vxe_all.append(veii)
                xyz_all.append(xcii)
                vel_all.append(vcii)
                rp_all.append(rp)
                b=get_impact_parameter(rp,xei,vxe[i],xyz,vel_part)
                b_all.append(b)
                i+=1
                count_unbound+=1
            simi = None
    xyz_all = np.array(xyz_all)
    vel_all = np.array(vel_all)
    rp_all = np.array(rp_all)

    return xs_all,xe_all,xm_all, xyz_all, vxs_all, vxe_all, vxm_all, vel_all, rp_all, b_all



def run_forward_ias15(xs,xe,xm,xc,vxs,vxe,vxm,vxc):
    sim = rebound.Simulation()
    sim.units = ('yr', 'AU', 'Msun')
    sim.integrator = "ias15"
#mjup = 0.0009547919099366768
    me = 5.9724e24 * u.kg
    me = me.to(u.Msun).value
    mm = 7.34767309e22 * u.kg
    mm = mm.to(u.Msun).value

    sim.add(m=1.0,x=xs[0],y=xs[1],z=xs[2],vx=vxs[0],vy=vxs[1],vz=vxs[2])
    sim.add(m=me,x=xe[0],y=xe[1],z=xe[2],vx=vxe[0],vy=vxe[1],vz=vxe[2])
    sim.add(m=me,x=xm[0],y=xm[1],z=xm[2],vx=vxm[0],vy=vxm[1],vz=vxm[2])
    sim.add(m=me,x=xc[0],y=xc[1],z=xc[2],vx=vxc[0],vy=vxc[1],vz=vxc[2])
    sim.move_to_com()
    sim.integrate(5) #placeholder value for now

    xsi = np.array([sim.particles[0].x,sim.particles[0].y,sim.particles[0].z])
    xei = np.array([sim.particles[1].x,sim.particles[1].y,sim.particles[1].z])
    xmi = np.array([sim.particles[2].x,sim.particles[2].y,sim.particles[2].z])
    xci = np.array([sim.particles[3].x,sim.particles[3].y,sim.particles[3].z])
    vci = np.array([sim.particles[3].vx,sim.particles[3].vy,sim.particles[3].vz])
    vsi = np.array([sim.particles[0].vx,sim.particles[0].vy,sim.particles[0].vz])
    vei = np.array([sim.particles[1].vx,sim.particles[1].vy,sim.particles[1].vz])
    vmi = np.array([sim.particles[2].vx,sim.particles[2].vy,sim.particles[2].vz])
    ene = get_energy(sim.G,me,mm,xsi,xei,xmi,xci,vci)
    isunbound=True
    if(ene < 0):
        isbound=False
    sim = None
    return ([isunbound,xsi,xei,xmi,xci,vsi,vei,vmi,vci])

def run_forward_kepler(xs,xe,xm,xc,vxs,vxe,vxm,vxc,t_end):
    xs_c = (c_double * 3)(*xs)
    xe_c = (c_double * 3)(*xe)
    xm_c = (c_double * 3)(*xm)
    xc_c = (c_double * 3)(*xc)

    vxs_c = (c_double * 3)(*vxs)
    vxe_c = (c_double * 3)(*vxe)
    vxm_c = (c_double * 3)(*vxm)
    vxc_c = (c_double * 3)(*vxc)

    collb_flag = c_int(0)

    libkepler.integrate(xs_c,xe_c,xm_c,xc_c,vxs_c,vxe_c,vxm_c,vxc_c,c_double(t_end),byref(collb_flag))


    #pack the data properly and return it


    #if(collb_flag.value == 4):
        #print("aye!")
    #    isunbound = False

    #DEBUG START
#     if(collb_flag.value == 0):
#         sim = rebound.Simulation()
#         sim = rebound.Simulation()
#         sim.units = ('yr', 'AU', 'Msun')

#         me = 5.9724e24 * u.kg
#         me = me.to(u.Msun).value
#         mm = 7.34767309e22 * u.kg
#         mm = mm.to(u.Msun).value

#         sim.add(m=1.0,x=xs[0],y=xs[1],z=xs[2],vx=vxs[0],vy=vxs[1],vz=vxs[2])
#         sim.add(m=me,x=xe[0],y=xe[1],z=xe[2],vx=vxe[0],vy=vxe[1],vz=vxe[2])
#         sim.add(m=me,x=xm[0],y=xm[1],z=xm[2],vx=vxm[0],vy=vxm[1],vz=vxm[2])
#         sim.add(m=0.0,x=xc[0],y=xc[1],z=xc[2],vx=vxc[0],vy=vxc[1],vz=vxc[2])
#         sim.move_to_com()
#         sim.integrate(t_end)
#         comet = sim.particles[3]
#         print('-------------------')
#         print((comet.x-xc_c[0])/comet.x,(comet.y-xc_c[1])/comet.y,(comet.z-xc_c[2])/comet.z)
#         print((comet.vx-vxc_c[0])/comet.vx,(comet.vy-vxc_c[1])/comet.vy,(comet.vz-vxc_c[2])/comet.vz)
#         #orbit = sim.particles[3].calculate_orbit(primary=sim.particles[0])
#         #print(orbit.a)
#         print('-------------------')
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





#@stopit_after_timeout(20, raise_exception=True)
def run_scattering_experiment(args):
    #set BATCHSIZE to be 100
    #generate BATCHSIZE random earth-moon posvels first
    v_a = args[0]
    BATCHSIZE=int(args[1])

    xe,xm,vxe,vxm = generate_em_posvel(BATCHSIZE)


    sim = rebound.Simulation()
    sim.units = ('yr', 'AU', 'Msun')

    me = 5.9724e24 * u.kg
    me = me.to(u.Msun).value
    mm = 7.34767309e22 * u.kg
    mm = mm.to(u.Msun).value

    Re = 6378.1 * u.km
    Re = Re.to(u.AU).value

    am = 384748 * u.km
    am = am.to(u.AU).value
    v_asym = v_a #km/s
    v_asym_kms = v_asym * (u.km/u.s)#km/s
    #print(v_asym_kms)
    v_asym_auyr = v_asym_kms.to(u.AU/u.yr).value
    s = np.zeros(BATCHSIZE) + v_asym

    #now run the simulation in two steps
    #1. run it backwards in time
    #2. if it is unbound, get a new random position for the moon
    xs_all, xe_all, xm_all, xyz_all, vxs_all, vxe_all, vxm_all, vel_all, rp_all, b_all = generate_unbound_orbits(BATCHSIZE,v_asym_auyr,sim.G,xe,vxe,xm,vxm)
    #3. run it forward
    rp_debug = []
    res_batch = []
    for i in range(0,BATCHSIZE):
#         xxb = [xb[i],yb[i],zb[i]]
#         vvxb = [vxb[i],vyb[i],vzb[i]]
#         xxc = [xyz_tang[i][0],xyz_tang[i][1],xyz_tang[i][2]]
#         vvxc = [vel_parts[i][0],vel_parts[i][1],vel_parts[i][2]]
#         cp = run_kepler_mod(xxb,xxc,vvxb,vvxc,t_end)

#         sim.add(m=1,x=0.0,y=0.0,z=0.0,vx=0.0,vy=0.0,vz=0.0)
#         sim.add(m=0.0,x=cp[0],y=cp[1],z=cp[2],vx=cp[3],vy=cp[4],vz=cp[5])
#         sim.move_to_com()
#         orbit = sim.particles[1].calculate_orbit(primary=sim.particles[0])
#         qi = (orbit.a*(1-orbit.e))
#         #rp_debug.append(qi)
#         t_end_s = 2*orbit.T
#         del sim.particles
#         xc = [cp[0],cp[1],cp[2]]
#         vxc = [cp[3],cp[4],cp[5]]
        t_end_s = 2 #### IMPORTANT!!! CHANGE IT ####
        res_s = run_forward_kepler(xs_all[i],xe_all[i],xm_all[i],xyz_all[i],vxs_all[i],vxe_all[i],vxm_all[i],vel_all[i],t_end_s)
        qi = rp_all[i]
        debug_e = np.array(res_s[1])
        debug_c = np.array(res_s[3])
        debug_d = np.sqrt(np.sum((debug_c-debug_e)**2,axis=0))
        rp_debug.append(debug_d)
        if(res_s[-1]>0):
            #xs,xe,xm,xc,vxs,vxe,vxm,vxc,initial_b,initial_rp,collb_flag
            res_batch.append([res_s[0],res_s[1],res_s[2],res_s[3],res_s[4],res_s[5],res_s[6],res_s[7],b_all[i],qi,res_s[-1]])
    #print("finished batch!")
    print(rp_all/am)
    print(rp_debug/am)
    print(rp_all/am - rp_debug/am)
    #print(((np.array(rp_debug)-rp_all)/rp_all))
    return res_batch
    #print((np.array(rp_debug)-np.array(rp_debug2))/np.array(rp_debug))
    #print(rp_debug)
if __name__ == "__main__":
    #multiprocessing.set_start_method('forkserver')
    BATCHSIZE = 100
    v_a = 1.0
    args = (v_a,BATCHSIZE)
    res = run_scattering_experiment(args)
    print(len(res))
    #NBATCHES = 2000
#     BATCHSIZE,NBATCHES,v_a = get_arguments()
#     args = (v_a, BATCHSIZE)
#     #run_scattering_experiment(args)
# #     nsample, v_a = get_arguments()
#     start_time = time.time()
#     res_data = []
#     count = 0
#     count_capbound = 0

#     print("RAPTOR v1.0. Made by Diptajyoti Mukherjee. Jan 13, 2023",flush=True)
#     print("System: Sun-Earth-Moon",flush=True)
#     print("BATCHSIZE: ",BATCHSIZE,flush=True)
#     print("NBATCHES: ", NBATCHES,flush=True)
#     print("Total experiments: ", BATCHSIZE*NBATCHES,flush=True)
#     print("V_asym (km/s): ", v_a,flush=True)

#     sim = rebound.Simulation()
#     sim.units = ('yr', 'AU', 'Msun')

#     me = 5.9724e24 * u.kg
#     me = me.to(u.Msun).value
#     mm = 7.34767309e22 * u.kg
#     mm = mm.to(u.Msun).value

#     v_asym = v_a #km/s
#     v_asym_kms = v_asym * (u.km/u.s)#km/s
#     #print(v_asym_kms)
#     v_asym_auyr = v_asym_kms.to(u.AU/u.yr).value
#     upper_b = np.sqrt(1 + 2 * sim.G*(me+mm+1)/(v_asym_auyr**2))
#     print("Max b (AU): {:21.16f}".format(upper_b),flush=True)
#     #print("generated earth, moon pos")
#     sim = None

#     UP = "\x1B[3A"
#     CLR = "\x1B[0K"
#     print("\n\n")

#     #for debug purposes only
#     rrr=[]
#     rrre=[]

#     #crude solution to take care of RAM overflows
#     #divide the experiments into batches
#     #currently hardcoded. fix it in the future
#     NDIVISIONS=100
#     DIVISION_NBATCH = int(NBATCHES/NDIVISIONS)
#     args_rep = ((args, ) * DIVISION_NBATCH)
#     for j in range(0,NDIVISIONS):

#         #for saving the data
#         xs = []
#         ys = []
#         zs = []
#         vxs = []
#         vys = []
#         vzs = []

#         xe = []
#         ye = []
#         ze = []
#         vxe = []
#         vye = []
#         vze = []

#         xm = []
#         ym = []
#         zm = []
#         vxm = []
#         vym = []
#         vzm = []

#         xc = []
#         yc = []
#         zc = []
#         vxc = []
#         vyc = []
#         vzc = []

#         q = []
#         b = []
#         collb_flag = []

#         count_division = 0


#         with multiprocessing.Pool() as pool:
#             iterator = pool.imap_unordered(run_scattering_experiment, args_rep)

#             while True:
#                 try:
#                     count+=1
#                     count_division+=1
#                     results = iterator.next()
#                     count_capbound += len(results)
#                     processed_percent = count/NBATCHES * 100
#                     print(f"{UP}processed batch: {processed_percent}%{CLR}\ncaptured/bound: {count_capbound}{CLR}\n",flush=True)

#                     if(len(results)>0):
#                         for result in results:
#                             xs.append(result[0][0])
#                             ys.append(result[0][1])
#                             zs.append(result[0][2])

#                             xe.append(result[1][0])
#                             ye.append(result[1][1])
#                             ze.append(result[1][2])

#                             xm.append(result[2][0])
#                             ym.append(result[2][1])
#                             zm.append(result[2][2])

#                             xc.append(result[3][0])
#                             yc.append(result[3][1])
#                             zc.append(result[3][2])

#                             vxs.append(result[4][0])
#                             vys.append(result[4][1])
#                             vzs.append(result[4][2])

#                             vxe.append(result[5][0])
#                             vye.append(result[5][1])
#                             vze.append(result[5][2])

#                             vxm.append(result[6][0])
#                             vym.append(result[6][1])
#                             vzm.append(result[6][2])

#                             vxc.append(result[7][0])
#                             vyc.append(result[7][1])
#                             vzc.append(result[7][2])

#                             b.append(result[8])
#                             q.append(result[9])
#                             collb_flag.append(result[10])

#                     #print(len(result))
#                     #rrr.append(result)

#                 except StopIteration:
#                     break
#                 except Exception as e:
#                     # do something
#                     rrre.append(e)
#             pool.terminate()

#         count-=1
#         count_division-=1
#         #print("\n\n")
#         #print(rrre)
#         #print("writing to file...")

#         hf = h5py.File('captured_bound_data_{}.h5'.format(j), 'w')
#         hf.create_dataset('xs',data=xs)
#         hf.create_dataset('ys',data=ys)
#         hf.create_dataset('zs',data=zs)
#         hf.create_dataset('vxs',data=vxs)
#         hf.create_dataset('vys',data=vys)
#         hf.create_dataset('vzs',data=vzs)


#         hf.create_dataset('xe',data=xe)
#         hf.create_dataset('ye',data=ye)
#         hf.create_dataset('ze',data=ze)
#         hf.create_dataset('xm',data=xm)
#         hf.create_dataset('ym',data=ym)
#         hf.create_dataset('zm',data=zm)

#         hf.create_dataset('vxe',data=vxe)
#         hf.create_dataset('vye',data=vye)
#         hf.create_dataset('vze',data=vze)
#         hf.create_dataset('vxm',data=vxm)
#         hf.create_dataset('vym',data=vym)
#         hf.create_dataset('vzm',data=vzm)

#         hf.create_dataset('xc',data=xc)
#         hf.create_dataset('yc',data=yc)
#         hf.create_dataset('zc',data=zc)
#         hf.create_dataset('vxc',data=vxc)
#         hf.create_dataset('vyc',data=vyc)
#         hf.create_dataset('vzc',data=vzc)

#         hf.create_dataset('b',data=b)
#         hf.create_dataset('q',data=q)
#         hf.create_dataset('collb_flag',data=collb_flag)

#         info_dset = hf.create_dataset("info",dtype="f8")
#         info_dset.attrs.create("v_asym",data=v_a)
#         info_dset.attrs.create("upper_b",data=upper_b)
#         info_dset.attrs.create("total_experiments",data=float(count_division*BATCHSIZE))
#         hf.close()
#     print("\n\n")
#     print(rrre)
#     print("total time taken: ",time.time()-start_time)

#     print("writing file completed")
