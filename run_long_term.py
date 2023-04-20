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
libkepler = cdll.LoadLibrary('./solarsystemPy.so')

libkepler.integrate.argtypes = [POINTER(c_double), POINTER(c_double), POINTER(c_double),
                                POINTER(c_double), POINTER(c_double), POINTER(c_double),
                                POINTER(c_double), POINTER(c_double),c_int, c_double,POINTER(c_int)]



planet_data = np.loadtxt("planet_data")
planet_names = ["Venus","Earth","Mars","Jupiter","Saturn","Uranus","Neptune"]
semi_planets = {}
mass_planets = {}
ecc_planets = {}
radii_planets = {}
for kk in range(0,len(planet_names)):
    mass_planets[planet_names[kk]] = planet_data[kk][0]
    semi_planets[planet_names[kk]] = planet_data[kk][1]
    ecc_planets[planet_names[kk]] = planet_data[kk][2]
    radii_planets[planet_names[kk]] = planet_data[kk][3]

radii_planets["Sun"] = 0.00465047 # sorry I am too lazy to change the name so the Sun is a planet :)
mass_planets["Sun"] = 1.0

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


def generate_planet_posvel():
    sim = rebound.Simulation()
    sim.units = ('yr', 'AU', 'Msun')
    rng = np.random.default_rng()
    sim.add(m=1.0)

    sim.add(primary=sim.particles[0],m=mass_planets["Venus"],a=semi_planets["Venus"],e=ecc_planets["Venus"],f=rng.random() * 2 * np.pi)
    sim.add(primary=sim.particles[0],m=mass_planets["Mars"],a=semi_planets["Mars"],e=ecc_planets["Mars"],f=rng.random() * 2 * np.pi)
    sim.add(primary=sim.particles[0],m=mass_planets["Jupiter"],a=semi_planets["Jupiter"],e=ecc_planets["Jupiter"],f=rng.random() * 2 * np.pi)
    sim.add(primary=sim.particles[0],m=mass_planets["Saturn"],a=semi_planets["Saturn"],e=ecc_planets["Saturn"],f=rng.random() * 2 * np.pi)
    sim.add(primary=sim.particles[0],m=mass_planets["Uranus"],a=semi_planets["Uranus"],e=ecc_planets["Uranus"],f=rng.random() * 2 * np.pi)
    sim.add(primary=sim.particles[0],m=mass_planets["Neptune"],a=semi_planets["Neptune"],e=ecc_planets["Neptune"],f=rng.random() * 2 * np.pi)
    x_p = []
    y_p = []
    z_p = []
    vx_p = []
    vy_p = []
    vz_p = []
    for i in range(1,7):
        x_p.append(sim.particles[i].x)
        y_p.append(sim.particles[i].y)
        z_p.append(sim.particles[i].z)
        vx_p.append(sim.particles[i].vx)
        vy_p.append(sim.particles[i].vy)
        vz_p.append(sim.particles[i].vz)

    return x_p, y_p, z_p, vx_p, vy_p, vz_p

#this function is going to be different for earth-moon vs. jupiter
def clean_data(xs,ys,zs,vxs,vys,vzs,
               xe,ye,ze,vxe,vye,vze,
               xm,ym,zm,vxm,vym,vzm,
               x_p,y_p,z_p,vx_p,vy_p,vz_p):

    me = 5.9724e24 * u.kg
    me = me.to(u.Msun).value
    mm = 7.34767309e22 * u.kg
    mm = mm.to(u.Msun).value
    mtot = me + mm
    #compute earth-moon barycenter
    xb = (me*xe + mm *xm) / mtot
    yb = (me*ye + mm *ym) / mtot
    zb = (me*ze + mm *zm) / mtot

    vxb = (me*vxe + mm *vxm) / mtot
    vyb = (me*vye + mm *vym) / mtot
    vzb = (me*vze + mm *vzm) / mtot

    #converting to heliocentric
    xb -= xs
    yb -= ys
    zb -= zs
    vxb -= vxs
    vyb -= vys
    vzb -= vzs

    #output the final position, vel, radii arrays including sun
    #has to be in the hierarchical order (for now)
    #1. Sun
    #2. Venus
    #3. Earth
    #4. Mars
    #5. Jupiter
    #6. Saturn
    #7. Uranus
    #8. Neptune

    x_all = [0.0,x_p[0],xb,x_p[1],x_p[2],x_p[3],x_p[4],x_p[5]]
    y_all = [0.0,y_p[0],yb,y_p[1],y_p[2],y_p[3],y_p[4],y_p[5]]
    z_all = [0.0,z_p[0],zb,z_p[1],z_p[2],z_p[3],z_p[4],z_p[5]]
    vx_all = [0.0,vx_p[0],vxb,vx_p[1],vx_p[2],vx_p[3],vx_p[4],vx_p[5]]
    vy_all = [0.0,vy_p[0],vyb,vy_p[1],vy_p[2],vy_p[3],vy_p[4],vy_p[5]]
    vz_all = [0.0,vz_p[0],vzb,vz_p[1],vz_p[2],vz_p[3],vz_p[4],vz_p[5]]
    object_names = ["Sun","Venus","Earth","Mars","Jupiter","Saturn","Uranus","Neptune"]
    r_all = []
    m_all = []
    for obj in object_names:
        r_all.append(radii_planets[obj])
        m_all.append(mass_planets[obj])

    m_all[2] = mtot #use the barycenter mass calculated from the data


    return x_all,y_all,z_all,vx_all,vy_all,vz_all,m_all,r_all


def integrate_system(x,y,z,vx,vy,vz,m,r,t_end):
    N = len(x)
    x_c = (c_double * N) (*x)
    y_c = (c_double * N) (*y)
    z_c = (c_double * N) (*z)
    vx_c = (c_double * N) (*vx)
    vy_c = (c_double * N) (*vy)
    vz_c = (c_double * N) (*vz)
    m_c = (c_double * N) (*m)
    r_c = (c_double * N) (*r)
    collb_flag = c_int(0)
    libkepler.integrate(x_c,y_c,z_c,vx_c,vy_c,vz_c,m_c,r_c,c_int(N),c_double(t_end),byref(collb_flag))
    for i in range(0,N):
        x[i] = x_c[i]
        y[i] = y_c[i]
        z[i] = z_c[i]
        vx[i] = vx_c[i]
        vy[i] = vy_c[i]
        vz[i] = vz_c[i]

    return collb_flag.value


@stopit_after_timeout(1*60*60, raise_exception=True)
def run_scattering_experiment(args):
    xs = args[0]
    ys = args[1]
    zs = args[2]
    vxs = args[3]
    vys = args[4]
    vzs = args[5]

    xe = args[6]
    ye = args[7]
    ze = args[8]
    vxe = args[9]
    vye = args[10]
    vze = args[11]

    xm = args[12]
    ym = args[13]
    zm = args[14]
    vxm = args[15]
    vym = args[16]
    vzm = args[17]

    xc = args[18]
    yc = args[19]
    zc = args[20]
    vxc = args[21]
    vyc = args[22]
    vzc = args[23]

    x_p, y_p, z_p, vx_p, vy_p, vz_p = generate_planet_posvel()
    x_p, y_p, z_p, vx_p, vy_p, vz_p, m_p, r_p = clean_data(xs,ys,zs,vxs,vys,vzs,
                                                 xe,ye,ze,vxe,vye,vze,
                                                 xm,ym,zm,vxm,vym,vzm,
                                                 x_p,y_p,z_p,vx_p,vy_p,vz_p)

    #add the comet now
    x_p.append(xc-xs)
    y_p.append(yc-ys)
    z_p.append(zc-zs)
    vx_p.append(vxc-vxs)
    vy_p.append(vyc-vys)
    vz_p.append(vzc-vzs)
    m_p.append(0.0)
    r_p.append(0.0)
        
    #DEBUG
    #x_p.append(x_p[-1])
    #y_p.append(y_p[-1])
    #z_p.append(z_p[-1])
    #vx_p.append(vx_p[-1])
    #vy_p.append(vy_p[-1])
    #vz_p.append(vz_p[-1])
    #m_p.append(0.0)
    #r_p.append(0.0)

        
        
    t_end = 1e3 #this is the length of each interval
    t_total = 1e7 #this is the total length
    t_current = 0.0
    t = [0.0]
    sim = rebound.Simulation()
    sim.units = ('yr','AU','Msun')
    sim.add(m=1.0)
    sim.add(m=0.0,x=x_p[-1],y=y_p[-1],z=z_p[-1],vx=vx_p[-1],vy=vy_p[-1],vz=vz_p[-1])

    orbit = sim.particles[1].calculate_orbit(primary=sim.particles[0])
    a = [orbit.a]
    ecc = [orbit.e]
    inc = [orbit.inc*180./np.pi]
    del sim.particles
    flag = 0

    while(t_current < t_total):
        flag = integrate_system(x_p,y_p,z_p,vx_p,vy_p,vz_p,m_p,r_p,t_end)
        t_current += t_end
        t.append(t_current)
        sim.add(m=1.0)
        sim.add(m=0.0,x=x_p[-1],y=y_p[-1],z=z_p[-1],vx=vx_p[-1],vy=vy_p[-1],vz=vz_p[-1])
        orbit = sim.particles[1].calculate_orbit(primary=sim.particles[0])
        a.append(orbit.a)
        ecc.append(orbit.e)
        inc.append(orbit.inc*180./np.pi)
        if(orbit.a >= 1e3):
            flag = 5 #use this value to indicate that the particle was not lost but semi-major axis went above 1000 au
        del sim.particles
        #print(t_current, flag, orbit.a, orbit.e) 
        if(flag > 0):
            #print("sim ended!")
            #print(flag)
            break
            
    N_datapoints = int(t_total) // int(t_end) + 1
    if (len(a) <= N_datapoints):
        a += [-999] * (N_datapoints - len(a))
        ecc += [-999] * (N_datapoints - len(ecc))
        inc += [-999] * (N_datapoints - len(inc))

    return [a,ecc,inc,flag]

def read_file(filename):
    data = h5py.File(filename,'r')
    xs = np.array(data['xs'])
    ys = np.array(data['ys'])
    zs = np.array(data['zs'])
    vxs = np.array(data['vxs'])
    vys = np.array(data['vys'])
    vzs = np.array(data['vzs'])

    xe = np.array(data['xe'])
    ye = np.array(data['ye'])
    ze = np.array(data['ze'])
    vxe = np.array(data['vxe'])
    vye = np.array(data['vye'])
    vze = np.array(data['vze'])

    xm = np.array(data['xm'])
    ym = np.array(data['ym'])
    zm = np.array(data['zm'])
    vxm = np.array(data['vxm'])
    vym = np.array(data['vym'])
    vzm = np.array(data['vzm'])

    xc = np.array(data['xc'])
    yc = np.array(data['yc'])
    zc = np.array(data['zc'])
    vxc = np.array(data['vxc'])
    vyc = np.array(data['vyc'])
    vzc = np.array(data['vzc'])

    args = []
    for i in range(0,len(xs)):
        arg = (xs[i],ys[i],zs[i],vxs[i],vys[i],vzs[i],
               xe[i],ye[i],ze[i],vxe[i],vye[i],vze[i],
               xm[i],ym[i],zm[i],vxm[i],vym[i],vzm[i],
               xc[i],yc[i],zc[i],vxc[i],vyc[i],vzc[i])
        args.append(arg)
    return args

if __name__ == "__main__":
    #multiprocessing.set_start_method('forkserver')
    #BATCHSIZE = 1000
    #v_a = 1.0
    #args = (v_a,BATCHSIZE)
    #NBATCHES = 2000
    #start_time = time.time()
    filename = "long_term_ic.h5"
    args = read_file(filename)#(v_a, BATCHSIZE)
    #run_scattering_experiment(args[-1])
    #print(time.time()-start_time)
#     nsample, v_a = get_arguments()
    start_time = time.time()
    res_data = []
    count = 0
    count_capbound = 0

    print("RAPTOR-long v1.0. Made by Diptajyoti Mukherjee. Apr 19, 2023",flush=True)
    print("Total experiments: ", len(args),flush=True)

    UP = "\x1B[3A"
    CLR = "\x1B[0K"
    print("\n\n")

    #for debug purposes only
    rrr=[]
    rrre=[]

    #crude solution to take care of RAM overflows
    #divide the experiments into batches
    #currently hardcoded. fix it in the future
    NDIVISIONS=1
    
    #DIVISION_NBATCH = int(NBATCHES/NDIVISIONS)
    #args_rep = ((args, ) * DIVISION_NBATCH)
    for j in range(0,NDIVISIONS):
        hf = h5py.File('long_term_data_{}.h5'.format(j), 'w')

        #a = []
        #ecc = []
        #inc = []
        #collb_flag = []

        count_division = 0


        with multiprocessing.Pool() as pool:
            iterator = pool.imap_unordered(run_scattering_experiment, args)

            while True:
                try:
                    count+=1
                    results = iterator.next()
                    #count_capbound += len(results)
                    processed_percent = count/len(args) * 100
                    print(f"{UP}processed batch: {processed_percent}%{CLR}\n{CLR}\n",flush=True)

                    if(True):
                        part_dset = hf.create_dataset('part_{}'.format(count_division),
                                                      data=np.c_[results[0],results[1],results[2]])
                        part_dset.attrs.create("collb_flag",data=results[3])
                        hf.flush()
                        count_division+=1

                        #a.append(results[0])
                        #ecc.append(results[1])
                        #inc.append(results[2])
                        #collb_flag.append(results[3])

                    #print(len(result))
                    #rrr.append(result)

                except StopIteration:
                    break
                except Exception as e:
                    # do something
                    rrre.append(e)
            pool.terminate()

        count-=1
        count_division-=1
        #print("\n\n")
        #print(rrre)
        #print("writing to file...")
        
        
#         hf = h5py.File('long_term_data_{}.h5'.format(j), 'w')
#         hf.create_dataset('a',data=a)
#         hf.create_dataset('ecc',data=ecc)
#         hf.create_dataset('inc',data=inc)
#         hf.create_dataset('collb_flag',data=collb_flag)

        hf.close()
    print("\n\n")
    print(rrre)
    print("total time taken: ",time.time()-start_time)

    print("writing file completed")
