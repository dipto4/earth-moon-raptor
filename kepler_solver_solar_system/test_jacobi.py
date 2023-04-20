import numpy as np


def cart2jacobi(mass,x,v):
    eta = np.cumsum(mass)
    #x_star = []
    #v_star = []
    x_star = np.zeros((len(mass),3))
    v_star = np.zeros((len(mass),3))
    rcm = mass[0] * x[0]
    vcm = mass[0] * v[0]

    for i in range(1,len(mass)):
        xi_star = x[i] - rcm / eta[i-1]
        vi_star = v[i] - vcm / eta[i-1]
        rcm = rcm*(1+mass[i]/eta[i-1]) + mass[i] * xi_star
        vcm = vcm*(1+mass[i]/eta[i-1]) + mass[i] * vi_star
        x_star[i] = xi_star
        v_star[i] = vi_star

    x_star[0] = rcm / eta[-1]
    v_star[0] = vcm / eta[-1]

    return eta, x_star, v_star, rcm,vcm
    pass



def jacobi2cart(eta,x_star,v_star,mass):
    xyz = np.zeros((len(eta),3))
    vel = np.zeros((len(eta),3))
    rcm = eta[-1] * x_star[0]
    vcm = eta[-1] * v_star[0]
    for i in range(len(eta)-1,0,-1):
        rcm = (rcm - mass[i]*x_star[i]) / eta[i]
        vcm = (vcm - mass[i]*v_star[i]) / eta[i]
        xyz[i] = x_star[i] + rcm
        vel[i] = v_star[i] + vcm
        rcm = eta[i-1]*rcm
        vcm = eta[i-1]*vcm

    xyz[0] = rcm / mass[0]
    vel[0] = vcm / mass[0]

    return xyz,vel
    pass


if __name__ == "__main__":
    #for testing purposes
    #mass = np.array([5.,4.,3.,2.])
    #x = np.array([0.,1.,2.,4.])
    #y = np.array([0.,1.,-1.,3.])
    #z = np.array([0.,0.,0.,0.])

    #vx = np.array([0.,1.,-1.,-1.])
    #vy = np.array([0.,2.,0.,2.])
    #vz = np.array([0.,0.,0.,0.])
    data = np.loadtxt("solar_system.dat")
    mass = data[:,0]
    x = data[:,1]
    y = data[:,2]
    z = data[:,3]
    vx = data[:,4]
    vy = data[:,5]
    vz = data[:,6]
    xyz = np.array((x,y,z))

    xyz = xyz.T
    vel = np.array((vx,vy,vz))
    vel = vel.T

    eta, x_star, v_star, rcm, vcm = cart2jacobi(mass,xyz,vel)
    #print(eta)
    print(x_star)
    print(v_star)

    xyz_r, vel_r = jacobi2cart(eta,x_star,v_star,mass)
    print(xyz_r)
    print(vel_r)
