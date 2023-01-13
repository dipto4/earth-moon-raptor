#include "universal.h"
#define M_EARTH 3.0036060914862233e-6
#define M_MOON 3.695250762067745e-8
#define M_SUN 1.0
#define R_SUN 0.00465047
#define R_EARTH 0.0000426342966684
#define R_MOON 0.0000116184808779
#define VSCALE 6.2830666414874994
#define ETA_SYM 0.00625
#define ETA 0.0001

#define max(a,b) \
    ({ __typeof__ (a) _a = (a); \
     __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

#define min(a,b) \
    ({ __typeof__ (a) _a = (a); \
     __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })

void force_and_jerk(double pos[2][3], double vel[2][3], double *tpos, double *tvel,
        double *ac, double *jc) {

    double accx = 0.0;
    double accy = 0.0;
    double accz = 0.0;

    double jerkx = 0.0;
    double jerky = 0.0;
    double jerkz = 0.0;

    double m = M_SUN + M_EARTH + M_MOON;
    double dx = tpos[0] - pos[0][0];
    double dy = tpos[1] - pos[0][1];
    double dz = tpos[2] - pos[0][2];
    double dist = sqrt(dx*dx + dy*dy + dz*dz );

    double dist2 = dist*dist;
    double dist3 = dist*dist*dist;

    double mass_over_dist3 = m/dist3;
    accx += -mass_over_dist3*dx;
    accy += -mass_over_dist3*dy;
    accz += -mass_over_dist3*dz;

    double dvx = tvel[0] - vel[0][0];
    double dvy = tvel[1] - vel[0][1];
    double dvz = tvel[2] - vel[0][2];
    double rrdot = (dx*dvx+dy*dvy+dz*dvz);


    jerkx += -m*(dvx/dist3 - 3*rrdot*dx/(dist3*dist2));
    jerky += -m*(dvy/dist3 - 3*rrdot*dy/(dist3*dist2));
    jerkz += -m*(dvz/dist3 - 3*rrdot*dz/(dist3*dist2));


    ac[0] = accx;
    ac[1] = accy;
    ac[2] = accz;
    jc[0] = jerkx;
    jc[1] = jerky;
    jc[2] = jerkz;

}

void timestep(double pos[2][3], double vel[2][3], double *dt) {
    double stepi;
    stepi = 1.e38;
    double mass = M_SUN + M_EARTH + M_MOON; 
    
    double dx = pos[1][0] - pos[0][0];
    double dy = pos[1][1] - pos[0][1];
    double dz = pos[1][2] - pos[0][2];
    double r2 = dx*dx + dy*dy + dz*dz;
    double r = sqrt(r2);

    double dvx = vel[1][0] - vel[0][0];
    double dvy = vel[1][1] - vel[0][1];
    double dvz = vel[1][2] - vel[0][2];

    double v2 = dvx * dvx + dvy * dvy + dvz * dvz ; 


    double vdotdr2 = (dx*dvx+dy*dvy+dz*dvz) / r2;

    double tau = ETA_SYM / M_SQRT2 * r * sqrt(r/(mass));
    double dtau = 3 * tau * vdotdr2 / 2;

    dtau = (dtau < 1) ? dtau : 1;
    tau /= (1 - dtau / 2);

    stepi = min(tau, stepi);
    tau = ETA_SYM * r / sqrt(v2);
    dtau = tau * vdotdr2 * (1 + (mass) / (v2 * r));
    dtau = (dtau < 1)? dtau: 1;

    tau /= (1 - dtau / 2);
    stepi = min(tau, stepi);

    *dt = stepi;
}

void initial(double pos[2][3],double vel[2][3], double *ac, double *jc, double* a2c, double *a3c, double *dt)  {
    double tpos[3] = {pos[1][0],pos[1][1],pos[1][2]};
    double tvel[3] = {vel[1][0],vel[1][1],vel[1][2]};
    force_and_jerk(pos,vel,tpos,tvel,ac,jc);
    timestep(pos,vel,dt);
    for(int d=0;d<3;d++) {
        a2c[d]=0;
        a3c[d]=0;
    }
}

void predictor(double pos[2][3], double vel[2][3], double *ac, double *jc, double *a2c, double *a3c, double dt,
        double *tpos, double *tvel) {
    
    double dt2_over_2 = dt*dt/2;
    double dt3_over_6 = dt*dt*dt/6;
    double dt4_over_24 = dt3_over_6*dt/4;
    double dt5_over_120 = dt4_over_24*dt/5;
    for(int d=0;d<3;d++) {
        tpos[d] = pos[1][d] + vel[1][d]*dt + ac[d]*dt2_over_2 + jc[d]*dt3_over_6 + a2c[d] * dt4_over_24 + a3c[d] * dt5_over_120;
        tvel[d] = vel[1][d] + ac[d]*dt + jc[d]*dt2_over_2 + a2c[d] * dt3_over_6 + a3c[d] * dt4_over_24;
    }
}

void corrector (double pos[2][3], double vel[2][3], double *ac, double *jc, double* a2c, double* a3c,
        double *tpos, double *tvel, double *tposc, double *tvelc,double *aci, double *jci, double *a2ci, double *a3ci,double dt) {
    double a2i[3], a3i[3];
    double step_i2 = dt * dt;
    double step_i3 =dt * dt * dt;

    for(int d=0;d<3;d++){
        a2i[d] = 2*(-3*(ac[d]-aci[d]) - dt*(2*jc[d] + 1*jci[d]))/(step_i2);
        a3i[d] = 6*(2*(ac[d]-aci[d]) + dt*(jc[d] + jci[d]))/(step_i3);
    }

    for(int d=0;d<3;d++) {
        tposc[d] = tpos[d] + (step_i2*step_i2)*(a2i[d]-a2c[d])/24 + (step_i2*step_i3)*(a3i[d]-a3c[d])/120;
        tvelc[d] = tvel[d] + (step_i3)*(a2i[d]-a2c[d])/6 + (step_i2*step_i2)*(a3i[d]-a3c[d])/24;
    }
    for(int d=0;d<3;d++) {
        a2ci[d] = a2i[d];
        a3ci[d] = a3i[d];
    }

}

void timestep_aarseth(double *aci, double *jci, double *a2ci, double *a3ci, double *dt) {
     double mod_f, mod_f1, mod_f2, mod_f3;
     double temp_t_i;
    double old = *dt;
     mod_f = sqrt(aci[0] * aci[0] + aci[1] * aci[1] + aci[2] * aci[2]);
     mod_f1 = sqrt(jci[0] * jci[0] + jci[1] * jci[1] + jci[2] * jci[2]);
     mod_f2 = sqrt(a2ci[0] * a2ci[0] + a2ci[1] * a2ci[1] + a2ci[2] * a2ci[2]);
     mod_f3 = sqrt(a3ci[0] * a3ci[0] + a3ci[1] * a3ci[1] + a3ci[2] * a3ci[2]);

     //find an optimal timestep
     temp_t_i = sqrt(ETA* ((mod_f * mod_f2) + mod_f1*mod_f1 ) / (mod_f1*mod_f3 + mod_f2*mod_f2));

     //check if timestep is not zero

     *dt = min((1.2*old), temp_t_i);
}

void evolve_hermite(double pos[2][3], double vel[2][3], double *ac, double *jc, double *a2c, double *a3c, double *dt) {
    double tpos[3], tvel[3];
    double tposc[3], tvelc[3];
    double aci[3], jci[3];
    double a2ci[3], a3ci[3];
    predictor(pos,vel,ac,jc,a2c,a3c,*dt,tpos,tvel);

    force_and_jerk(pos,vel,tpos,tvel,aci,jci);
    corrector(pos,vel,ac,jc,a2c,a3c,tpos,tvel,tposc,tvelc,aci,jci,a2ci,a3ci,*dt);
    
    //force_and_jerk(pos,vel,tposc,tvelc,aci,jci);
    //corrector(pos,vel,ac,jc,a2c,a3c,tpos,tvel,tposc,tvelc,aci,jci,a2ci,a3ci,*dt);
    

    //force_and_jerk(pos,vel,tposc,tvelc,aci,jci);
    //corrector(pos,vel,ac,jc,a2c,a3c,tpos,tvel,tposc,tvelc,aci,jci,a2ci,a3ci,dt);
 
    //force_and_jerk(pos,vel,tposc,tvelc,aci,jci);
    //corrector(pos,vel,ac,jc,a2c,a3c,tpos,tvel,tposc,tvelc,aci,jci,a2ci,a3ci,*dt);
 


    /*force_and_jerk(pos,vel,tposc,tvelc,aci,jci);
    corrector(pos,vel,ac,jc,tpos,tvel,tposc,tvelc,aci,jci,dt);
    
    force_and_jerk(pos,vel,tposc,tvelc,aci,jci);
    corrector(pos,vel,ac,jc,tpos,tvel,tposc,tvelc,aci,jci,dt);

    force_and_jerk(pos,vel,tposc,tvelc,aci,jci);
    corrector(pos,vel,ac,jc,tpos,tvel,tposc,tvelc,aci,jci,dt);
        */
    /*for(int i=0;i<2;i++) {
        force_and_jerk(pos,vel,tpos,tvel,aci,jci);
        corrector(pos,vel,ac,jc,tpos,tvel,aci,jci,dt);
    }*/
    for(int d=0;d<3;d++) {
        pos[1][d] = tposc[d];
        vel[1][d] = tvelc[d];
        ac[d] = aci[d];
        jc[d] = jci[d];
        a2ci[d] = a2ci[d] + *dt * a3ci[d];
        a3ci[d] =  a3ci[d];
    }
    timestep_aarseth(aci, jci, a2ci, a3ci, dt);

}

int heartbeat(double pos[2][3]) {
    int flag = 0;
    double dx = pos[1][0] - pos[0][0];
    double dy = pos[1][1] - pos[0][1];
    double dz = pos[1][2] - pos[0][2];
    double d = sqrt(dx*dx + dy*dy + dz*dz);
    if(d <= 20)
        flag = 1;
    return flag;
}

double getEne(double pos[2][3], double vel[2][3]) {
    double v2 = vel[1][0] * vel[1][0] + vel[1][1] * vel[1][1] + vel[1][2] * vel[1][2];
    double dx = pos[1][0] - pos[0][0];
    double dy = pos[1][1] - pos[0][1];
    double dz = pos[1][2] - pos[0][2];
    double d = sqrt(dx*dx + dy*dy + dz*dz);
    return (0.5*v2 - (M_SUN+M_MOON+M_EARTH)/d);
}

void evolve_system(double pos[2][3], double vel[2][3], double *ac, double *jc, double *a2c, double *a3c,double *dt, double *t, double t_end) {
    double t_cur = *t;
    //printf("%.7e %.7e\n",t_cur, *dt);
    size_t counter=0;
    double e0 = getEne(pos,vel);
    int flag = 0;
    while(t_cur < t_end) {
        counter++;
        if((t_cur+*dt) > t_end)
            *dt  = t_end-t_cur;
        t_cur += *dt;

        evolve_hermite(pos,vel,ac,jc,a2c,a3c,dt);
        flag = heartbeat(pos);
        if(flag==1)
            break;
        //debug_print(pos,vel);
        //timestep(pos,vel,dt);
        //double e = getEne(pos,vel);
        //printf("%.16e\n",(e-e0)/e0);
    }
    //printf("t_cur:%.16e\n",t_cur/VSCALE);
    //printf("steps:%i\n",counter);

}

double scale2nbody(double vel[2][3]) {
    for(size_t d=0; d<3;d++) {
        vel[0][d] = vel[0][d] / VSCALE;
        vel[1][d] = vel[1][d] / VSCALE;
    }
}
double scale2physical(double vel[2][3]) {
    for(size_t d=0; d<3;d++) {
        vel[0][d] = vel[0][d] * VSCALE;
        vel[1][d] = vel[1][d] * VSCALE;
    }
}


void integrate2body(double *xb, double *xc, double *vxb, double *vxc, double t_end, double *xci, double *vxci) {
    double comx[3], comv[3];
    double mb = M_SUN + M_EARTH + M_MOON;
    double ac[3], jc[3], a2c[3], a3c[3];
    double pos[2][3];
    double vel[2][3];
    for(int d=0;d<3;d++) {
        comx[d] = (mb * xb[d])/mb;
        comv[d] = (mb * vxb[d])/mb;
    }
    
    for(int d=0; d<3; d++) {
        xb[d] -= comx[d];
        xc[d] -= comx[d];
        vxb[d] -= comv[d];
        vxc[d] -= comv[d];
    }
    
    for(int d=0;d<3;d++) {
        pos[0][d] = xb[d];
        pos[1][d] = xc[d];
        
        vel[0][d] = vxb[d]/VSCALE;
        vel[1][d] = vxc[d]/VSCALE;
    }
    
    //evolve_kepler_2body(pos,vel,t_end*VSCALE);
    double t=0;
    double dt = 0.0;
    initial(pos,vel,ac,jc,a2c,a3c,&dt);
    evolve_system(pos,vel,ac,jc,a2c,a3c,&dt,&t,(t+t_end*VSCALE));
    for(int d=0;d<3;d++) {
        xb[d] = pos[0][d];
        vxb[d] = vel[0][d]*VSCALE;
        xci[d] = pos[1][d];
        vxci[d] = vel[1][d]*VSCALE;
    }

}

