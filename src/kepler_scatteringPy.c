#include "universal.h"
#define M_EARTH 3.0036060914862233e-6
#define M_MOON 3.695250762067745e-8
#define M_SUN 1.0
#define R_SUN 0.00465047
#define R_EARTH 0.0000426342966684
#define R_MOON 0.0000116184808779
#define VSCALE 6.2830666414874994
#define ETA_SYM 0.0125
#define ETA 0.001
#define max(a,b) \
    ({ __typeof__ (a) _a = (a); \
     __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

#define min(a,b) \
    ({ __typeof__ (a) _a = (a); \
     __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })

void debug_print(double pos[4][3], double vel[4][3]) {
    /*debug */
    printf("%.16e %.16e %.16e %.16e %.16e %.16e ",pos[0][0],pos[0][1],pos[0][2],vel[0][0],vel[0][1],vel[0][2]);
    printf("%.16e %.16e %.16e %.16e %.16e %.16e ",pos[1][0],pos[1][1],pos[1][2],vel[1][0],vel[1][1],vel[1][2]);
    printf("%.16e %.16e %.16e %.16e %.16e %.16e ",pos[2][0],pos[2][1],pos[2][2],vel[2][0],vel[2][1],vel[2][2]);
    printf("%.16e %.16e %.16e %.16e %.16e %.16e \n",pos[3][0],pos[3][1],pos[3][2],vel[3][0],vel[3][1],vel[3][2]);
    /*debug ends*/

}

void evolve_kepler_2body(double pos[2][3], double vel[2][3], double dt) {
    double xb[3] = {pos[0][0], pos[0][1], pos[0][2]}; 
    double xc[3] = {pos[1][0], pos[1][1], pos[1][2]}; 
    
    double vxb[3] = {vel[0][0], vel[0][1], vel[0][2]}; 
    double vxc[3] = {vel[1][0], vel[1][1], vel[1][2]}; 
    
    double dpos[3],dpos0[3];
    double pos_cm[3];
    double dvel[3],dvel0[3];
    double vel_cm[3];
    double m1 = M_SUN+M_EARTH+M_MOON;
    double m2 = 0.0;
    double mtot = (m1+m2);
    double f1 = m2 / mtot;
    double f2 = m1 / mtot;

    for(int d=0; d<3; d++){
        dpos0[d] = xb[d]-xc[d];
        dvel0[d] = vxb[d]-vxc[d];
    }

    for(int d=0; d<3; d++){
        pos_cm[d] = (m1*xb[d] + m2*xc[d]) / mtot;
        vel_cm[d] = (m1*vxb[d] + m2*vxc[d]) / mtot;
    }
    // evolve center of mass for dt
    for(int d=0; d<3; d++)
        pos_cm[d] += vel_cm[d] * dt;

    // call the Kepler solver
    Binary state0 = {dpos0[0],dpos0[1],dpos0[2], dvel0[0],dvel0[1],dvel0[2], 0, 0, 0};
    Binary state1;

    kepler_step(mtot, dt, &state0, &state1);

    // translate coordinates from 2-body frame to original frame
    xb[0] = pos_cm[0]+f1*state1.x;
    xb[1] = pos_cm[1]+f1*state1.y;
    xb[2] = pos_cm[2]+f1*state1.z;
    vxb[0] = vel_cm[0]+f1*state1.xd;
    vxb[1] = vel_cm[1]+f1*state1.yd;
    vxb[2] = vel_cm[2]+f1*state1.zd;
    xc[0] = pos_cm[0]-f2*state1.x;
    xc[1] = pos_cm[1]-f2*state1.y;
    xc[2] = pos_cm[2]-f2*state1.z;
    vxc[0] = vel_cm[0]-f2*state1.xd;
    vxc[1] = vel_cm[1]-f2*state1.yd;
    vxc[2] = vel_cm[2]-f2*state1.zd;
    
    for(int d=0;d<3;d++) {
        pos[0][d] = xb[d];
        pos[1][d] = xc[d];
        vel[0][d]= vxb[d];
        vel[1][d]= vxc[d];
    }


}

void evolve_kepler(double pos[4][3], double vel[4][3], double dt) {
    /*scale velocities to G=1 units first*/
    /*in (Msun,AU,yr) to obtain G=1 units, divide velocity by 2*pi*/
    double xs[3] = {pos[0][0], pos[0][1], pos[0][2]}; 
    double xe[3] = {pos[1][0], pos[1][1], pos[1][2]}; 
    double xm[3] = {pos[2][0], pos[2][1], pos[2][2]}; 

    double vxs[3] = {vel[0][0], vel[0][1], vel[0][2]}; 
    double vxe[3] = {vel[1][0], vel[1][1], vel[1][2]}; 
    double vxm[3] = {vel[2][0], vel[2][1], vel[2][2]}; 


    double comx[3], comv[3];
    for(int d=0;d<3;d++) {
        comx[d] = (M_SUN*xs[d] + M_EARTH*xe[d] +M_MOON*xm[d])/(M_SUN+M_EARTH+M_MOON);
        comv[d] = (M_SUN*vxs[d] + M_EARTH*vxe[d] +M_MOON*vxm[d])/(M_SUN+M_EARTH+M_MOON);
    }

    for(int d=0; d<3; d++) {
        xs[d] -= comx[d];
        xe[d] -= comx[d];
        xm[d] -= comx[d];

        vxs[d] -= comv[d];
        vxe[d] -= comv[d];
        vxm[d] -= comv[d];
    }

    // sun-(earth-moon) system first
   double em_pos[3], em_vel[3];

    for(int d=0;d<3;d++){
        em_pos[d] = (M_EARTH*xe[d]+M_MOON*xm[d])/(M_EARTH+M_MOON);
        em_vel[d] = (M_EARTH*vxe[d]+M_MOON*vxm[d])/(M_EARTH+M_MOON);
    }

    double dpos[3],dpos0[3];
    double pos_cm[3];
    double dvel[3],dvel0[3];
    double vel_cm[3];
    double m1 = M_SUN;
    double m2 = (M_MOON+M_EARTH);
    double mtot = (m1+m2);
    double f1 = m2 / mtot;
    double f2 = m1 / mtot;

    for(int d=0; d<3; d++){
        dpos0[d] = xs[d]-em_pos[d];
        dvel0[d] = vxs[d]-em_vel[d];
    }

    for(int d=0; d<3; d++){
        pos_cm[d] = (m1*xs[d] + m2*em_pos[d]) / mtot;
        vel_cm[d] = (m1*vxs[d] + m2*em_vel[d]) / mtot;
    }
    // evolve center of mass for dt
    //for(int d=0; d<3; d++)
    //    pos_cm[d] += vel_cm[d] * dt;

    // call the Kepler solver
    Binary state0 = {dpos0[0],dpos0[1],dpos0[2], dvel0[0],dvel0[1],dvel0[2], 0, 0, 0};
    Binary state1;

    kepler_step(mtot, dt, &state0, &state1);

    // translate coordinates from 2-body frame to original frame
    xs[0] = pos_cm[0]+f1*state1.x;
    xs[1] = pos_cm[1]+f1*state1.y;
    xs[2] = pos_cm[2]+f1*state1.z;
    vxs[0] = vel_cm[0]+f1*state1.xd;
    vxs[1] = vel_cm[1]+f1*state1.yd;
    vxs[2] = vel_cm[2]+f1*state1.zd;
    em_pos[0] = pos_cm[0]-f2*state1.x;
    em_pos[1] = pos_cm[1]-f2*state1.y;
    em_pos[2] = pos_cm[2]-f2*state1.z;
    em_vel[0] = vel_cm[0]-f2*state1.xd;
    em_vel[1] = vel_cm[1]-f2*state1.yd;
    em_vel[2] = vel_cm[2]-f2*state1.zd;

    //earth-moon system next
    double dpos0_em[3], dvel0_em[3];
    for(int d=0; d<3; d++){
        dpos0_em[d] = xe[d]-xm[d];
        dvel0_em[d] = vxe[d]-vxm[d];
    }
    double m1em, m2em;
    m1em = M_EARTH;
    m2em = M_MOON;
    double mtotem = (m1em+m2em);
    double f1em = m2em / mtotem;
    double f2em = m1em / mtotem;
    //for(int d=0; d<3; d++)
    //    em_pos[d] += em_vel[d] * dt * VSCALE;

    Binary state0em = {dpos0_em[0],dpos0_em[1],dpos0_em[2], dvel0_em[0],dvel0_em[1],dvel0_em[2], 0, 0, 0};
    Binary state1em;

    kepler_step(mtotem, dt, &state0em, &state1em);


    xe[0] = em_pos[0]+f1em*state1em.x;
    xe[1] = em_pos[1]+f1em*state1em.y;
    xe[2] = em_pos[2]+f1em*state1em.z;
    vxe[0] = em_vel[0]+f1em*state1em.xd;
    vxe[1] = em_vel[1]+f1em*state1em.yd;
    vxe[2] = em_vel[2]+f1em*state1em.zd;
    xm[0] = em_pos[0]-f2em*state1em.x;
    xm[1] = em_pos[1]-f2em*state1em.y;
    xm[2] = em_pos[2]-f2em*state1em.z;
    vxm[0] = em_vel[0]-f2em*state1em.xd;
    vxm[1] = em_vel[1]-f2em*state1em.yd;
    vxm[2] = em_vel[2]-f2em*state1em.zd;


    /*for(int d=0; d<3; d++) {
      comx[d] += dt * comv[d];
      }*/

    for(int d=0; d<3; d++) {
        xs[d] += comx[d];
        xe[d] += comx[d];
        xm[d] += comx[d];

        vxs[d] += comv[d];
        vxe[d] += comv[d];
        vxm[d] += comv[d];
    }


    /*rescale back to physical units*/

    for(int d=0;d<3;d++) {
        pos[0][d] = xs[d];
        pos[1][d] = xe[d];
        pos[2][d] = xm[d];
        vel[0][d]= vxs[d];
        vel[1][d]= vxe[d];
        vel[2][d]= vxm[d];
    }


}



void force_and_jerk(double pos[4][3], double vel[4][3], double *tpos, double *tvel,
        double *ac, double *jc) {

    double accx = 0.0;
    double accy = 0.0;
    double accz = 0.0;

    double jerkx = 0.0;
    double jerky = 0.0;
    double jerkz = 0.0;


    double mass[4] = {M_SUN, M_EARTH, M_MOON, 0.0};
    for(int j=0; j<3;j++) {
        double dx = tpos[0] - pos[j][0];
        double dy = tpos[1] - pos[j][1];
        double dz = tpos[2] - pos[j][2];
        double dist = sqrt(dx*dx + dy*dy + dz*dz );

        double dist2 = dist*dist;
        double dist3 = dist*dist*dist;

        double mass_over_dist3 = mass[j]/dist3;
        accx += -mass_over_dist3*dx;
        accy += -mass_over_dist3*dy;
        accz += -mass_over_dist3*dz;

        double dvx = tvel[0] - vel[j][0];
        double dvy = tvel[1] - vel[j][1];
        double dvz = tvel[2] - vel[j][2];
        double rrdot = (dx*dvx+dy*dvy+dz*dvz);


        jerkx += -mass[j]*(dvx/dist3 - 3*rrdot*dx/(dist3*dist2));
        jerky += -mass[j]*(dvy/dist3 - 3*rrdot*dy/(dist3*dist2));
        jerkz += -mass[j]*(dvz/dist3 - 3*rrdot*dz/(dist3*dist2));


    }
    ac[0] = accx;
    ac[1] = accy;
    ac[2] = accz;
    jc[0] = jerkx;
    jc[1] = jerky;
    jc[2] = jerkz;

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


void timestep(double pos[4][3], double vel[4][3], double *dt) {
    double stepi;
    stepi = 1.e38;
    double mass[4] = {M_SUN,M_EARTH,M_MOON,0.0};
    for (int j=0; j<3; j++) {
        double dx = pos[3][0] - pos[j][0];
        double dy = pos[3][1] - pos[j][1];
        double dz = pos[3][2] - pos[j][2];
        double r2 = dx*dx + dy*dy + dz*dz;
        double r = sqrt(r2);

        double dvx = vel[3][0] - vel[j][0];
        double dvy = vel[3][1] - vel[j][1];
        double dvz = vel[3][2] - vel[j][2];

        double v2 = dvx * dvx + dvy * dvy + dvz * dvz ; 


        double vdotdr2 = (dx*dvx+dy*dvy+dz*dvz) / r2;

        double tau = ETA_SYM / M_SQRT2 * r * sqrt(r/(mass[3]+mass[j]));
        double dtau = 3 * tau * vdotdr2 / 2;

        dtau = (dtau < 1) ? dtau : 1;
        tau /= (1 - dtau / 2);

        stepi = min(tau, stepi);
        tau = ETA_SYM * r / sqrt(v2);
        dtau = tau * vdotdr2 * (1 + (mass[3]+mass[j]) / (v2 * r));
        dtau = (dtau < 1)? dtau: 1;

        tau /= (1 - dtau / 2);
        stepi = min(tau, stepi);

    }
    *dt = stepi;
}

void initial(double pos[4][3],double vel[4][3], double *ac, double *jc, double *dt)  {
    double tpos[3] = {pos[3][0],pos[3][1],pos[3][2]};
    double tvel[3] = {vel[3][0],vel[3][1],vel[3][2]};
    force_and_jerk(pos,vel,tpos,tvel,ac,jc);
    timestep(pos,vel,dt);
}

void predictor(double pos[4][3], double vel[4][3], double *ac, double *jc, double dt,
        double *tpos, double *tvel) {
    
    evolve_kepler(pos,vel,dt);

    double dt2_over_2 = dt*dt/2;
    double dt3_over_6 = dt*dt*dt/6;
    
    for(int d=0;d<3;d++) {
        tpos[d] = pos[3][d] + vel[3][d]*dt + ac[d]*dt2_over_2 + jc[d]*dt3_over_6;
        tvel[d] = vel[3][d] + ac[d]*dt + jc[d]*dt2_over_2;
    }
}

void corrector (double pos[4][3], double vel[4][3], double *ac, double *jc,
        double *tpos, double *tvel, double *tposc, double *tvelc,double *aci, double *jci, double *a2ci, double *a3ci,double dt) {
    double a2i[3], a3i[3];
    double step_i2 = dt * dt;
    double step_i3 =dt * dt * dt;

    for(int d=0;d<3;d++){
        a2i[d] = 2*(-3*(ac[d]-aci[d]) - dt*(2*jc[d] + 1*jci[d]))/(step_i2);
        a3i[d] = 6*(2*(ac[d]-aci[d]) + dt*(jc[d] + jci[d]))/(step_i3);
    }

    for(int d=0;d<3;d++) {
        tposc[d] = tpos[d] + (step_i2*step_i2)*(a2i[d])/24 + (step_i2*step_i3)*(a3i[d])/120;
        tvelc[d] = tvel[d] + (step_i3)*(a2i[d])/6 + (step_i2*step_i2)*(a3i[d])/24;
    }
    for(int d=0;d<3;d++) {
        a2ci[d] = a2i[d];
        a3ci[d] = a3i[d];
    }


}

void evolve_hermite(double pos[4][3], double vel[4][3], double *ac, double *jc, double *dt) {
    double tpos[3], tvel[3];
    double tposc[3], tvelc[3];
    double aci[3], jci[3];
    double a2ci[3], a3ci[3];

    predictor(pos,vel,ac,jc,*dt,tpos,tvel);

    force_and_jerk(pos,vel,tpos,tvel,aci,jci);
    corrector(pos,vel,ac,jc,tpos,tvel,tposc,tvelc,aci,jci,a2ci,a3ci,*dt);
    
    force_and_jerk(pos,vel,tposc,tvelc,aci,jci);
    corrector(pos,vel,ac,jc,tpos,tvel,tposc,tvelc,aci,jci,a2ci,a3ci,*dt);
    
    force_and_jerk(pos,vel,tposc,tvelc,aci,jci);
    corrector(pos,vel,ac,jc,tpos,tvel,tposc,tvelc,aci,jci,a2ci,a3ci,*dt);

    //force_and_jerk(pos,vel,tposc,tvelc,aci,jci);
    //corrector(pos,vel,ac,jc,tpos,tvel,tposc,tvelc,aci,jci,dt);

    /*for(int i=0;i<2;i++) {
        force_and_jerk(pos,vel,tpos,tvel,aci,jci);
        corrector(pos,vel,ac,jc,tpos,tvel,aci,jci,dt);
    }*/
    for(int d=0;d<3;d++) {
        pos[3][d] = tposc[d];
        vel[3][d] = tvelc[d];
        ac[d] = aci[d];
        jc[d] = jci[d];
        //a2ci[d] = a2ci[d] + *dt * a3ci[d];
    }
    //timestep_aarseth(aci,jci,a2ci,a3ci,dt);

}

int check_collision_and_bound(double pos[4][3],double vel[4][3]) {
    double x = pos[3][0];
    double y = pos[3][1];
    double z = pos[3][2];

    double v2 = vel[3][0] * vel[3][0] + vel[3][1] * vel[3][1] + vel[3][2] * vel[3][2];
    
    double ene = 0.5 * v2;

    int flag = 0; //0 = no collision, 1 = collision with sun, 2 = collision with earth, 3 = collision with moon, 4 = bound
    double mass[3] = {M_SUN,M_EARTH,M_MOON}; 
    double thresh[3] = {R_SUN, R_EARTH, R_MOON};
    int colcode[3] = {1,2,3};
    for(int j=0;j<3;j++) {
        double dx = pos[j][0] - x;
        double dy = pos[j][1] - y;
        double dz = pos[j][2] - z;
        double d = sqrt(dx*dx + dy*dy + dz*dz);
        ene -= (mass[j] / d);
        if(d <= thresh[j]) {
            flag = colcode[j];
            break;
        }
    }
    if(flag==0 && ene < 0) {
        flag = 4;
    }
    return flag;
}

int evolve_system(double pos[4][3], double vel[4][3], double *ac, double *jc, double *dt, double *t, double t_end) {
    double t_cur = *t;
    int flag = 0;
    //printf("%.7e %.7e\n",t_cur, *dt);
    while(t_cur < t_end) {
        if((t_cur+*dt) > t_end)
            *dt  = t_end-t_cur;
        t_cur += *dt;

        evolve_hermite(pos,vel,ac,jc,dt);
        //debug_print(pos,vel);
        flag = check_collision_and_bound(pos,vel);
        if(flag > 0 && flag < 4) {
            break;
        }
        timestep(pos,vel,dt);
    }
    return flag;
}

double scale2nbody(double vel[4][3]) {
    for(size_t d=0; d<3;d++) {
        vel[0][d] = vel[0][d] / VSCALE;
        vel[1][d] = vel[1][d] / VSCALE;
        vel[2][d] = vel[2][d] / VSCALE;
        vel[3][d] = vel[3][d] / VSCALE;
    }
}
double scale2physical(double vel[4][3]) {
    for(size_t d=0; d<3;d++) {
        vel[0][d] = vel[0][d] * VSCALE;
        vel[1][d] = vel[1][d] * VSCALE;
        vel[2][d] = vel[2][d] * VSCALE;
        vel[3][d] = vel[3][d] * VSCALE;
    }
}


void integrate(double *xs, double *xe, double *xm, double *xc, 
        double *vxs, double *vxe, double *vxm, double *vxc,
        double t_end, int* collb_flag) {
   
    double comx[3], comv[3];

    double pos[4][3];
    double vel[4][3];
    double ac[3], jc[3];
    for(int d=0;d<3;d++) {
        comx[d] = (M_SUN*xs[d] + M_EARTH*xe[d] +M_MOON*xm[d])/(M_SUN+M_EARTH+M_MOON);
        comv[d] = (M_SUN*vxs[d] + M_EARTH*vxe[d] +M_MOON*vxm[d])/(M_SUN+M_EARTH+M_MOON);
    }
    
    for(int d=0; d<3; d++) {
        xs[d] -= comx[d];
        xe[d] -= comx[d];
        xm[d] -= comx[d];
        xc[d] -= comx[d];
        vxs[d] -= comv[d];
        vxe[d] -= comv[d];
        vxm[d] -= comv[d];
        vxc[d] -= comv[d];
    }
    


    for(int d=0;d<3;d++) {
        pos[0][d] = xs[d];
        pos[1][d] = xe[d];
        pos[2][d] = xm[d];
        pos[3][d] = xc[d];
        
        vel[0][d] = vxs[d];
        vel[1][d] = vxe[d];
        vel[2][d] = vxm[d];
        vel[3][d] = vxc[d];
    }
    double t=0;
    //double t_end = 5.0;
    //double odt = 0.01;
    double dt = 0.0;
    //initial(pos,vel,ac,jc,&dt);

    scale2nbody(vel);

    initial(pos,vel,ac,jc,&dt);
    int flag = evolve_system(pos,vel,ac,jc,&dt,&t,(t+t_end*VSCALE));
    flag = check_collision_and_bound(pos,vel);
    scale2physical(vel);
    *collb_flag = flag;

    for(int d=0;d<3;d++) {
        xs[d] = pos[0][d];
        xe[d] = pos[1][d];
        xm[d] = pos[2][d];
        xc[d] = pos[3][d];
        
        vxs[d] = vel[0][d];
        vxe[d] = vel[1][d];
        vxm[d] = vel[2][d];
        vxc[d] = vel[3][d];
    }

}

void integrate2body(double *xb, double *xc, double *vxb, double *vxc, double t_end, double *xci, double *vxci) {
    double comx[3], comv[3];
    double mb = M_SUN + M_EARTH + M_MOON;
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
    evolve_kepler_2body(pos,vel,t_end*VSCALE);
    
    for(int d=0;d<3;d++) {
        xb[d] = pos[0][d];
        vxb[d] = vel[0][d]*VSCALE;
        xci[d] = pos[1][d];
        vxci[d] = vel[1][d]*VSCALE;
    }

}
/*int main() {
    double xs[3] = {0.0,0.0,0.0};
    double xe[3] = {0.9833000000000001, 0.0, 0.0};
    double xm[3] = {0.9857306852303347, 0.0, 0.0};
    double xc[3] = {0.999, 0.0, 0.0};
    //double xc[3] = {0.001,0.0,0.0};
    
    //double xc[3] = {0.64, 0.0, 0.0};
    //double xc[3] = {0.040000000000000015, 0.0, 0.0}; 


    double vxs[3] = {0.0, 0.0, 0.0};
    double vxe[3] = {0.0, 6.388894409407078, 0.0};
    double vxm[3] = {0.0, 6.61713329602396, 0.0};
    //double vxc[3] = {0.0, 280.28393385731573, 0.0};
    //double vxc[3] = {0.0, 8.603443324627317, 0.0};
    double vxc[3] = {0.0, 6.286210532923412, 0.0};
    //double vxc[3] = {0.0, 43.869125396364275, 0.0};
    double comx[3], comv[3];

    double pos[4][3];
    double vel[4][3];
    double ac[3], jc[3];
    for(int d=0;d<3;d++) {
        comx[d] = (M_SUN*xs[d] + M_EARTH*xe[d] +M_MOON*xm[d])/(M_SUN+M_EARTH+M_MOON);
        comv[d] = (M_SUN*vxs[d] + M_EARTH*vxe[d] +M_MOON*vxm[d])/(M_SUN+M_EARTH+M_MOON);
    }
    
    for(int d=0; d<3; d++) {
        xs[d] -= comx[d];
        xe[d] -= comx[d];
        xm[d] -= comx[d];
        xc[d] -= comx[d];
        vxs[d] -= comv[d];
        vxe[d] -= comv[d];
        vxm[d] -= comv[d];
        vxc[d] -= comv[d];
    }
    


    for(int d=0;d<3;d++) {
        pos[0][d] = xs[d];
        pos[1][d] = xe[d];
        pos[2][d] = xm[d];
        pos[3][d] = xc[d];
        
        vel[0][d] = vxs[d];
        vel[1][d] = vxe[d];
        vel[2][d] = vxm[d];
        vel[3][d] = vxc[d];
    }
    double t=0;
    double t_end = 5.0;
    double odt = 0.01;
    double dt = 0.0;
    initial(pos,vel,ac,jc,&dt);
    while (t < t_end) {
        scale2nbody(vel);
        evolve_system(pos,vel,ac,jc,&dt,&t,(t+odt*VSCALE)); 
        t+=odt;
        scale2physical(vel);
        debug_print(pos,vel);
    }


}*/
