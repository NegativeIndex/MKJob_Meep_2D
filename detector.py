import meep as mp
import argparse
import math
import cmath
import sys,os
import numpy as np
import datetime

folder=os.path.join(os.path.expanduser('~'),'python-study')
sys.path.append(folder)
import MeepFunctions.my_helper as my
import MeepFunctions.my_meep as mymeep
import meep.materials as mat

#########################################
#### common parameters
#########################################
class common:
    pols=["Ez",]
    geoms=["Empty","RC"]
    fcens=np.arange(0.245,0.31,0.03)
    dfs=np.ones(fcens.size)*0.03

    pol="Ez"
    geom="RC"
    fcen=0.245
    df=0.03

    sig='{}_{}_F{:5.3f}'.format(geom,pol,fcen)


#########################################
#### material
#########################################
class mymat:
    n_GaAs=3.72   
    n_AlAsSb=3.16 
    n_air=1
  
    # materials
    # define InAsSb
    wl0=3.72
    f0=1/wl0
    alpha=4000
    nreal=3.6
    nimag=alpha*wl0/(4*math.pi)*1e-4
    epsi=(nreal+1j*nimag)**2
    InAsSb=mp.Medium(epsilon=epsi.real,
                       D_conductivity=2*math.pi*f0*epsi.imag/epsi.real)    

 
#########################################
#### simulation function
#########################################
def simulation_fun():
    resolution=200
    geom=common.geom
    sig=common.sig
    
    # material 
    mat_src=mp.Medium(index=mymat.n_air)
    mat_low=mp.Medium(index=mymat.n_AlAsSb)
    mat_high=mp.Medium(index=mymat.n_GaAs)
    mat_active=mymat.InAsSb
  
    t_low=0.285
    t_high=0.236
    nDBR1=12  # bottom layer is low, top layer is high
    nDBR2=5  
    t_a1=0.096  # active layer
    t_a2=0.03
    t_b1=0.211  # barrier layer, low
    t_b2=0.211

    t_sub=1  # GaSb substrate
    t_air=3  # top air

    # build the model
    dpml=2.0            
    sx=0.5
    sy_DBR=(nDBR1+nDBR2)*(t_low+t_high)+t_high  # DBR
    sy_other=t_a1+t_a2+t_b1+t_b2                  # active and barrier
    sy=sy_DBR+sy_other+t_sub+t_air
    sxx=sx
    syy=sy+2*dpml

    cell_size = mp.Vector3(sxx,syy,0)
    pml_layers = [mp.PML(thickness=dpml,direction=mp.Y)]

    fcen=common.fcen
    df=common.df
    df2=df*1.4             # source frequency width
    nfreq = 201            # number of frequency bins

    # geometry
    geometry=[mp.Block(size=mp.Vector3(sxx,syy,mp.inf),
                       center=mp.Vector3(0,0,0),
                       material=mat_src)]

    if geom != "Empty":
        # all the material region
        y0=-syy/2                       # bottom y
        t0=dpml+t_sub+sy_DBR+sy_other   # thickness
        geometry.append(mp.Block(size=mp.Vector3(sxx,t0,mp.inf),
                                 center=mp.Vector3(0,y0+t0/2,0),
                                 material=mat_high))

        # bottom DBR
        y0=-sy/2+t_sub  # bottom y
        yc_low=[y0+i*(t_low+t_high)+t_low/2 for i in range(nDBR1)]
        geom_DBR1=[mp.Block(size=mp.Vector3(sxx,t_low,mp.inf),
                            center=mp.Vector3(0,yc,0),
                            material=mat_low) for yc in yc_low]
        geometry+=geom_DBR1

        # barrier, active 
        y0=-sy/2+t_sub+nDBR1*(t_high+t_low)
        t0=t_b1
        geometry.append(mp.Block(size=mp.Vector3(sxx,t0,mp.inf),
                                 center=mp.Vector3(0,y0+t0/2,0),
                                 material=mat_low))
        y0+=t0
        t0=t_a1
        geometry.append(mp.Block(size=mp.Vector3(sxx,t0,mp.inf),
                                 center=mp.Vector3(0,y0+t0/2,0),
                                 material=mat_active))

        y0+=t0
        t0=t_b2
        geometry.append(mp.Block(size=mp.Vector3(sxx,t0,mp.inf),
                                 center=mp.Vector3(0,y0+t0/2,0),
                                 material=mat_low))
        y0+=t0
        t0=t_a2
        geometry.append(mp.Block(size=mp.Vector3(sxx,t0,mp.inf),
                                 center=mp.Vector3(0,y0+t0/2,0),
                                 material=mat_active))
        # top DBR
        y0+=t0+t_high
        yc_low=[y0+i*(t_low+t_high)+t_low/2 for i in range(nDBR2)]
        geom_DBR2=[mp.Block(size=mp.Vector3(sxx,t_low,mp.inf),
                            center=mp.Vector3(0,yc,0),
                            material=mat_low) for yc in yc_low]
        geometry+=geom_DBR2
  
    # rotation angle (in degrees) of source
    theta=0
    theta_r=math.radians(theta)

    # bloch k vector
    k=mp.Vector3(math.sin(theta_r),math.cos(theta_r),0).scale(
        fcen*mymeep.get_refractive_index(fcen,mat_src).real)
    amp_func=lambda x: cmath.exp(1j*2*math.pi*k.dot(x))

    pol=common.pol
    src_cmpt=eval('mp.{}'.format(pol))

    # oblique source
    y0=sy/2-t_air/4
    sources = [mp.Source(mp.GaussianSource(fcen,fwidth=df2),
                         component=src_cmpt,
                         center=mp.Vector3(0,y0),
                         size=mp.Vector3(sx,0,0),
                         amp_func=amp_func
                     )]

    
    # setup simulations
    sim = mp.Simulation(cell_size=cell_size,
                        geometry=geometry,
                        boundary_layers=pml_layers,
                        sources=sources,
                        k_point=k,
                        filename_prefix=sig,
                        resolution=resolution)

    # flux 
    y_frs=[-0.5*sy+t_sub/2, 0.5*sy-t_air/4*3]

    frs=[mp.FluxRegion(center=mp.Vector3(0,y,0),
                       size=mp.Vector3(sx,0,0),
                       weight=-1) for y in y_frs]
    trans=[sim.add_flux(fcen,df,nfreq,fr) for fr in frs]

    # define step function to display fluxes
    my_display_count=0
    step_time=40
    def my_display_fluxes(sim):
        nonlocal my_display_count
        nonlocal step_time
        freqs=mp.get_flux_freqs(trans[0])
        fluxes=[mp.get_fluxes(tran) for tran in trans]
        data=np.column_stack((freqs, *fluxes,))
        my.matrix_output(None,data,"{:10.3e}","flux")
        my_display_count+=1
        print('No. {} display at t={}'.format(
            my_display_count,step_time))
        mymeep.my_flush_step(sim)


    # run simulations
    pt_field=mp.Vector3(0.1*sx,-0.05*sy,0) # monitor point
    sim.run(mp.at_beginning(mp.output_epsilon),
            mp.at_every(10, mymeep.my_flush_step),
            mp.at_every(step_time, my_display_fluxes),
            until_after_sources=mp.stop_when_fields_decayed(
                step_time, src_cmpt, pt_field, 1e-9))
    sys.stdout.flush()

    # collect data
    freqs = mp.get_flux_freqs(trans[0])
    fluxes=[mp.get_fluxes(tran) for tran in trans]
    data=np.column_stack((freqs,*fluxes,))


    # output
    fname=sig+".dat"
    my.matrix_output(fname,data,"{:10.3e}","flux")

######################### 
# main function 
##########################
if __name__ == "__main__":
    time0=datetime.datetime.now()
    print('-'*50)
    print(common.sig)
    simulation_fun()

    time1=datetime.datetime.now()
    dtime=time1-time0
    print("Simulation used {:0.2f} seconds".format(
        dtime.total_seconds()))
    sys.stdout.flush()
