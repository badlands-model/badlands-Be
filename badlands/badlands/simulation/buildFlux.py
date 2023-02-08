##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
##                                                                                   ##
##  This file forms part of the Badlands surface processes modelling application.    ##
##                                                                                   ##
##  For full license and copyright information, please refer to the LICENSE.md file  ##
##  located at the project root, or contact the authors.                             ##
##                                                                                   ##
##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
"""
This file is the main entry point to compute flow network and associated sedimentary fluxes.
"""

import sys
import time
import numpy as np
from matplotlib import path

import os
if 'READTHEDOCS' not in os.environ:
    from badlands import (elevationTIN,buildMesh)

def streamflow(input, FVmesh, recGrid, force, hillslope, flow, elevation, \
                 lGIDs, rain, tNow, melt=None, ice=None, verbose=False):
    """
    Compute stream flow.

    Args:
        input: class containing XML input file parameters.
        FVmesh: class describing the finite volume mesh.
        recGrid: class describing the regular grid characteristics.
        force: class describing the forcing parameters.
        hillslope: class describing hillslope processes.
        flow: class describing stream power law processes.
        elevation: numpy array containing the elevations for the domain.
        lGIDs: numpy 1D array containing the node indices.
        rain: numpy 1D array containing rainfall precipitation values.
        tNow: simulation time step.
        melt : numpy array containing meltwaters height (to be added to discharge)
        ice: numpy array containing ice thickness
        verbose : (bool) when :code:`True`, output additional debug information (default: :code:`False`).

    Returns
    -------
    fillH
        numpy 1D array containing the depression-less elevation values.
    elev
        numpy 1D array containing the elevations.
    """

    # Update sea-level
    walltime = time.time()
    ref_elev = buildMesh.get_reference_elevation(input, recGrid, elevation)
    force.getSea(tNow,input.udw,ref_elev)
    fillH = None

    # Update river input
    force.getRivers(tNow)
    riverrain = rain+force.rivQw

    # Build an initial depression-less surface at start time if required
    if input.tStart == tNow and input.nopit == 1 :
        fillH = elevationTIN.pit_stack(elevation,input.nopit,force.sealevel)
        elevation = fillH
    else:
        fillH = elevationTIN.pit_stack(elevation,0,force.sealevel)

    if verbose and input.spl:
        print(" -   depression-less algorithm PD with stack", time.time() - walltime)

    # Compute stream network
    walltime = time.time()
    flow.SFD_receivers(fillH, elevation, FVmesh.neighbours,
                       FVmesh.vor_edges, FVmesh.edge_length,
                       lGIDs)

    if verbose:
        print(" -   compute receivers parallel ", time.time() - walltime)

    # Distribute evenly local minimas to processors on filled surface
    walltime = time.time()
    flow.localbase = flow.base
    flow.ordered_node_array_filled()
    if verbose:
        print(" -   compute stack order locally for filled surface", time.time() - walltime)

    walltime = time.time()
    flow.stack = flow.localstack
    if verbose:
        print(" -   send stack order for filled surface globally ", time.time() - walltime)

    # Distribute evenly local minimas on real surface
    walltime = time.time()
    flow.localbase1 = flow.base1
    flow.ordered_node_array_elev()
    if verbose:
        print(" -   compute stack order locally for real surface", time.time() - walltime)

    walltime = time.time()
    flow.stack1 = flow.localstack1
    if verbose:
        print(" -   send stack order for real surface globally ", time.time() - walltime)

    # Compute a unique ID for each local depression and their downstream draining nodes
    flow.compute_parameters_depression(fillH,elevation,FVmesh.control_volumes,force.sealevel)

    # Compute discharge
    walltime = time.time()
    flow.compute_flow(force.sealevel, elevation, FVmesh.control_volumes, riverrain, melt, ice)
    if verbose:
        print(" -   compute discharge ", time.time() - walltime)

    return fillH, elevation

def sediment_flux(input, recGrid, hillslope, FVmesh, tMesh, flow, force, rain, lGIDs, applyDisp, straTIN, \
                  mapero, cumdiff, cumhill, cumfail, bedrock_Be, qprop, ice, fillH, disp, inGIDs, elevation, \
                  tNow, tEnd, verbose=False):
    """
    Compute sediment fluxes.

    Args:
        input: class containing XML input file parameters.
        recGrid: class describing the regular grid characteristics.
        hillslope: class describing hillslope processes.
        FVmesh: finite volume mesh.
        tMesh: TIN mesh.
        flow: flow parameters.
        force:  forcing parameters.
        rain: rain values.
        lGIDs: global nodes indices.
        applyDisp: applying displacement boolean.
        straTIN: class for stratigraphic TIN grid.
        mapero: imposed erodibility maps.
        cumdiff: cumulative total erosion/deposition changes.
        cumhill: cumulative hillslope erosion/deposition changes.
        cumfail: cumulative failure induced erosion/deposition changes.
        bedrock_Be: 10Be concentration (at of Qtz) in bedrock or old sediments.
        qprop : quartz proportion in bedrock or old sediments (1 for 100% quartz).
        ice : ice thickness.
        fillH: filled elevation mesh.
        disp: displacement values .
        inGIDs: nodes indices.
        elevation: elevation mesh
        tNow: simulation time step.
        tEnd: simulation end time.
        verbose : (bool) when :code:`True`, output additional debug information (default: :code:`False`).

    Returns
    -------
    tNow
        simulation current time.
    elevation
        updated numpy 1D array containing the elevations.
    cumdiff
        updated cumulative total erosion/deposition changes
    cumhill
        updated cumulative hillslope erosion/deposition changes
    cumfail
        updated cumulative failure induced erosion/deposition changes
    bedrock_Be
        updated 10Be concentration
    qprop
        updated Quartz proportion

    """

    flow_time = time.time()

    # Get active layer
    if straTIN is not None:
        walltime = time.time()
        flow.activelay[flow.activelay<1.] = 1.
        flow.activelay[flow.activelay>straTIN.activeh] = straTIN.activeh
        straTIN.get_active_layer(flow.activelay,verbose)
        activelay = straTIN.alayR
        flow.straTIN = 1
        # Set the average erodibility based on rock types in the active layer
        flow.erodibility = np.sum(straTIN.rockCk*activelay/flow.activelay.reshape(len(elevation),1),axis=1)
        eroCk = straTIN.rockCk
        if verbose:
            print(" -   Get active layer ", time.time() - walltime)
    else:
        activelay = None
        eroCk = 0.

    # Find border/inside nodes
    if flow.domain is None:
        ids = np.arange(len(FVmesh.control_volumes))
        tmp1 = np.where(FVmesh.control_volumes>0.)[0]
        xyMin = [recGrid.regX.min()-1., recGrid.regY.min()-1.]
        xyMax = [recGrid.regX.max()+1., recGrid.regY.max()+1.]
        flow.domain = path.Path([(xyMin[0],xyMin[1]),(xyMax[0],xyMin[1]), (xyMax[0],xyMax[1]), (xyMin[0],xyMax[1])])
        tmp2 = flow.domain.contains_points(flow.xycoords)
        flow.insideIDs = np.intersect1d(tmp1,ids[tmp2])
        flow.borders = np.zeros(len(FVmesh.control_volumes),dtype=int)
        flow.borders[flow.insideIDs] = 1
        flow.outsideIDs = np.where(flow.borders==0)[0]
        xyMin2 = [recGrid.regX.min()+recGrid.resEdges, recGrid.regY.min()+recGrid.resEdges]
        xyMax2 = [recGrid.regX.max()-recGrid.resEdges, recGrid.regY.max()-recGrid.resEdges]
        xyMin2 = [recGrid.regX.min()+1, recGrid.regY.min()+1]
        xyMax2 = [recGrid.regX.max()-1, recGrid.regY.max()-1]
        domain = path.Path([(xyMin2[0],xyMin2[1]),(xyMax2[0],xyMin2[1]), (xyMax2[0],xyMax2[1]), (xyMin2[0],xyMax2[1])])
        tmp3 = domain.contains_points(flow.xycoords)
        flow.insideIDs2 = ids[tmp3]
        flow.borders2 = np.zeros(len(FVmesh.control_volumes),dtype=int)
        flow.borders2[flow.insideIDs2] = 1
        flow.outsideIDs2 = np.where(flow.borders2==0)[0]

    # Compute CFL condition
    walltime = time.time()
    if input.Hillslope and hillslope.updatedt == 0:
        if hillslope.Sc == 0:
            hillslope.dt_stability(FVmesh.edge_length[inGIDs,:tMesh.maxNgbh])
        else:
            hillslope.dt_stabilityCs(elevation, FVmesh.neighbours, FVmesh.edge_length,
                    lGIDs, flow.borders2)
            if hillslope.CFL < input.minDT:
                print('Decrease your hillslope diffusion coefficients to ensure stability.')
                sys.exit(0)
        hillslope.dt_stability_ms(FVmesh.edge_length[inGIDs,:tMesh.maxNgbh])
        hillslope.dt_stability_fail(FVmesh.edge_length[inGIDs,:tMesh.maxNgbh])
    elif hillslope.CFL is None:
        hillslope.CFL = tEnd-tNow

    flow.dt_stability(fillH, inGIDs)
    CFLtime = min(flow.CFL, hillslope.CFL)
    if CFLtime>1.:
        CFLtime = float(round(CFLtime-0.5,0))
    if verbose:
        print('CFL for hillslope and flow ',hillslope.CFL,flow.CFL,CFLtime)
    CFLtime = min(CFLtime, tEnd - tNow)
    if input.minDT > 1:
        if CFLtime < input.minDT:
            if input.minDT > tEnd - tNow:
                CFLtime = tEnd - tNow
            else:
                CFLtime = max(input.minDT, CFLtime)
        else:
            CFLtime = max(input.minDT, CFLtime)
    else:
        CFLtime = max(input.minDT, CFLtime)
    CFLtime = min(input.maxDT, CFLtime)
    if verbose:
        print(" -   Get CFL time step ", time.time() - walltime)

    # Compute sediment fluxes
    if (input.erolays and input.erolays >= 0):
        oldelev = np.copy(elevation)

    # Initial cumulative elevation change
    walltime = time.time()

    timestep, sedchange, erosion, deposition, sed_be, qprop, slopeTIN = flow.compute_sedflux(FVmesh.control_volumes, elevation, qprop, 
                                          bedrock_Be, rain, fillH, CFLtime, activelay, eroCk, force.rivQs, 
                                          force.sealevel, input.perc_dep, input.slp_cr, FVmesh.neighbours, verbose=False)

    if timestep < CFLtime:
        if input.minDT > tEnd - tNow:
            CFLtime = tEnd - tNow
        else:
            CFLtime = max(input.minDT, CFLtime)
    else:
        CFLtime = max(input.minDT, CFLtime)

    if verbose:
        print(" -   Get stream fluxes ", time.time() - walltime)

    ed = np.sum(sedchange,axis=1)
    
    elevation += ed
    cumdiff += ed

    # update 10Be concentration in newly deposited sediments
    sedBeIDs             = np.where(ed>0.)[0]
    bedrock_Be[sedBeIDs] = sed_be[sedBeIDs]

    # Compute marine sediment diffusion
    if hillslope.CDriver > 0.:
        walltime = time.time()

        # Initialise marine sediments diffusion array
        it = 0
        sumdep = np.sum(deposition,axis=1)
        maxth = 0.1
        diffstep = timestep
        diffcoeff = hillslope.sedfluxmarine(force.sealevel, elevation, FVmesh.control_volumes)

        # Perform river related sediment diffusion
        while diffstep > 0. and it < 1000:
            # Define maximum time step
            maxstep = min(hillslope.CFLms,diffstep)
            # Compute maximum marine fluxes and maximum timestep to avoid excessive diffusion erosion
            diffmarine, diffQz, diffBe, mindt = flow.compute_marine_diffusion(elevation, qprop, bedrock_Be, sumdep, FVmesh.neighbours, FVmesh.vor_edges,
                                            FVmesh.edge_length, diffcoeff, lGIDs, force.sealevel, maxth, maxstep)
            diffmarine[flow.outsideIDs] = 0.
            maxstep = min(mindt,maxstep)

            # Update Quartz and 10Be content in deposited sediments
            sedBeIDs             = np.where(np.logical_and(diffQz>0., diffmarine>0.))[0]
            qprop[sedBeIDs]      = diffQz[sedBeIDs]
            bedrock_Be[sedBeIDs] = diffBe[sedBeIDs]

            # if maxstep < input.minDT:
            #    print 'WARNING: marine diffusion time step is smaller than minimum timestep:',maxstep
            #    print 'You will need to decrease your diffusion coefficient for criver'
            #    stop

            # Update diffusion time step and total diffused thicknesses
            diffstep -= maxstep

            # Distribute rock based on their respective proportions in the deposited columns
            if straTIN is not None:
                # Compute multi-rock diffusion
                sedpropflux, difftot, diffQz, diffBe = flow.compute_sediment_marine(elevation, qprop, bedrock_Be, deposition, sumdep,
                                                diffcoeff*maxstep, FVmesh.neighbours, force.sealevel,
                                                maxth, FVmesh.vor_edges, FVmesh.edge_length, lGIDs)
                difftot[flow.outsideIDs] = 0.
                sedpropflux[flow.outsideIDs,:] = 0.

                # Update deposition for each rock type
                # Update Quartz and 10Be content in deposited sediments
                deposition += sedpropflux
                deposition[deposition<0] = 0.
                sedBeIDs             = np.where(np.logical_and(diffQz > 0., difftot > 0.))[0]
                qprop[sedBeIDs]      = diffQz[sedBeIDs]
                bedrock_Be[sedBeIDs] = diffBe[sedBeIDs]

                # Update elevation, erosion/deposition
                sumdep += difftot
                elevation += difftot
                cumdiff += difftot
            else:
                # Update elevation, erosion/deposition
                sumdep += diffmarine*maxstep
                elevation += diffmarine*maxstep
                cumdiff += diffmarine*maxstep

            it += 1

    # Compute slope failures
    if hillslope.Sfail > 0.:
        walltime = time.time()

        # Initialise sediments diffusion array
        it = 0
        walltime = time.time()
        erofail = flow.compute_failure(elevation, hillslope.Sfail)

        # Add slope failure erosion
        slumpID = np.where(erofail<0)[0]
        sumdep = -erofail
        maxth = 0.1
        diffstep = timestep
        diffcoeff = hillslope.sedfluxfailure(FVmesh.control_volumes)

        # Perform river related sediment diffusion
        if len(slumpID) > 0:
            while diffstep > 0. and it < 2000:
                # Define maximum time step
                maxstep = min(hillslope.CFLfail,diffstep)
                # Compute maximum marine fluxes and maximum timestep to avoid excessive diffusion erosion
                difffail, diffQz, diffBe, mindt = flow.compute_failure_diffusion(elevation, qprop, bedrock_Be, sumdep, FVmesh.neighbours, FVmesh.vor_edges,
                                                FVmesh.edge_length, diffcoeff, lGIDs, maxth, maxstep)

                difffail[flow.outsideIDs] = 0.
                maxstep = min(mindt,maxstep)

                # Update diffusion time step and total diffused thicknesses
                diffstep -= maxstep

                # Update elevation, erosion/deposition
                sumdep += difffail*maxstep
                elevation += difffail*maxstep
                cumdiff += difffail*maxstep
                cumfail += difffail*maxstep
                it += 1

                # Update Quartz and 10Be content in deposited sediments
                sedBeIDs             = np.where(np.logical_and(diffQz > 0, difffail > 0))[0]
                qprop[sedBeIDs]      = diffQz[sedBeIDs]
                bedrock_Be[sedBeIDs] = diffBe[sedBeIDs]

        if verbose:
            print(" -   Get slope failure sediment fluxes ", time.time() - walltime)

    # Compute hillslope processes
    dtype = 1
    if straTIN is None:
        dtype = 0
    walltime = time.time()
    area = np.copy(FVmesh.control_volumes)
    area[flow.outsideIDs2] = 0.
    diffcoeff = hillslope.sedflux(force.sealevel, ice, elevation, FVmesh.control_volumes)
    diffcoeff[flow.outsideIDs2] = 0.
    diff_flux, diffQz, diffBe = flow.compute_hillslope_diffusion(elevation, qprop, bedrock_Be, FVmesh.neighbours, FVmesh.vor_edges,
                       FVmesh.edge_length, lGIDs, dtype, hillslope.Sc)
    diff_flux[flow.outsideIDs2] = 0.
    cdiff = diffcoeff*diff_flux*timestep

    # Update Quartz and 10Be content in deposited sediments
    sedBeIDs             = np.where(np.logical_and(diffQz > 0., cdiff > 0.))[0]
    qprop[sedBeIDs]      = diffQz[sedBeIDs]
    bedrock_Be[sedBeIDs] = diffBe[sedBeIDs]

    if straTIN is None:
        if input.btype == 'outlet':
            cdiff[flow.insideIDs[0]] = 0.
        # Update dataset
        elevation[flow.insideIDs] += cdiff[flow.insideIDs]
        cumdiff[flow.insideIDs] += cdiff[flow.insideIDs]
        cumhill[flow.insideIDs] += cdiff[flow.insideIDs]
    else:
        straTIN.update_layers(erosion, deposition, elevation, verbose)

        # Get the active layer thickness to erode using diffusion
        maxlayh = -cdiff
        maxlayh[maxlayh<1.] = 1.
        straTIN.get_active_layer(maxlayh)
        # Compute multi-rock diffusion
        tdiff, diffQz, diffBe, erosion, deposition = flow.compute_sediment_hillslope(elevation, qprop, bedrock_Be, straTIN.alayR,
                                        diffcoeff*timestep, FVmesh.neighbours, FVmesh.vor_edges,
                                        maxlayh, FVmesh.edge_length, lGIDs)
        if input.btype == 'outlet':
            tdiff[flow.insideIDs[0],:] = 0.

        # Update dataset
        elevation += tdiff
        cumdiff += tdiff
        cumhill += tdiff
                
        # Update Quartz and 10Be content in deposited sediments (CP)
        sedBeIDs             = np.where(np.logical_and(diffQz > 0, tdiff > 0))[0]
        qprop[sedBeIDs]      = diffQz[sedBeIDs]
        bedrock_Be[sedBeIDs] = diffBe[sedBeIDs]

        # Update active layer
        straTIN.update_layers(erosion, deposition, elevation, verbose)

    if input.btype == 'slope':
        elevation[:len(flow.parentIDs)] = elevation[flow.parentIDs]-0.1
    elif input.btype == 'flat':
        elevation[:len(flow.parentIDs)] = elevation[flow.parentIDs]
    elif input.btype == 'wall':
        elevation[:len(flow.parentIDs)] = elevation[flow.parentIDs]+100.
    elif input.btype == 'outlet':
        elevation[1:len(flow.parentIDs)] = elevation[flow.parentIDs[1:]]+100.
    elif input.btype == 'wall1':
        elevation[:len(flow.parentIDs)] = elevation[flow.parentIDs]-0.1
        elevation[:recGrid.nx+1] = elevation[flow.parentIDs[:recGrid.nx+1]]+100.

    if verbose:
        print(" -   Get hillslope fluxes ", time.time() - walltime)

    # Update erodibility values
    if (input.erolays and input.erolays >= 0):
        mapero.getErodibility(elevation-oldelev)
        flow.erodibility = mapero.erodibility

    if applyDisp:
        elevation += disp * timestep

    tNow += timestep

    if verbose:
        print(" - Flow computation ", time.time() - flow_time)
    return tNow, elevation, cumdiff, cumhill, cumfail, bedrock_Be, qprop, slopeTIN
