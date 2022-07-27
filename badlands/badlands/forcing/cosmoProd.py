##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
##                                                                                   ##
##  This file forms part of the Badlands surface processes modelling application.    ##
##                                                                                   ##
##  For full license and copyright information, please refer to the LICENSE.md file  ##
##  located at the project root, or contact the authors.                             ##
##                                                                                   ##
##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
"""
This module contains parameters for 10Be production rate computation depending on 
geographical and geometrical parameters as in Norton et al. 2009 and Stone et al. 2000
(Carole Petit)
"""
import math
import numpy

class cosmoProd:
    
    """
    This class  contains tools needed for computing TCN (10Be) evolution
    N_dic stores 10Be production from 3 different mechanisms
    Nn is for neutron spallation
    Nf is for fast muon capture
    Ns is for slow muon capture
    """
    N_dic = {'Nnval':[],'Nfval':[],'Nsval':[]}

    def __init__(self,lati):
        """
        Initialization.
        Be_lambda    = radioactive decay constant of 10Be in yr-1
        rho          = rock density in g. cm-2
        erosion_rate = in g. yr-1
        a,b,c,d,e are atmospheric pressure coefficients depending on latitude
        Production rates for spallation, slow and fast muons capture
        L_n,L_f and Ls are attenuation lengths for neutrons, fast and slow muons in g. cm-2
        see parameters in Stone et al. 2000 and Braucher et al. 2011
        Coeff is used to calibrate production rates according to Mariotti et al., 2019
        """
        self.spall_prod_rate = 4.11
        self.slow_prod_rate  = 0.011
        self.fast_prod_rate  = 0.039
        self.abcde = {0:[31.8518,250.3193,-0.083393,7.426e-5,-2.2397e-8,0.587],10:[34.3699,258.4759,-0.089807,\
                     7.9457e-5,-2.3697e-8,0.6],20:[40.3153,308.9894,-0.106248,9.4508e-5,-2.8234e-8,0.678],\
                     30:[42.0983,512.6857,-0.120551,1.1752e-4,-3.8809e-8,0.833],40:[56.7733,649.1343,-0.160859,1.5463e-4,-5.0330e-8,0.933],\
                     50:[69.0720,832.4566,-0.199252,1.9391e-4,-6.3653e-8,1],60:[71.8733,863.1927,-0.207069,2.0127e-4,-6.6043e-8,1]}
        self.f_slow  = 0.012
        self.f_fast  = 0.0065
        self.f_spall = 1 - self.f_slow - self.f_fast

        self.L_n = 160.
        self.L_f = 4320.
        self.L_s = 1500.
        self.Be_lambda = 4.9867E-7
        self.rho = 2.7
        self.lat = lati 
        self.coeff = 1.2

        return

    def Be_production_rate(self, elevation, tshield, qprop):
        """
        Computes 10Be production rate depending on latitude and altitude 
        and topographic shielding
        TCN production is in atoms per gram of Qtz per year
        """
        numpts = len(elevation)
        p_atm = numpy.zeros(numpts)
        prodn = numpy.zeros(numpts)
        prods = numpy.zeros(numpts)
        prodf = numpy.zeros(numpts)
        P_rate = numpy.zeros(numpts)

        lat10 = self.lat//10*10
        if lat10 > 60:
           lat10 = 60

        a = float(self.abcde[lat10][0])
        b = float(self.abcde[lat10][1])
        c = float(self.abcde[lat10][2])
        d = float(self.abcde[lat10][3])
        e = float(self.abcde[lat10][4])
        m = float(self.abcde[lat10][5])

        p_atm = 1013.25*numpy.exp(-0.03417/0.0065*(numpy.log(288.15)-numpy.log(288.15-0.0065*elevation)))
        p_atm[numpy.where(p_atm>1013.25)[0]]=1013.25
        slambda = a + b*numpy.exp(-p_atm/150) + c*p_atm + d*p_atm**2 + e*p_atm**3
        mlambda = m*numpy.exp((1013.25-p_atm)/242)

        notQIDs= numpy.where(qprop<=0)[0]
        prodn = self.f_spall*slambda*self.spall_prod_rate * tshield * self.coeff
        prods = self.f_slow*mlambda*self.slow_prod_rate   * tshield * self.coeff
        prodf = self.f_fast*mlambda*self.fast_prod_rate   * tshield * self.coeff
        prodn[notQIDs]=0
        prods[notQIDs]=0
        prodf[notQIDs]=0
        P_rate = prodn + prodf + prods

        return  prodn,prods,prodf,P_rate

    def topography_shielding(self, fact, iceT, sealevel, elevation, neighbours, edge_length):
        
        """
        Shielding function computes the angle of topographic obstruction
        and the corresponding shielding factor for cosmonuclide production rate.

        Parameters
        ----------
        real : fact
            Parameter (<=1) that corrects shielding underestimation due to DEM grid smoothing

        variable : elevation
            Numpy arrays containing the elevation of the TIN nodes.

        variable : iceT
            Numpy arrays containing the ice thickness.

        variable : neighbours
            Numpy integer-type array with the neighbourhood IDs.

        variable : edges
            Numpy real-type array with the voronoi edges length for each neighbours of the TIN nodes.

        """

        # compute angle and shielding correction for the topography
        numpts      = len(elevation)
        shield_mean = numpy.zeros(numpts)

        # points with non-null shielding
        IDs_shield = numpy.where(numpy.logical_and(iceT<20,elevation>sealevel))[0]

        # for id in range(numpts):
        for id in IDs_shield:

            distance = edge_length[id,:]
            ngbhs    = neighbours[id,:]

            # find number of neighbours
            nIDs   = numpy.where (distance > 0.1)[0]
            lngbhs = len(nIDs)
            shield_angle = 0.
            # assumes equal repartition of neighbours in space
            mean_phi     = 360./lngbhs
            shielding    = 0.

            for n in range(lngbhs):
                if elevation[id] < elevation [ngbhs[n]] :
                    shield_angle = math.atan((elevation[ngbhs[n]]-elevation[id])/edge_length[id,n])
                    shielding   += mean_phi/360.*(math.sin(shield_angle))**3.3

            # correction to shielding due to DEM smoothing (amplifies shielding)
            temp            = ((1.-shielding)*fact)**3
            shield_mean[id] = temp

        return shield_mean   

    def compute_Be_concentration (self,seal,elev,tshield,qprop,nbe,cero,tstep):
        """
        Function computing  10Be concentration at the surface due to cosmic
        ray exposure duration in at. g-1

        Parameters
        ----------
        elev
            Numpy arrays containing the elevation of the TIN nodes.

        tshield
            Numpy array with the topographic shielding factor 

        qprop
            Numpy float array indicating the quartz proportion.

        nbe
            Numpy float array with the number of atoms of 10Be per gram of Qtz

        cero
            Numpy real-type array with the amount of erosion at the given timestep.

        """
        numpts  = len(elev)
        # computes topographic shielding angle and 10Be production rate
        prod_raten,prod_rates,prod_ratef,prod_rate = self.Be_production_rate(elev,tshield,qprop)

        # computes 10Be concentration at the surface on each node
        # use Eulerian formulation see paper by Knudsen et al. 2019 Quater. Geol.
        # TCN concentration is in atoms per gram of Qtz so we convert erosion rate
        # taking into account the Qtz concentration

        Nn = self.N_dic['Nnval']
        Nf = self.N_dic['Nfval']
        Ns = self.N_dic['Nsval']

        for i in range(1,numpts):
            erosion = -cero[i]*self.rho*100*qprop[i]/tstep
            An = (self.Be_lambda + erosion/self.L_n)*tstep
            Af = (self.Be_lambda + erosion/self.L_f)*tstep
            As = (self.Be_lambda + erosion/self.L_s)*tstep
            if erosion > 0 and elev[i] >= seal :
               Nn[i] = numpy.exp(-An) * (Nn[i] + prod_raten[i] * 1/(self.Be_lambda + erosion/self.L_n) * (numpy.exp(An)-1))
               Nf[i] = numpy.exp(-Af) * (Nf[i] + prod_ratef[i] * 1/(self.Be_lambda + erosion/self.L_f) * (numpy.exp(Af)-1))
               Ns[i] = numpy.exp(-As) * (Ns[i] + prod_rates[i] * 1/(self.Be_lambda + erosion/self.L_s) * (numpy.exp(As)-1))
            else:
               Nn[i] = nbe[i] * numpy.exp(-self.Be_lambda*tstep) 
               Nf[i] = 0
               Ns[i] = 0

        N_Be =  Nn + Ns + Nf
        self.N_dic['Nnval']=Nn
        self.N_dic['Nfval']=Nf
        self.N_dic['Nsval']=Ns

        return N_Be, prod_rate

    def compute_init_Be(self,seal,elev,tshield,qprop,erorate):
        """
        Function computing  10Be concentration at the surface due to cosmic
        ray exposure duration in at. g-1 at steady state for initial values

        Parameters
        ----------
        elev
            Numpy arrays containing the elevation of the TIN nodes.

        tshield
            Numpy array with the topographic shielding factor 

        qprop
            Numpy float array indicating the concentration of quartz.

        erorate        
            float = steady-state erosion rate in m.yr-1.
            It is converted into mass (of Qtz) loss per year (erosion_rate)

        N_Be
            output 10Be concentration in atoms per gram of quartz.

        prod_rate
            output 10Be production rate in atoms per gram of quartz per year.

        """
        numpts  = len(elev)
        N_Be    = numpy.zeros(numpts)
        prod_raten,prod_rates,prod_ratef,prod_rate = self.Be_production_rate(elev,tshield,qprop)

        # computes 10Be concentration at the surface on each node
        # in g.cm-2.a-1

        # eliminates interpolation effect on the borders of the quartz map
        qmax          = numpy.max(qprop)
        erosion_rate  = erorate*100*self.rho*qmax

        Nn = prod_raten / ( erosion_rate/self.L_n + self.Be_lambda)
        Nf = prod_ratef / ( erosion_rate/self.L_f + self.Be_lambda)
        Ns = prod_rates / ( erosion_rate/self.L_s + self.Be_lambda)
        N_Be = Nn + Nf + Ns
        self.N_dic['Nnval']=Nn
        self.N_dic['Nfval']=Nf
        self.N_dic['Nsval']=Ns

        return N_Be, prod_rate

    def update_Qz_Be(self,nbe,qpr):
       """ Update Quartz map removing nodes under a quartz concentration threshold
       """
       noQIDs     = numpy.where(qpr < 0.005)
       qpr[noQIDs]= 0.
       nbe[noQIDs]= 0.
       return nbe,qpr
