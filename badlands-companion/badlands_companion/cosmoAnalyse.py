#--------------------------------------------------------#
# extract data on rock type and computes isotopic signature
# on a vertical line
#--------------------------------------------------------#
import h5py
import numpy
from badlands import cosmoProd

class cosmoAnalyse:
    """ Tools to extract information on erosion and deposition rates, geochemical 10Be signature in relation
    with surface exposure ages or marine sediment analyses (CP) """
    def __init__(self,folder,numsteps,dt):
        self.folder = folder
        self.nt = numsteps
        self.dt = dt

    def extract_incis_tot(self,point_dic):
        """ Extracts total incision on a point of the grid through time for comparison with surface exposure ages.
        Returns a dictionnary with total incision height (value) for each measurement point (key).

        point_incis
            output 2D array with time and incision height
        """
        point_incis = {}
        for point in point_dic.keys():
            print(point)
            point_incis[point] = numpy.zeros((self.nt,2))
        for i in range (0,self.nt):
            erotot    = numpy.zeros((self.nt,2))
            t         = i*self.dt
            filename  = self.folder+'tin.time'+str(i)+'.hdf5'
            f         = h5py.File(filename,'r')
            my_data   = f.get('cumdiff')
            my_data   = numpy.array(my_data)
            for point in point_dic.keys():
                pointnum    = point_dic[point]
                point_ero   = my_data[pointnum][0]
                erotot[i,0] = t
                erotot[i,1] = point_ero
                point_incis[point][i,:] = erotot[i,:]
            if f.__bool__():
                f.close()
        return point_incis

    def extract_rock_type(self,point):
        """ Extracts the proportion of each rock type in sediments for comparison with source estimates (epsilon Nd).
        To be used if different surface rocks have been defined. Made for 3 different rock sources.
        rt is a 4 columns output array with time, proportion of rock 1, rock2 and rock3 """
        rt         = numpy.zeros((self.nt,4))
        filename   = self.folder+'stratal.time'+str(self.nt)+'.hdf5'
        filename0  = self.folder+'stratal.time0'+'.hdf5'
        f          = h5py.File(filename,'r')
        f0         = h5py.File(filename0,'r')
        my_data1   = f.get('depoThickRock0')
        my_data1_0 = f0.get('depoThickRock0')
        my_data2   = f.get('depoThickRock1')
        my_data3   = f.get('depoThickRock2')
        my_data1   = numpy.array(my_data1)
        my_data1_0 = numpy.array(my_data1_0)
        my_data2   = numpy.array(my_data2)
        my_data3   = numpy.array(my_data3)
        r1         = my_data1[point,:]
        r1_0       = my_data1_0[point]
        r2         = my_data2[point,:]
        r3         = my_data3[point,:]
        n = -1
        for i in range (len(r1_0),len(r1)):
            n       = n+1
            t       = n*self.dt
            rt[n,0] = t
            rt[n,1] = r1[i]
            rt[n,2] = r2[i]
            rt[n,3] = r3[i]
        if f.__bool__():
            f.close()
        return rt

    def extract_Be_concentration(self,xmin,xmax,ymin,ymax,list_of_points=None):
        """ Extracts average 10Be concentration in deposited sediments with time on
        a square grid or on a list of points.

        NBe
            2D output array with time and average 10Be concentration per gram of Qz
        NBevar
            2D output array with time and standard deviation of 10Be concentration per gram of Qz
        """
        NBe    = numpy.zeros((self.nt,2))
        NBevar = numpy.zeros((self.nt,2))
        for i in range (0,self.nt):
            filename  = self.folder+'tin.time'+str(i)+'.hdf5'
            f         = h5py.File(filename,'r')
            my_data   = f.get('nbe')
            my_data   = numpy.array(my_data)
            cord      = f.get('coords')
            cord      = numpy.array(cord)
            qpr       = f.get('qprop')
            qpr       = numpy.array(qpr)
            x_data    = cord[:,0]
            y_data    = cord[:,1] 
            t         = i*self.dt
            NBe[i,0]  = t
            NBevar[i,0] = t
            NBetot = 0
            nbepts=0.
            if list_of_points is not None:
                pointIDs = numpy.array(list_of_points)
            else:
                pointIDs = numpy.where((x_data>xmin)&(x_data<xmax)&(y_data>ymin)&(y_data<ymax))[0]

            #isBe    = numpy.where(my_data[pointIDs]>0.)[0]
            isQ    = numpy.where(qpr[pointIDs]>0.)[0]

            if len(pointIDs[isQ])>0:
                #nbemean = numpy.mean(my_data[pointIDs[isQ]])
                #nbestd  = numpy.std(my_data[pointIDs[isQ]])
                nbemean = numpy.mean(my_data[pointIDs])
                nbestd  = numpy.std(my_data[pointIDs])
            else:
                nbemean = 0.
                nbestd = 0.

            NBe[i,1]    = nbemean
            NBevar[i,1] = nbestd

            if f.__bool__():
                f.close()

        return NBe,NBevar

    def extract_Be_erosion(self,xmin,xmax,ymin,ymax,list_bedrock=None):
        """Extract surface concentration in 10Be, erosion rate, quartz proportion, surface area
        and elevation with time on a square grid or on a list of bedrock points. 
        
        NBe
            2D output array with time and average 10Be concentration per gram of Qz
        erosion
            2D output array with time and total erosion (only) in meters
        quartz 
            2D output array with time and qz proportion in eroded rocks
        surface
            2D output array with time and the surface area
        topo 
            2D output array with time and elevation

        """
        NBe     = numpy.zeros((self.nt,2))
        erosion = numpy.zeros((self.nt,2))
        surface = numpy.zeros((self.nt,2))
        quartz  = numpy.zeros((self.nt,2))
        topo    = numpy.zeros((self.nt,2))
        for i in range (0,self.nt):
            NBe[i,0]     = i*self.dt
            erosion[i,0] = i*self.dt
            surface[i,0] = i*self.dt
            quartz[i,0]  = i*self.dt
            topo[i,0]    = i*self.dt
            filename  = self.folder+'tin.time'+str(i)+'.hdf5'
            f      = h5py.File(filename,'r')
            qpr    = f.get('qprop')
            cord   = f.get('coords')
            ero    = f.get('cumdiff')
            diffus = f.get('cumhill')
            Be     = f.get('nbe')
            areapx = f.get('area')
            qpr    = numpy.array(qpr)
            cord   = numpy.array(cord)
            ero    = numpy.array(ero)
            diffus = numpy.array(diffus)
            Be     = numpy.array(Be)
            areapx = numpy.array(areapx)
            x_data = cord[:,0]
            y_data = cord[:,1] 
            elev   = cord[:,2]
            if list_bedrock is not None:
                pointIDs = numpy.array(list_bedrock)
            else:
                pointIDs = numpy.where((x_data>xmin)&(x_data<xmax)&(y_data>ymin)&(y_data<ymax))[0]
            if len(pointIDs) > 0:
                NBe[i,1]       = numpy.mean(Be[pointIDs])
                erosion[i,1]   = numpy.mean(ero[pointIDs] + diffus[pointIDs])
                surface[i,1]   = numpy.mean(areapx[pointIDs])
                quartz[i,1]    = numpy.mean(qpr[pointIDs])
                topo[i,1]      = numpy.mean(elev[pointIDs])

            if f.__bool__():
                f.close()

        return NBe,erosion,quartz,surface,topo

    def extract_Be_deposition(self,lat,xmin,xmax,ymin,ymax,xminb,xmaxb,yminb,ymaxb,list_bedrock=None,list_of_points=None):
        """Extract the average 10Be concentration in deposited sediments (per g of Qtz) and the apparent erosion rate
        of the catchment with time (converted to mm.yr-1)

        depoBe 
            2D array with time and average 10Be concentration in deposited sediments (at per gram of Qz)
        apparentE
            2D array with time and apparent denudation rate from depoBe assuming steady state (in m.yr-1 of total rock)
        """
        depoBe = numpy.zeros(self.nt)
        apparentE = numpy.zeros(self.nt)
        cosmo     = cosmoProd.cosmoProd(lat)
        filename   = self.folder+'tin.time0.hdf5'
        f          = h5py.File(filename,'r')
        cord       = f.get('coords')
        cord       = numpy.array(cord)
        x_data     = cord[:,0]
        y_data     = cord[:,1] 

        if f.__bool__():
            f.close()

        # sediments points where to sample 10Be
        if list_of_points is not None:
            pointbasinIDs = numpy.array(list_of_points)
        else:
            pointbasinIDs = numpy.where((x_data>xminb)&(x_data<xmaxb)&(y_data>yminb)&(y_data<ymaxb))[0]

        # bedrock points where to compute production
        if list_bedrock is not None:
            pointbedrockIDs = numpy.array(list_bedrock)
        else:
            pointbedrockIDs = numpy.where((x_data>xmin)&(x_data<xmax)&(y_data>ymin)&(y_data<ymax))[0]

        for i in range(0,self.nt):
            filename   = self.folder+'tin.time'+str(i)+'.hdf5'
            f      = h5py.File(filename,'r')
            ero    = f.get('cumdiff')
            diffus = f.get('cumhill')
            Be     = f.get('nbe')
            qpr    = f.get('qprop')
            prod   = f.get('prate')
            ero    = numpy.array(ero)
            diffus = numpy.array(diffus)
            Be     = numpy.array(Be)
            qpr    = numpy.array(qpr)
            prod   = numpy.array(prod)
            x_data = cord[:,0]
            y_data = cord[:,1] 
            erotot = ero+diffus
            
            # extract mean 10Be concentration in the basin points
            isQ = numpy.where(qpr[pointbasinIDs]>0.)[0]
            if len(isQ) > 0:
                NBe = Be[pointbasinIDs[isQ]]
            else:
                NBe = 0
            mean_Be    = numpy.mean(NBe) 

            # extract mean 10Be concentration in bedrock points
            isQb = numpy.where(qpr[pointbedrockIDs]>0.)[0]
            if len(isQb) > 0:
                mean_prod = numpy.mean(prod[pointbedrockIDs[isQb]])
                mean_qp   = numpy.mean(qpr[pointbedrockIDs[isQb]])
            else:
                mean_prod = 0.
                mean_qp   = 0.5

            # Mean erosion rate (g.cm-2.yr-1)
            # Here production is per gram of Qtz per year and 10Be concentration also
            if mean_Be > 0:
                Be_ero = (mean_prod/mean_Be - cosmo.Be_lambda)*cosmo.L_n
                Be_ero = Be_ero / cosmo.rho / 100
            else:
                Be_ero = 0.
        
            apparentE[i] = Be_ero / mean_qp
            depoBe[i]    = mean_Be

            if f.__bool__():
                f.close()

        return depoBe,apparentE

    def extract_production (self,filename,xmin,xmax,ymin,ymax,list_of_points=None):
        """ Extract production rate (total of spallation and muon capture) on a square grid
        or on a list of points for a given time step.

        prodrate
            1D array for the list of points with their production rate (in grams of Qz per year)
        """
        f     = h5py.File(filename,'r')
        cord  = f.get('coords')
        prate = f.get('prate')
        cord  = numpy.array(cord)
        prate = numpy.array(prate)
        x_data = cord[:,0]
        y_data = cord[:,1]  

        if list_of_points is not None:
            pointIDs = numpy.array(list_of_points)
        else:
            pointIDs = numpy.where((x_data>xmin)&(x_data<xmax)&(y_data>ymin)&(y_data<ymax))[0]
        
        prodrate = prate[pointIDs]

        if f.__bool__():
            f.close()

        return prodrate

    def extract_erodedBe (self,lat,filename,xmin,xmax,ymin,ymax,list_of_points=None):
        """Extract the mean number of eroded and deposited 10Be atoms over a given area (square or list of points)
        for a given timestep.

        meanBe
            mean value of the 10Be concentration (atoms per gram of Qz) for the given points
        """
        cosmo     = cosmoProd.cosmoProd(lat)
        f         = h5py.File(filename,'r')
        cord      = f.get('coords')
        nBe       = f.get('nbe')
        eros      = f.get('cumdiff')
        diffus    = f.get('cumhill')
        qtz       = f.get('qprop')
        cord      = numpy.array(cord)
        nBe       = numpy.array(nBe)
        eros      = numpy.array(eros)
        diffus    = numpy.array(diffus)
        qtz       = numpy.array(qtz)
        x_data    = cord[:,0]
        y_data    = cord[:,1]  

        erotot    = numpy.zeros(len(nBe))
        depotot   = numpy.zeros(len(nBe))

        # eroded part
        erotot    = eros+diffus
        erotot[erotot > 0.] = 0.
        erodedQ  = -erotot * qtz
        erodedBe = erodedQ * nBe * 100 * cosmo.rho

        #d eposited part
        depotot   = eros+diffus
        depotot[depotot < 0.] = 0.
        depositedQ  = depotot * qtz
        depositedBe = depositedQ * nBe * 100 * cosmo.rho

        if list_of_points is not None:
            pointIDs = numpy.array(list_of_points)
        else:
            pointIDs = numpy.where((x_data>xmin)&(x_data<xmax)&(y_data>ymin)&(y_data<ymax))[0]

        isQ     = numpy.where(erodedQ[pointIDs]>0.)[0]
        Be      = erodedBe[pointIDs[isQ]]
        Qz      = erodedQ[pointIDs[isQ]]
        coord_ero = cord[pointIDs[isQ]]
        isdepoQ = numpy.where(depositedQ[pointIDs]>0.)[0]
        depoBe  = depositedBe[pointIDs[isdepoQ]]
        depoQz   = depositedQ[pointIDs[isdepoQ]]
        coord_depo = cord[pointIDs[isdepoQ]]

        if numpy.sum(Qz) > 0.:
            meanBe = numpy.sum(Be) / (numpy.sum(Qz)*100*cosmo.rho)
        else:
            meanBe = 0.

        if numpy.sum(depoQz) > 0.:
            meandepoBe = numpy.sum(depoBe) / (numpy.sum(depoQz)*100*cosmo.rho)
        else:
            meandepoBe = 0.    

        if f.__bool__():
            f.close()

        return Be,Qz,depoBe,depoQz,meanBe,meandepoBe,coord_ero,coord_depo,cord,eros+diffus

    def extract_shielding(self,filename):
        """Extract the shielding value for a given timestep.
        """
        f      = h5py.File(filename,'r')
        shield = f.get('topshield')
        shield = numpy.array(shield)
        if f.__bool__():
            f.close()
        return shield

    def extract_denud_Be_time(self,lat,xmin,xmax,ymin,ymax,list_of_points=None):
        """
        Extract the real erosion rate (m.yr-1) and the erosion rate from in-situ produced sediments
        in bedrock sources.
        realE
            Real average erosion rate (m.yr-1) with time computed from two consecutive time steps (1D array).
        apparentE
            Erosion rate (m.yr-1) computed from the 10BE concentration in eroded sediments with time (1D array).

        """
        cosmo     = cosmoProd.cosmoProd(lat)   
        realE     = numpy.zeros(self.nt)
        apparentE = numpy.zeros(self.nt)

        filename   = self.folder+'tin.time0.hdf5'
        f          = h5py.File(filename,'r')
        cord       = f.get('coords')
        cord       = numpy.array(cord)
        x_data     = cord[:,0]
        y_data     = cord[:,1] 

        if list_of_points is not None:
            pointIDs = numpy.array(list_of_points)
        else:
            pointIDs = numpy.where((x_data>xmin)&(x_data<xmax)&(y_data>ymin)&(y_data<ymax))[0]

        if f.__bool__():
            f.close()

        for i in range (self.nt):
            filename = self.folder+'tin.time'+str(i)+'.hdf5'
            f      = h5py.File(filename,'r')
            eros   = f.get('cumdiff')
            diffus = f.get('cumhill')
            Be     = f.get('nbe')
            qtz    = f.get('qprop')
            prod   = f.get('prate')
            eros   = numpy.array(eros)
            diffus = numpy.array(diffus)
            Be     = numpy.array(Be)
            qtz    = numpy.array(qtz)
            prod   = numpy.array(prod)
            x_data = cord[:,0]
            y_data = cord[:,1] 
            erotot = eros+diffus
            erotot[erotot > 0.] = 0.
    
            # how much has been eroded
            if i == 0:
                eroded = erotot
            else:
                eroded = erotot - erotot_old
    
            erotot_old = numpy.copy(erotot)
    
            # erosion rate in mm/yr
            erosion_rate = eroded/self.dt
    
            # mass of eroded sediments in g
            m_eroded = -eroded * 100 * cosmo.rho
    
            # mass of eroded quartz in g
            mqtz_eroded = m_eroded * qtz
    
            # number of atoms of 10Be eroded
            nBe_eroded  = mqtz_eroded * Be
    
            # mean production rate of 10Be per gram of quartz
            isQ     = numpy.where(qtz[pointIDs] > 0.)[0]
            meanQ   = numpy.mean(mqtz_eroded[pointIDs[isQ]])
            meanqpr = numpy.mean(qtz[pointIDs[isQ]])
            if meanQ > 0.:
                mean_Be   = numpy.mean(nBe_eroded[pointIDs[isQ]])/meanQ
            else:
                mean_Be = 0.
    
            #isprod    = numpy.where(prod[pointIDs]>0.)[0]
            mean_prod = numpy.mean(prod[pointIDs[isQ]])
    
            # mean average erosion rate
            mean_ero  = numpy.mean(erosion_rate[pointIDs])
    
            # mean erosion rate (g.cm-2.yr-1)
            if mean_Be > 0:
                Be_ero = (mean_prod/mean_Be - cosmo.Be_lambda)*cosmo.L_n
                Be_ero = Be_ero / cosmo.rho / 100
            else:
                Be_ero = 0.
        
            realE[i]     = -mean_ero
            apparentE[i] = Be_ero/meanqpr

            if f.__bool__():
                f.close()

        return realE, apparentE
