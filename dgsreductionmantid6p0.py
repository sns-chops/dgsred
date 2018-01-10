#!/usr/bin/env python

import sys
import os,stat
import imp
from ARLibrary import Load_2_Monitors
from XMLparser import XMLparser
from mantid import *
from mantid.simpleapi import *
from numpy import *
import mantid.api
from string import *


def GetPathFromRunNumber(instrument,run):
    if instrument in ['ARCS','CNCS','SEQUOIA','HYSPEC']:
        try:
            alg = mantid.api.AlgorithmManager.createUnmanaged("Load")
            alg.initialize()
            alg.setPropertyValue('Filename',instrument+str(run))
            fname = alg.getProperty('Filename').value
            return fname
        except:
            raise ValueError("Event Nexus file not found:  " + instrument + str(run))
    else:
        raise ValueError("Instrument not yet implemented")    
      
      
class dgsreduction(object):

    def __init__(self, XMLfile=None):
        
        
        if XMLfile!=None:
            if not(os.path.isfile(XMLfile)):
                raise IOError ("data text file "+ XMLfile+ " not found")
            # These three text variables are for the summary file.
            self.loadtext=''
            self.vantext=''
            self.reductiontext=''
            self.RunFromXML(XMLfile)    
        
      
      
    def loadvan(self,parsed):        

        if parsed.calibdict['VanPath'] == "":
            parsed.calibdict['VanPath'] = os.getcwd()
    
        vanfileout = os.path.join(parsed.calibdict['VanPath'],parsed.calibdict['SaveProcDetVanFilename'])
        
        #print "*****"
        #print vanfileout
        
        #After parsing, and before treating the data, need to examine the vanadium sensitivity correction
        if os.path.isfile(vanfileout):
            Load(Filename=vanfileout,OutputWorkspace="__VAN")
            dictvan={'UseProcessedDetVan':'1','DetectorVanadiumInputWorkspace':'__VAN'}
            self.vantext+='Using processed vanadium file '+vanfileout + '\n'

            
        else:
            #add the vanpath to the data search path
            #parse the string of VanRuns ans instrument_run1+instrument_run2+etc.
            #  
            # Make the van file
            # add van runs to the path
            config.appendDataSearchDir(parsed.calibdict['VanPath'])
            
            #Transform the list of ints (for the VanRuns) into a list of strings of the format instrument_run1+instrument_run2+etc.
            vanrunstring = parsed.calibdict['Instrument']+"_"+join(str(parsed.calibdict['VanRuns']).strip("[").strip("]").replace(" ","").split(","),"+"+parsed.calibdict['Instrument']+"_")
            Load(Filename=vanrunstring,OutputWorkspace="__VAN")
            self.vantext+='Loading vanadium files '+vanrunstring+ '\n'

            if parsed.datadictsother[0]['FilterBadPulses']:
                FilterBadPulses(InputWorkspace = '__VAN', OutputWorkspace = '__VAN')
                self.vantext+='Bad pulses filtered from vanadium files\n'


            #Do the Masking.
            for elem in parsed.calibdict['Mask']:
                if elem.has_key('algorithm') and elem['algorithm'].lower()=='maskbtp':
                    del elem['algorithm']
                    MaskBTP(Workspace="__VAN",**elem)
                    self.vantext+='Vanadium file masked with MaskBTP('+str(elem)+')\n'


                if elem.has_key('algorithm') and elem['algorithm'].lower()=='maskangle':
                    del elem['algorithm']
                    MaskAngle(Workspace='__VAN',**elem)
                    self.vantext+='Vanadium file masked with MaskAngle('+str(elem)+')\n'

             
            dictvan={'SaveProcessedDetVan':'1','DetectorVanadiumInputWorkspace':'__VAN','SaveProcDetVanFilename':vanfileout}
            self.vantext+='Vanadium calibration file saved as '+vanfileout+'\n'

        return dictvan

    def RunFromXML(self,filename):      
        parsed=XMLparser(filename)
        
        #The output of XMLparse contains at least two arrays of dictionarys, one for dgsreduction keywords and one for other keywords
        # datadictsdgs
        # datadictsother    
        

        for d in range(len(parsed.datadictsdgs)):

            #initialize the text variables.
            self.loadtext=''
            self.vantext=''
            self.reductiontext=''


            #Deal with the vanadium sensitivity correction
            vandict = self.loadvan(parsed)


            #determine the scantype
            #single OR sweep
            if parsed.datadictsother[d]['ScanType'] == 'single' or  parsed.datadictsother[d]['ScanType'] == 'sweep':
                runno = int(parsed.datadictsother[d]['Runs'][0])
                FileName = parsed.datadictsother[d]['Instrument']+"_"+str(runno)
           
                #Load each file and fix the time-series to start at 'zero'
                data = Load(Filename=FileName)
                if runno in range(141471, 141475+1):
                    ChangeBinOffset(InputWorkspace="data",OutputWorkspace="data",Offset=882)
                self.loadtext+='Data files loaded '+FileName+'\n'
                
                path = data.getRun()['Filename'].value
                #do the correction for log times
                CorrectLogTimes('data')
                self.loadtext+='Data files corrected for log times\n'
                
                monitors = Load_2_Monitors(path)
                monitors = monitors.rename()

                if parsed.datadictsother[d]['FilterBadPulses']:
                    data = FilterBadPulses(InputWorkspace =data)
                    self.loadtext+="Bad pulses filtered from data files\n"
            
                #filter by additional log values.
                if parsed.datadictsother[d]['FilterNames'] != None:
                    for cntr,part in enumerate(parsed.datadictsother[d]['FilterNames']):
                        data = FilterByLogValue(InputWorkspace = 'data', LogName=part, 
                                MinimumValue=parsed.datadictsother[d]['FilterMin'][cntr],
                                MaximumValue=parsed.datadictsother[d]['FilterMax'][cntr],TimeTolerance=0,LogBoundary='Left')

                        self.loadtext += "Data filtered by "+part+" between "+str(MinimumValue)+" and "+str(MaximumValue)+".\n"
                        #print "Data filtered by "+part+" between "+str(parsed.datadictsother[d]['FilterMin'][cntr])+" and "+str(parsed.datadictsother[d]['FilterMax'][cntr])+".\n"
                


                #now deal with all the other runs, if there are more than one.
                if len(parsed.datadictsother[d]['Runs']) > 1:
                    for i in range(1,len(parsed.datadictsother[d]['Runs'])):
                        runno = int(parsed.datadictsother[d]['Runs'][i])
                        FileName = parsed.datadictsother[d]['Instrument']+"_"+str(runno)
                        #Load each file and fix the time-series to start at 'zero'
                        datatemp = Load(Filename=FileName)
                        if runno in range(141471, 141475+1):
                            ChangeBinOffset(InputWorkspace="datatemp",OutputWorkspace="datatemp",Offset=882)
                        path = datatemp.getRun()['Filename'].value
                        self.loadtext+='Data files loaded '+FileName+'\n'


                        #Fix all of the time series log values to start at the same time as the proton_charge
                        #do the correction for log times
                        CorrectLogTimes('datatemp')
                        self.loadtext+='Data files corrected for log times\n'
                        monitorstemp = Load_2_Monitors(path)
                        monitors += monitorstemp

                        if parsed.datadictsother[d]['FilterBadPulses']:
                            datatemp = FilterBadPulses(InputWorkspace =datatemp)
                            self.loadtext+="Bad pulses filtered from data files\n"
                            
                        #filter by additional log values.
                        if parsed.datadictsother[d]['FilterNames'] != None:
                            for cntr,part in enumerate(parsed.datadictsother[d]['FilterNames']):
                                datatemp = FilterByLogValue(InputWorkspace = 'datatemp', LogName=part, 
                                MinimumValue=parsed.datadictsother[d]['FilterMin'][cntr],
                                MaximumValue=parsed.datadictsother[d]['FilterMax'][cntr],TimeTolerance=0,LogBoundary='Left')
                                self.loadtext += "Data filtered by "+part+" between "+str(MinimumValue)+" and "+str(MaximumValue)+".\n"
                                
                        data += datatemp
                        
                        self.loadtext+='Data added to the previous workspace\n'



                #This is where the reduction is done.
                if  parsed.datadictsother[d]['ScanType'] == 'single':
                    self.ProcessWorkspace(data,monitors,parsed.datadictsdgs[d],parsed.datadictsother[d],vandict)
                else:
                    #split up the sweep by the sweep variable.
                    
                    logvalue = parsed.datadictsother[d]['LogValue']
                    logvaluemin = parsed.datadictsother[d]['LogValueMin']
                    logvaluemax = parsed.datadictsother[d]['LogValueMax'] 
                    logvaluestep = parsed.datadictsother[d]['LogValueStep']
                                      
                    #Check if the logvalue has been set                    
                    if logvalue == None or data.run().hasProperty(logvalue)==False:
                        raise ValueError("No LogValue given OR the given log value was not found in the file.")

                    #need to split the data by an independt variable , some log value.
                    #Create the array of logvalue BOUNDARIES
                    if logvaluemin == None:
                        logvaluemin= array(data.run().getProperty(logvalue).value).min()
                    if logvaluemax == None:
                        logvaluemax= array(dat.run().getProperty(logvalue).value).max()
                    if logvaluestep == None:
                        logvaluestep = logvaluemax - logvaluemin

                    bounds = arange(float(logvaluemin), float(logvaluemax)+float(logvaluestep), float(logvaluestep))

                    #Get the time correlation correct if you set the time correlation keyword.
                    #To first approximation, set the time to zero for the first.
                    for i in range(len(bounds)-1):
                        dataslice = FilterByLogValue(InputWorkspace=data, LogName= logvalue,MinimumValue=float(bounds[i]) ,MaximumValue = float(bounds[i+1]))
                        if dataslice.getNumberEvents()>0:
                            values=array(dataslice.run().getProperty(logvalue).value)
                            self.reductiontext= "Processing data for "+logvalue+" between "+str(bounds[i])+" and "+str(bounds[i+1])+", mean="+str(values.mean())+" std="+str(values.std())+"\n"
                            self.ProcessWorkspace(dataslice,monitors,parsed.datadictsdgs[d],parsed.datadictsother[d],vandict)

                  
            
            if parsed.datadictsother[d]['ScanType'] == 'step':

                for currentrun in parsed.datadictsother[d]['Runs']:
                
                    FileName = parsed.datadictsother[d]['Instrument']+"_"+str(currentrun)
           
                
                    #Load each file and fix the time-series to start at 'zero'
                    data = Load(Filename=FileName)
                    self.loadtext='Data files loaded '+FileName+'\n'
                    path = data.getRun()['Filename'].value
                    #do the correction for log times
                    CorrectLogTimes('data')
                    monitors = Load_2_Monitors(path)
                    self.loadtext+='Data files corrected for log times\n'

                    if parsed.datadictsother[d]['FilterBadPulses']:
                        data = FilterBadPulses(InputWorkspace =data)
                        self.loadtext+="Bad pulses filtered from data files\n"
                    #filter by additional log values.
                    if parsed.datadictsother[d]['FilterNames'] != None:
                        for cntr,part in enumerate(parsed.datadictsother[d]['FilterNames']):
                            data = FilterByLogValue(InputWorkspace = 'data', LogName=part, 
                                    MinimumValue=parsed.datadictsother[d]['FilterMin'][cntr],
                                    MaximumValue=parsed.datadictsother[d]['FilterMax'][cntr],TimeTolerance=0,LogBoundary='Left')
                            self.loadtext += "Data filtered by "+part+" between "+str(MinimumValue)+" and "+str(MaximumValue)+".\n"
                            
                    self.reductiontext = ''
                    self.ProcessWorkspace(data,monitors,parsed.datadictsdgs[d],parsed.datadictsother[d],vandict)

    
  
    



    def ProcessWorkspace(self,data,monitors,datadictsdgs,datadictsother,vandict):
        #merge the datadictsdgs with the vandict dictionary.
        datadictsdgs.update(vandict)

        if datadictsdgs['IncidentEnergyGuess'] == None:
            if data.getRun().hasProperty('EnergyRequest'):
                datadictsdgs['IncidentEnergyGuess'] = data.getRun().getProperty('EnergyRequest').getStatistics().mean
            else:
                raise ValueError("no IncidentEnergyGuess has been set, and no EnergyRequest was found in file.")

        #if we have event monitors we need to deal with them seperately
        monitorws=Rebin(monitors,"1",PreserveEvents=False)

        #Time-ind-bg subtraction.
        if (datadictsdgs['TimeIndepBackgroundSub']):
            #check if tibmin and tibmax have been defined.
            tibmin = datadictsdgs['TibTofRangeStart']
            tibmax = datadictsdgs['TibTofRangeEnd']
            if tibmin== None or tibmax == None:
                if data.getInstrument().getName() == 'CNCS':
                    [tibmin,tibmax] = SuggestTibCNCS(datadictsdgs['IncidentEnergyGuess'])
                elif data.getInstrument().getName() == 'HYSPEC':
                    [tibmin,tibmax] = SuggestTibHYSPEC(datadictsdgs['IncidentEnergyGuess'])
                else:                
                    raise ValueError("Time independent background subtraction selected, but no limits set.")


        #Generate the Grouping file
        if datadictsdgs['GroupingFile'] != None:
            datadictsdgs['GroupingFile'] = self.checkgrouping(data,datadictsdgs,datadictsother)


        #Do the DGSreduction command  **datadictsdgs fills in all the requsted keyword values
        #Here is where the reduction is actually done.
        out = DgsReduction(SampleInputWorkspace=data,SampleInputMonitorWorkspace=monitorws,**datadictsdgs)
        out = out[0]


        self.reductiontext += 'DgsReduction was called with the following parameters\n'
        for key, value in datadictsdgs.items():
            self.reductiontext += '\t' + key + ':  ' + str(value) + '\n'
            
            
            
        if datadictsdgs['UseIncidentEnergyGuess'] == True :
            tempstring="(Fixed)"
        else:
            tempstring="(Calculated)"
 
        self.reductiontext+='Incident energy is '+str(out.getRun()['Ei'].value)+' meV  '+tempstring+'\n'
        self.reductiontext+='Emision time is '+str(out.getRun()['CalculatedT0'].value)+' microseconds  '+tempstring+'\n'
        
        
        
        
        if datadictsdgs['EnergyTransferRange'] == '' or datadictsdgs['EnergyTransferRange']==None:
            Eguess = out.getRun()['EnergyRequest'].value
            if datadictsdgs['UseIncidentEnergyGuess'] == True:
                Eguess = out.getRun()['Ei'].value
            self.reductiontext += 'EnergyTransferRange automatically chosen (min, step, max):  ' + '-' +str(0.5*Eguess)+', '+ str(0.01*Eguess)+', '+ str(0.99*Eguess) +'\n'
        else:
            self.reductiontext += 'EnergyTransferRange chosen by user (min, step, max):  ' +datadictsdgs['EnergyTransferRange'] +'\n'
        
        totalmuAhr = out.run().getProtonCharge()
        totalcoul  = totalmuAhr/1000*3.6
        
        self.reductiontext += 'Accumulated proton charge: '+str(totalmuAhr) + ' (micro-Amp-hours), '+str(totalcoul)+'  (Coulombs)\n'


	    #deal with angles
        [psiangle, angletext]= definegoniometer(datadictsother['GoniometerMotor'],datadictsother['GoniometerMotorOffset'],datadictsother['GoniometerMotorDirection'],datadictsother['GoniometerMotorAxis'],out)


        self.reductiontext += angletext + '\n'


        #This if statement will normalize the vanadium sensitivity correction to fluctuate around 1.0
        if datadictsother['NormalizedCalibration'] and not (datadictsdgs.has_key('UseProcessedDetVan') and datadictsdgs['UseProcessedDetVan'] == '1'):   
            LoadNexus(Filename=datadictsdgs['SaveProcDetVanFilename'],OutputWorkspace="__VAN")
            datay = mtd['__VAN'].extractY()
            meanval = float(datay[datay>0].mean())
            CreateSingleValuedWorkspace(OutputWorkspace='__meanval',DataValue=meanval)
            Divide(LHSWorkspace='__VAN',RHSWorkspace='__meanval',OutputWorkspace='__VAN')  #Divide the vanadium by the mean
            out  = Multiply(LHSWorkspace=out,RHSWorkspace='__meanval') #multiple by the mean of vanadium Normalized data = Data / (Van/meanvan) = Data *meanvan/Van
            SaveNexus(InputWorkspace="__VAN", Filename= datadictsdgs['SaveProcDetVanFilename'])
            self.vantext+='Vanadium normalization file has been normalized to fluctuate about 1.0\n'

            


#filetypesexist = ['phx-DONE', 'spe-DONE', 'nxspe-DONE', 'par-DONE', 'jpg-DONE', 'nxs-DONE', 'mdnxs-halfway', 'iofq','iofe',
#                            'iofphiecolumn','iofphiearray','iofqecolumn','iofqearray','summary',
#                           'sqw','vannorm']
  
         #Now deal with saving files
        if datadictsother['Save'] != None:
            #create a friendly name
            friendlynamebase = self.CreateFriendlyFilename(out,datadictsother)
            if 'summary' in datadictsother['Save']:
                summaryfilename = friendlynamebase+"_summary.txt"
                parentdir = os.path.dirname(summaryfilename)    
                if os.path.isdir(parentdir) == False:
                    os.mkdir(parentdir,0755)
                sumfile = open(summaryfilename, 'w')
                sumfile.write("-----------VANADIUM CALIBRATION AND MASKING-----------\n")
                sumfile.write(self.vantext)
                sumfile.write("\n-----------DATA-----------------------------\n")
                sumfile.write(self.loadtext)
                sumfile.write(self.reductiontext)
                sumfile.close()
                changepermissions(summaryfilename)

            #Parse the filetypes to save
            if 'nxspe' in datadictsother['Save']:
                efixed = out.getRun()['Ei'].value
                #if there is a powder maping file in the current directory, then use it with the .nxspe file
                if datadictsdgs['GroupingFile'] == os.path.join(os.path.abspath(os.curdir),'powder.xml'):
                    SaveNXSPE(Filename=friendlynamebase+".nxspe", InputWorkspace=out, Efixed=str(efixed),Psi=0.0, KiOverKfScaling=datadictsdgs['CorrectKiKf'], ParFile=os.path.join(os.path.abspath(os.curdir),'powder.par'))
                else:
                    SaveNXSPE(Filename=friendlynamebase+".nxspe", InputWorkspace=out, Efixed=str(efixed),Psi=str(psiangle), KiOverKfScaling=datadictsdgs['CorrectKiKf'])


                print friendlynamebase

                #change permissions of the directory and file
                changepermissions(friendlynamebase+".nxspe")
 
            if 'nxs' in datadictsother['Save']:
                #save the nxs
                SaveNexus(Filename=friendlynamebase+".nxs", InputWorkspace=out)
                #change permissions of the directory and file
                changepermissions(friendlynamebase+".nxs")


            if 'par' in datadictsother['Save']:
                if datadictsdgs['GroupingFile'] == os.path.join(os.path.abspath(os.curdir),'powder.xml'):
                    shutil.copy(os.path.abspath(os.curdir) +"/powdergroup.par",friendlynamebase+".par")
                else:
                    SavePAR(Filename=friendlynamebase+".par", InputWorkspace=out)
                #change permissions of the directory and file
                changepermissions(friendlynamebase+".par")

            if 'phx' in datadictsother['Save']:
                SavePHX(Filename=friendlynamebase+".phx",InputWorkspace=out)
                changepermissions(friendlynamebase+".phx")

            if 'spe' in datadictsother['Save']:
                SaveSPE(Filename=friendlynamebase+".spe",InputWorkspace=out)
                changepermissions(friendlynamebase+".spe")

            if 'jpg' in datadictsother['Save']:
                #plot sqw
                sys.path.insert(0,"/mnt/software/lib/python2.6/site-packages/matplotlib-1.2.0-py2.6-linux-x86_64.egg/")
                import  matplotlib
                matplotlib.use("agg")
                import  matplotlib.pyplot as plt



                #plots
                # Update ConvertToMDHelper to new algorithm name per mandtid changeset 9396 - Ricardo 2015-06-25
                # minvals,maxvals=ConvertToMDHelper(out,'|Q|','Direct')
                minvals,maxvals=ConvertToMDMinMaxGlobal(out,'|Q|','Direct')
                xmin=minvals[0]
                xmax=maxvals[0]
                xstep=(xmax-xmin)*0.01
                ymin=minvals[1]
                ymax=maxvals[1]
                ystep=(ymax-ymin)*0.01
                x=arange(xmin,xmax,xstep)[0:100]
                y=arange(ymin,ymax,ystep)[0:100]
                X,Y=meshgrid(x,y)

                MD=ConvertToMD(out,QDimensions='|Q|',dEAnalysisMode='Direct',MinValues=minvals,MaxValues=maxvals)
                ad0='|Q|,'+str(xmin)+','+str(xmax)+',100'
                ad1='DeltaE,'+str(ymin)+','+str(ymax)+',100'
                MDH=BinMD(InputWorkspace=MD,AlignedDim0=ad0,AlignedDim1=ad1)
                d=MDH.getSignalArray()
                ne=MDH.getNumEventsArray()
                dne=d/ne

                ## Save the plot in raw in adition to the jpg - Ricardo 2015-06-25
                dne_no_nan = np.nan_to_num(dne)
                f = open(str(friendlynamebase+"_sqw_2d.dat"),'w')
                f.write('#X\tY\tZ\tE\n')
                for xidx, xi in enumerate(x):
                    for yidx, yi in enumerate(y):
                        f.write("%f\t%f\t%f\t%f\n"%(xi,yi,dne_no_nan[xidx,yidx],sqrt(dne_no_nan[xidx,yidx]) ) )
                f.close()

                ## Save the plot in 1D raw in adition to the jpg - Ricardo 2015-06-25
                out2 = SumSpectra(InputWorkspace='out')
                s1d=SumSpectra(out2)
                x1d=s1d.readX(0)
                y1d=s1d.readY(0)
                f = open(str(friendlynamebase+"_sqw_1d.dat"),'w')
                f.write('#X\tY\tE\n')
                for xi,yi in zip(x1d[1:],y1d):
                    f.write("%e\t%e\t%e\n"%(xi,yi,sqrt(yi) ) )
                f.close()
                DeleteWorkspace(out2)


                Zm=ma.masked_where(ne==0,dne)
                plt.pcolormesh(X,Y,log(Zm),shading='gouraud')
                if matplotlib.__version__>'1.1':
                    plt.xlabel('|Q| ($\AA^{-1}$)')
                    plt.ylabel('$\hbar\omega$ (meV)')
                else:
                    plt.xlabel('|Q| (inverse Angstroms)')
                    plt.ylabel('E (meV)')
                plt.title(str(friendlynamebase))

                plt.savefig(str(friendlynamebase+"_sqw.png"),bbox_inches='tight')
                changepermissions(friendlynamebase+"_sqw.png")
                DeleteWorkspace(MD)
                DeleteWorkspace(MDH)
                
            if 'mdnxs' in datadictsother['Save']:
                #todo: setgoniometer and setub matrix based upon user input.
                minval,maxval=ConvertToMDHelper(out,'Q3D','Direct', 'AutoSelect')

                outMD = ConvertToMD(InputWorkspace = out, QDimensions = "Q3D", Q3DFrames='AutoSelect', QConversionScales="HKL", dEAnalysisMode="Direct", MinValues=minval , MaxValues=maxval,MaxRecursionDepth='1' )
                #save the Md workspace.
                SaveMD(InputWorkspace=outMD, Filename=friendlynamebase+"_MD.nxs")
                #self.datatext += "Data have been saved as a MD nexus file, FILENAME="+friendlynamebase+"_MD.nxs.\n" 
                #change permissions of the directory and file
                changepermissions(friendlynamebase+"_MD.nxs")
                
            if 'iofe' in datadictsother['Save']:
                #get the q binning.  'QTransferRange'
                #datadictsdgs['EnergyTransferRange'] == '' or datadictsdgs['EnergyTransferRange']==None:
                
                
                if datadictsdgs['QTransferRange'] == '' or datadictsdgs['QTransferRange']==None:
                    [qmin, qmax] = calqrangefromworkspace(datawsname)
                    qstep = (qmax-qmin)/150.0
                    qbinparams = str(qmin)+","+str(qstep)+","+str(qmax)
                    

                else:
                    qbinparams = datadictsdgs['QTransferRange']
                
                #Also get the Q binning if the data were all in a single Q bin
                #********************WORKING HERE *******************
                qfullstep



                if self.qmax == None:
                    [tempqmin, qmax] = calqrangefromworkspace(datawsname)
                else:
                    qmax = self.qmax

                #check that qbining has been set
                if self.qstep == None:
                    #if not set then set it to 150 bins for the full range
                    qstep = (qmax-qmin)/150.0
                else:
                    qstep = self.qstep

                qbinparams = str(qmin)+","+str(qstep)+","+str(qmax)
                qfullstep = (qmax-qmin)*1.01
                qfullstepbinparams = str(qmin)+","+str(qfullstep)+","+str(qmax)

                SofQW3(InputWorkspace=datawsname,OutputWorkspace='SofQWdata',QAxisBinning=qbinparams,Emode="Direct",Efixed=efixed)
                Transpose(InputWorkspace='SofQWdata',OutputWorkspace='SofQWdata')
                wsofe = Rebin2D(InputWorkspace='SofQWdata',Axis1Binning=qfullstepbinparams,Axis2Binning=Erange,UseFractionalArea=True,Transpose=True)
                SaveAscii(Filename=friendlynamebase+"_iofe.dat",InputWorkspace='wsofe')
                self.datatext += "Data have been saved as a iofe.dat file, FILENAME="+friendlynamebase+"_iofe.dat.\n" 
                #change permissions of the directory and file
                changepermissions(friendlynamebase+"_iofe.dat")
 
                


    def CreateFriendlyFilename(self,ws,datadictsother):
        #access the sample log to generate the friendlyname
        friendlyfilename = os.path.abspath(os.curdir)+'/'+datadictsother['FriendlyName']+'/'+datadictsother['FriendlyName']
        if datadictsother['FriendlyNameLogs'] != None:
            #get the handle to the run
            run = ws.run()
            for part in datadictsother['FriendlyNameLogs']:
                if run.hasProperty(part):
                    value = run.getProperty(part).value
                 
                    try:
                        friendlyfilename += "_"+part+"_"+value  #This will only be done if value is a string
                    except:
                        #splitting up the value in case of decimal points.
                        #Typically, decimal points in filenames
                        #cause issues with Horace.
                        value=array(value).mean()

                        roundedvalue = "%.2f" % value
                        valuestringwithoutdot = str(roundedvalue).replace('.', 'p')
                        friendlyfilename += "_"+part+"_"+valuestringwithoutdot

        return friendlyfilename      




        
    # Returns a string that is the path and name of the grouping file.
    # if no grouping file can be created or exists, then it returns an error.    
    #If it is powder it will always overwrite
    #If it is NXM (N and M appropraite integers) it checks if it exists, and will only create it once
    def checkgrouping(self,ws,datadictsdgs,datadictsother):
        cwd = os.path.abspath(os.curdir)
        #Always generate the powder file.  This is in case the powder angle step CHANGES for the different data being reduced.
        if ((datadictsdgs['GroupingFile']== 'powder') and (datadictsother['PowderAngleStep'] != None)):
            #generate the powder grouping file
            GenerateGroupingPowder(InputWorkspace=ws,GroupingFilename=cwd+"/powder.xml",AngleStep=datadictsother['PowderAngleStep'])
            datadictsdgs['GroupingFile'] =  os.path.join(cwd,'powder.xml')
            changepermissions(os.path.join(cwd,'powder.xml'))
            return os.path.join(cwd,'powder.xml')
        elif ((datadictsdgs['GroupingFile']== 'powder') and (datadictsother['PowderAngleStep'] == None)):
            raise ValueError("PowderAngleStep not given for powder grouping file.")
        
        elif os.path.isfile(datadictsdgs['GroupingFile']):
            return datadictsdgs['GroupingFile']
        elif os.path.isfile(os.path.join(cwd,datadictsdgs['GroupingFile'])):
            return os.path.join(cwd,datadictsdgs['GroupingFile'])
        elif os.path.isfile(os.path.join(cwd,datadictsdgs['GroupingFile']+'.xml')):
            return os.path.join(cwd,datadictsdgs['GroupingFile']+'.xml')
        else:
            twoelem = datadictsdgs['GroupingFile'].split('X')
            if len(twoelem) == 2:
                #good, now check that they are the correct value
                try:
                    pixely = int(twoelem[0])
                    pixelx = int(twoelem[1])
                except:
                    pixelx = 0
                    pixely = 0
                if (pixely in [1,2,4,8,16,32,64,128]) and (pixelx in [1,2,4,8]):
                    GenerateGroupingSNSInelastic(AlongTubes=str(pixely),AcrossTubes=str(pixelx),Instrument=datadictsother['Instrument'],Filename=os.path.join(cwd,datadictsdgs['GroupingFile']+'.xml'))
                    changepermissions(os.path.join(cwd.datadictsdgs['GroupingFile']+'.xml'))
                    return os.path.join(cwd,datadictsdgs['GroupingFile']+'.xml')
                else:
                    raise ValueError("Can not generate grouping file with current pixel choice")
            else:
                raise ValueError("Grouping File Not Found or created")
                    
                    

def changepermissions(filename):


    """
    change permissions of the directory and file of filename
    to read, write for everyone
    directory also allows execute.
    """
    parentdir = os.path.dirname(filename)
    try:
        os.chmod(parentdir, stat.S_IRWXG | stat.S_IRWXO | stat.S_IRWXU)
    except:
        print "Not able to change permissions of " + parentdir

    try:
        os.chmod(filename, stat.S_IRUSR | stat.S_IRGRP | stat.S_IROTH |stat.S_IWUSR | stat.S_IWGRP | stat.S_IWOTH )
    except:
       print "Not able to change permissions of " + filename




def definegoniometer(names, offsets, directions, axes, ws):
    """
    Function for defining the goniometer
    Returns the psi value which is written to the nxspe file.
    workspace is a string
    returns [psi, text string of goniometer settings]
    """


    outputstring = ""
    psivalue     = 0


    #Transform all inputs to LISTS.
    if str(type(offsets)) != "<type 'list'>" and str(type(offsets)) != "<type 'NoneType'>":
        offsets = [offsets]

    if str(type(names)) != "<type 'list'>" and str(type(names)) != "<type 'NoneType'>":
        names = [names]

    if str(type(directions)) != "<type 'list'>" and str(type(directions)) != "<type 'NoneType'>":
        directions = [directions]

    if str(type(axes)) != "<type 'list'>" and str(type(axes)) != "<type 'NoneType'>":
        axes = [axes]


    #first case, no motor name given (i.e. None)
    if names == None:
        if offsets == None:
            #no motor name, and no offset, just exit
            outputstring = "No goniometer set.\n"
            return [psivalue, outputstring]
        else:
           #check the number of offsets = N axes = Ndirections
            if (len(offsets)==len(directions) and len(directions)==len(axes) and len(offsets)<7):
                #because no motor names given, we need to make a fake log.
                #for loop over the angles
                #list of 6 empty strings that will be filled in
                anglelist = ["","","","","",""]
                try:
                    for i in range(len(offsets)):
                        AddSampleLog(Workspace=ws,LogName="angle"+str(i),
                                    LogText=str(offsets[i]),LogType="Number Series")
                        anglelist[i] = "angle"+str(i)+","+axes[i]+","+str(directions[i])
                        #print anglelist[i]
                except:
                    raise RuntimeError("Could not find goniometer axis(axes)")
                SetGoniometer(Workspace=ws,Axis0=anglelist[0],Axis1=anglelist[1],Axis2=anglelist[2],
                              Axis3=anglelist[3],Axis4=anglelist[4],Axis5=anglelist[5])
                outputstring = "The following axes have been set:\n"

                for i in range(len(offsets)):
                    tempstr = "CCW"
                    if directions[i] == -1:
                        tempstr = "CW"
                    outputstring += "   Axis"+str(i)+" along the "+ axes[i] +" direction, "+tempstr+", rotation angle="+str(offsets[i])+"\n"

                psivalue = offsets[0]
                return [psivalue, outputstring]
            else:
                raise ValueError("Number of angle offsets, directions and axes do not match.")

    #other big case, Motor name is given
    else:
        #Are there any offsets listed
        if offsets==None:
            #create offsets = 0 for ALL motor names listed
            offsets = zeros(len(names))
        #check if the noffsets = Naxes = Ndirectiosn = Nnames
        if (len(offsets)==len(directions) and len(directions)==len(axes) and len(axes)==len(names) and len(offsets)<7):
            #everything is ready
            anglelist = ["","","","","",""]
            anglevalues = []
            try:
                for i in range(len(offsets)):
                    #get the correct log from the workspace
                    angle = mean(ws.run().get(names[i]).value) + offsets[i]
                    anglevalues.append(angle)
                    AddSampleLog(Workspace=ws,LogName="angle"+str(i),
                                LogText=str(angle),LogType="Number Series")
                    anglelist[i] = "angle"+str(i)+","+axes[i]+","+str(directions[i])
            except:
                raise RuntimeError("Could not find goniometer axis "+names[i])
            SetGoniometer(Workspace=ws,Axis0=anglelist[0],Axis1=anglelist[1],Axis2=anglelist[2],
                              Axis3=anglelist[3],Axis4=anglelist[4],Axis5=anglelist[5])
            outputstring = "The following axes have been set:\n"

            for i in range(len(offsets)):
                tempstr = "CCW"
                if directions[i] == -1:
                    tempstr = "CW"
                outputstring += "   Axis"+str(i)+" along the "+ axes[i] +" direction, "+tempstr+", rotation angle="+names[i]+"+" +str(offsets[i])+"="+str(anglevalues[i])+"\n"

            return [anglevalues[0],outputstring]
        else:
            raise ValueError("Number of angle names, offsets, directions and axes do not match.")





      
if __name__ == "__main__":
    #check number of arguments
    if (len(sys.argv) != 2): 
        print "reduction code requires a datatext file"
        sys.exit()
    if not(os.path.isfile(sys.argv[1])):
        print "data text file ", sys.argv[1], " not found"
        sys.exit()
    dgsreduction(XMLfile=sys.argv[1])       
