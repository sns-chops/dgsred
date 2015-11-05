'''
Created on Feb 20, 2012
Modified for use with dgsreductionmantid v6.0 October, 2013
    v6.0 uses mantid dgsreduction alg.

@author: andrei
'''
from lxml import etree
#from string import lower
from copy import deepcopy,copy
import os
import imp
from mantid.kernel import logger
from ARLibrary import _parseBTPlist

class XMLparser(object):
    '''
    classdocs
    '''

    def __init__(self,filename):
        '''
        Constructor
        '''

        #KeyWords are listed here in alphabetical order
        self.OtherAllowedKeyWords = ['DataPath','FilterBadPulses','FilterMax','FilterMin','FilterNames','FriendlyName',
        'FriendlyNameLogs','GoniometerMotor','GoniometerMotorAxis','GoniometerMotorDirection','GoniometerMotorOffset',
        'Instrument','LatticeParams','LogValue','LogValueMax','LogValueMin', 'LogValueStep','Mask',
        'NormalizedCalibration','PowderAngleMax','PowderAngleMin','PowderAngleStep','QTransferRange','Runs', 'Save',
        'ScanType','UB','UVector','VanPath','VanRuns','VVector',]                             
      
        #KeyWords are listed here in alphabetical order
        self.DGSAllowedKeyWords = ['AbsUnitsMaximumEnergy','AbsUnitsMinimumEnergy','CorrectKiKf','DetVanIntRangeHigh',
        'DetVanIntRangeLow','DetVanIntRangeUnits','EnergyTransferRange','ErrorBarCriterion','GroupingFile',
        'HardMaskFile','HighCounts','HighOutlier','IncidentBeamNormalization','IncidentEnergyGuess','LowCounts',
        'LowOutlier','MedianTestCorrectForSolidAngle','MedianTestHigh','MedianTestLevelsUp','MedianTestLow',
        'MonitorIntRangeHigh','MonitorIntRangeLow','Monitor1SpecId','Monitor2SpecId','SaveProcDetVanFilename',
        'SofPhiEIsDistribution','TimeIndepBackgroundSub','TibTofRangeEnd','TibTofRangeStart','TimeZeroGuess',
        'UseIncidentEnergyGuess']


        if filename!=None:
            doc=etree.parse(filename)
            root=doc.getroot()
            if root.tag not in ['dgsreduction']:
                raise RuntimeError("This is not a dgsreduction xml file")
            self.root=root
            defaults=root.find('defaults')
            if defaults==None:
                raise RuntimeError("The defaults section is missing")
            temp={}
            self.recursiveGetAttributes(defaults,temp)
            #First verify that there is an instrument in the tempdefaultdict.
            if not temp.has_key('Instrument'):
                raise RuntimeError("The Instrument is not defined")
            if temp['Instrument'] not in ['ARCS','SEQUOIA','HYSPEC','CNCS']:
                raise ValueError("Unknown Instrument "+str(temp['Instrument']))    
            else:
                #if there is an instrument value, then try to load the
                #defaults for that instrument.
                self.loadinstrumentdict(temp['Instrument'])
            self.defaultdict.update(temp)

            self.calibdict=deepcopy(self.defaultdict)
            calibration=root.find('calibration')
            self.recursiveGetAttributes(calibration,self.calibdict)                      
            self.datadicts=[]
            self.datadictsdgs=[]
            self.datadictsother=[]
            listruns = root.findall('scan')
        
            for i in range(len(listruns)):
                #for some reason deepcopy did not work here.
                #issued a seg fault
                tempi = copy(self.calibdict)
                self.recursiveGetAttributes(listruns[i],tempi)
                if tempi['Runs']!=[]:
                    self.datadicts.append(tempi)
                    #now split up the datadicts to be TWO arrays.  datadictsdgs and datadictsothers
                    tempidgs = {}
                    tempiother = {}
                    for k,v in tempi.items():
                        if k in self.DGSAllowedKeyWords:
                            tempidgs[k] = v
                        else:
                            if k != 'Mask':
                                tempiother[k] = v
                            else:
                                maskdictother=[]
                                for maskdicti in v:
                                    if maskdicti.has_key('algorithm'):
                                        maskdictother.append(maskdicti)
                                    else:
                                        for u,w in maskdicti.items():
                                            tempidgs[u]=w
                                if len(maskdictother) != 0:
                                    tempiother['Mask'] = maskdictother
                                                        
                    self.datadictsdgs.append(tempidgs)
                    self.datadictsother.append(tempiother)
                                  

                
    def recursiveGetAttributes(self,node,dictionary):
        if node.tag!=etree.Comment:
            if node.tag=='Mask':
                if not dictionary.has_key('Mask'):
                    dictionary['Mask']=[]
                if node.attrib!={}:        
                    dictionary['Mask'].append(node.attrib)
            else:
                for k in node.attrib.keys():
                    if (k not in self.DGSAllowedKeyWords) and (k not in self.OtherAllowedKeyWords):
                        raise ValueError("Unknown Keyword "+k+" in the XML file")
                    dictionary[k]=self.parseto(k,node.attrib[k])
                if (node.text!=None) and (len(node.text.strip())!=0):
                    if (node.tag not in self.DGSAllowedKeyWords) and (node.tag not in self.OtherAllowedKeyWords):
                        raise ValueError("Unknown Keyword "+node.tag+" in the XML file")
                    dictionary[node.tag]=self.parseto(node.tag,node.text)
                for c in node.iterchildren():
                    self.recursiveGetAttributes(c,dictionary)

    def parseto(self,tag,value):
        if tag == 'ScanType':
            tstring = value[0:2].lower()
            if tstring == 'si':
                return 'single'
            elif tstring == 'st':
                return 'step'
            elif tstring == 'sw':
                return 'sweep'
            else:
                raise ValueError("scantype not understood (single, step, sweep or vannorm).")
                
        #We choose not to convert to floats for the variables, but keep as strings.
        #if tag in ['t0','efixed','emin','emax','ebin','qmin','qmax','qstep','tibgstart','tibgstop','vanmin','vanmax',
        #           'logvaluemin','logvaluemax','logvaluestep','powderanglestep','powderanglemin','powderanglemax',
        #           'vanmine','vanmaxe','vanminangle','vanmaxangle','vanmass','samplemass','samplermm','scalefactor']:
        #    return float(value)

        if tag in ['Runs','VanRuns']:
            return _parseBTPlist(value)

        if tag in ['UseIncidentEnergyGuess','CorrectKiKf','TimeIndepBackgroundSub','NormalizedCalibration','FilterBadPulses','MedianTestCorrectForSolidAngle']:
            return (value.lower() in ("yes", "true", "t", "1", "tru", "tr", "y"))

        if tag=='DetVanIntRangeUnits':
            #replace all whitespace with nothing
            value.replace(' ', '')
            #this is a list of strings of the supported filetypes
            unittypesexist = ['DeltaE', 'DeltaE_inWavenumber', 'Energy', 'Energy_inWavenumber',
                              'Momentum', 'MomentumTransfer', 'QSquared', 'TOF', 'Wavelength',
                              'dspacing']
            #get a string all lowercase.
            unitlower = []
            for u in unittypesexist:
                unitlower.append(u.lower())
            #Find the matching value of the unit list.
            i = unitlower.index(value.lower())
            return unittypesexist[i]

        if tag=='Save':
            #replace all whitespace with commas
            value.replace(' ', ',')
            #split up parts by comma delimited values
            parts = value.split(',')
            #this is a list of strings of the supported filetypes
            filetypesexist = ['phx', 'spe', 'nxspe', 'par', 'jpg', 'nxs', 'mdnxs', 'iofq','iofe',
                            'iofphiecolumn','iofphiearray','iofqecolumn','iofqearray','summary',
                            'sqw','vannorm']
            #loop over the elements and fill up a list of the
            #filetypes to return.
            filestoprocess = []
            for p in parts:
                p = p.strip(' "')
                if p.lower() in filetypesexist:
                    filestoprocess.append(p.lower())
                else:
                    if len(p) > 0:
                        logger.warning("*****FILETYPE "+p+" DOES NOT EXIST. This file will be ignored.*****")
            return filestoprocess

        if tag=='FriendlyNameLogs':
            #convert to an arry of strings, split by the commas#
            parts = value.replace(' ', '').split(',')
            return parts

        if tag=='FilterNames':
            parts = value.replace(' ', '').split(',')
            return parts

        if tag=='FilterMin':
            parts = value.replace(' ', '').split(',')
            return parts

        if tag=='FilterMax':
            parts = value.replace(' ', '').split(',')
            return parts

        if tag in ['GoniometerMotor','GoniometerMotororAxis','GoniometerMotorOffset','GoniometerMotorDirection']:
            try:
                return eval(value)
            except:
                return value            

        return value



    def loadinstrumentdict(self,instrumentname):
        """
        instrumentname is a string
        returns a dictionary of instrument defaults
        """


        modulename = instrumentname.lower()+'default'

        #We have the instrument name.
        #Try to import the default file for that instrument
        try:
            m = imp.load_source(modulename,os.path.abspath(os.curdir)+"/"+modulename+'.py')

        except:
            try:
                m = imp.load_source(modulename,os.path.dirname(os.path.realpath(__file__))+"/"+modulename+'.py')
            except:
                raise RuntimeError("Could not find instrument definition module "+os.path.abspath(os.curdir)+"/"+modulename+'.py')

        tempdict = m.instrumentparameters()
        for k in tempdict.keys():
            if (k not in self.DGSAllowedKeyWords) and (k not in self.OtherAllowedKeyWords):
                raise ValueError("Unknown Keyword "+k+" in the instrument default file")

        self.defaultdict=tempdict




if __name__ == "__main__":
    filename='newstylexml.xml'
    parsed=XMLparser(filename)
    
    print "-----------default--------------"
    print parsed.defaultdict 
    print "-----------calibration--------------"
    print parsed.calibdict
    print "-----------data--------------"
    print parsed.datadicts
    print "----------dictdgs-----------"
    print parsed.datadictsdgs
    print "------------dictother--------"
    print parsed.datadictsother
