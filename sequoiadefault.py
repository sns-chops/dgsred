#this is a dictionary of all possible things which have to be parsed and dealt with
#that belong to the instrument

#ALL KEYWORDS FOR THE DEFAULT DICTIONARY SHOULD BE lowercase.

# RICARDO:
#
# DUMMY FILE
# I did this for testing in FERMI
#
#

def instrumentparameters():

    globaldict = dict()


    #instrument items
    globaldict['Instrument']='SEQUOIA'
    globaldict['FilterBadPulses'] = True

    #calibration and mask items
    globaldict['VanRuns'] = None
    globaldict['DetVanIntRangeUnits']   = 'TOF'
    globaldict['DetVanIntRangeLow']  = None
    globaldict['DetVanIntRangeHigh']  = None
    globaldict['SaveProcDetVanFilename'] = 'vanadium.nx5'
    globaldict['HardMaskFile'] = None
    globaldict['Mask'] = [] #the mask parameters are stored as a list of dictionaries. (check this)
    globaldict['NormalizedCalibration'] = False
    globaldict['VanPath'] = ''


    #data items
    globaldict['DataPath'] = ''
    globaldict['Runs'] = []  #returns a list of ints
    globaldict['IncidentEnergyGuess']= None #number double
    globaldict['TimeZeroGuess'] = 0.0 #number double
    globaldict['UseIncidentEnergyGuess'] = True
    globaldict['Monitor1SpecId'] = 1
    globaldict['Monitor2SpecId'] = 2
    globaldict['EnergyTransferRange']=None #number double
    globaldict['QTransferRange']=None #number double
    globaldict['CorrectKiKf']=True #can only use boolean here, not 1 or 0.
    globaldict['TimeIndepBackgroundSub']=False
    globaldict['TibTofRangeStart']=None
    globaldict['TibTofRangeEnd']= None
    globaldict['GroupingFile']=None
    globaldict['PowderAngleStep'] = 0.5
    globaldict['GoniometerMotor']=None  #script takes care of default situations.
    globaldict['GoniometerMotorOffset']=None
    globaldict['GoniometerMotorAxis']=["0,1,0"]
    globaldict['GoniometerMotorDirection']=[1]
    globaldict['Save']=[]
    globaldict['FriendlyName']='ARCS'
    globaldict['FriendlyNameLogs'] = ['run_number']  #will get the run_number from the logs, first run number only. (must be in square brackets here)
    globaldict['FilterNames']=None
    globaldict['FilterMin']=None
    globaldict['FilterMax']=None
    return globaldict

