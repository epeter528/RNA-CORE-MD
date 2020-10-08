# import chimera
from ctypes import *
from ctypes.util import find_library


import numpy
from numpy.ctypeslib import ndpointer

from pyMMB import *

## Types of polymers. Indices match MMB types enum
PolyTypes = ('RNA','protein','DNA','Unassigned')

##################################################################
MMB.init.argtypes = [c_char_p, c_char_p]
MMB.printBiopolymerCoordinates.argtypes = [c_char_p, c_char_p]
MMB.printAllSettings.argtypes           = [c_char_p]

def initMMB(parametersFileName=''):
    call('init', parametersFileName)

## Initialize the particles coordinates and bonds.
def initializeMolecules():
    return call('initializeMolecules')
MMB.initializeMolecules.argtypes        = [c_char_p]

## Write a pdb file of the current state of the polymers.
#  @param outputName the output pdb file.
def writePdb(filename):
    return call('writePdb', filename)
MMB.writePdb.argtypes                   = [c_char_p, c_char_p]

## Write a pdb file of the default state of the polymers.
#  @param outputName the output pdb file.
def writeDefaultPdb(filename):
    return call('writeDefaultPdb', filename)
MMB.writeDefaultPdb.argtypes                   = [c_char_p, c_char_p]

## Postinitialization, to be call after all settings, polymers, 
#  forces and constraints are set.
def postInitialize():
    return call('postInitialize')
MMB.postInitialize.argtypes             = [c_char_p]

## Initialize the dynamics.
#  The last frame will be write in a file called last.1.pdb, 
#  the trajectory in trajectory.1.pdb and current frame in frame.pdb.
def initDynamics():
    return call('initDynamics')
MMB.initDynamics.argtypes               = [c_char_p]

## Run one step of the dynamics as defined by current parameters.
#  @return the remaining number of steps.
def runOneStep():
    return call('runOneStep')
MMB.runOneStep.argtypes                 = [c_char_p]
MMB.runOneStep.restype                 = c_int

## Run all steps of the dynamics as defined by current parameters.
def runAllSteps():
    return call('runAllSteps')
MMB.runAllSteps.argtypes                = [c_char_p]

## End dynamics.
#  Write the last frame and delete the integrator
def postDynamics():
    return call('postDynamics')
MMB.postDynamics.argtypes               = [c_char_p]


## Make all necessary initialization and run all steps of the dynamics as defined by current parameters.
#  Buggy...
#  Do not call another initialization if you use runDynamics.
def runDynamics():
    return call('runDynamics')
MMB.runDynamics.argtypes                = [c_char_p]


## Initialize the coordinates of a polymer identified by its chain ID.
#  Use coordinates frome pdbFileName or generate default coordinates if pdbFileName is empty.
#  The chainID and number of atoms of the polymer must match one of the chain in the pdb file.
#  @param chainID the chain id of the polymer to initialize.
#  @param pdbFileName the pdb file to extract the coordinates from
def initializePolymer(chainID, filename=""):
    return call('initializePolymer', chainID, filename)
MMB.initializePolymer.argtypes          = [c_char_p, c_char_p, c_char_p]

def initializeMoleculesAndBonds():
    return call('initializeMoleculesAndBonds')
MMB.initializeMoleculesAndBonds.argtypes           = [c_char_p]

## Initialize bonds, the multibody tree and state.
#  Necessary to get a correct structure.
def initializeBodies():
    return call('initializeBodies')
MMB.initializeBodies.argtypes           = [c_char_p]



## Find residues within the radius (nm) of resID
#  @return a list of ("chainID", resID (int))
def getResiduesWithin(chainID, resID, radius):
    resWithin = call('getResiduesWithin', chainID, resID, radius)
    sets = []
    for r in resWithin.split(",")[:-1]:
        r = r.strip().split()
        sets.append( ( r[0], int(r[1]) ) )
    return sets
MMB.getResiduesWithin.argtypes          = [c_char_p, c_int, c_double, c_char_p]
MMB.getResiduesWithin.restype          = c_char_p



####################################################################
## Store sequences information from MMB BiopolymerClass
class BioPolymerSeq(Structure):
    _fields_ = [('chainID', c_char_p),
                ('sequence', c_char_p),
                ('firstResNum', c_int),
                ('polyType', c_int),
                ('pdbFileName', c_char_p),
                ('loadFromPdb', c_bool),
                ('activePhysics', c_bool)]

    def getTypeName(self):
        return PolyTypes[self.polyType]

    def mmbUpdate(self):
        call("updatePolymer", byref(self))

    def __str__(self):
        typeName = PolyTypes[self.polyType]
        return typeName + " " + self.chainID + " " + str(self.firstResNum) + " " + self.sequence


BioPolymerSeq_ptr = POINTER(BioPolymerSeq)

MMB.deletePolymer.argtypes = [c_char_p, c_char_p]
def deletePolymer(chainID):
    call("deletePolymer", chainID)

MMB.updatePolymer.argtypes = [c_void_p, c_char_p]
# def updatePolymer(chainID, seq, pdbFile, loadFromPdb):
#     call("updatePolymer", chainID, seq, pdbFile, loadFromPdb)

## Get sequences currently in MMB
#  @param seqs array of sequences to fill for external use
#  @return a list of the sequences in MMB as BioPolymerSeq objects
def getSequences():
    seqs = POINTER(BioPolymerSeq)()
    nbSeqs = call('getSequences',byref(seqs))
    return seqs[0:nbSeqs]
MMB.getSequences.argtypes = [c_void_p, c_char_p]
MMB.getSequences.restype = c_int

## Load the sequences from the file to MMB.
#  initializePolymer will use the name of this file if needed
#  @param pdbFileName the name of the pdb file to extract the sequences from
#  @param seqs array of sequences to fill for external use
#  @return The number of sequences extracted
def loadSequencesFromPdb(filename):
    """
    Do not import coordinates. Just the sequences.
    """
    return call('loadSequencesFromPdb',filename)
MMB.loadSequencesFromPdb.argtypes       = [c_char_p, c_char_p]


def loadSequencesFromPdbContent(pdbStream):
    """
    Do not import coordinates. Just the sequences.
    """
    return call('loadSequencesFromPdb',pdbStream.getvalue())
MMB.loadSequencesFromPdbContent.argtypes       = [c_char_p, c_char_p]

####################################################################
class MobilizerStretch_wrapper(Structure):
    _fields_ = [('mmbID', c_int),
                ('chainID', c_char_p),
                ('mobility', c_char_p),
                ('resStart', c_int),
                ('resEnd', c_int)]

    def __str__(self):
        resStart = str(self.resStart)
        resEnd = str(self.resEnd)
        chain = str(self.chainID)

        if chain == "All":
            resStart = ""
            resEnd = ""
            chain = ""
        elif resStart == "0":
            resStart = ""
            resEnd = ""

        return "mobilizer %s %s %s %s" % (self.mobility,
                                          chain,
                                          resStart, resEnd)

    def mmbAdd(self):
        cmd(str(self))

    ## Update a Mobilizer
    def mmbUpdate(self):
        call("updateMobilizerStretch", self.mmbID, self.chainID, 
                                              self.resStart, self.resEnd, 
                                              self.mobility)

    ## Delete a Mobilizer stretch identified by its mmbID
    #  @param mmbID the id in MMB's vector
    def mmbDelete(self):
        return call('deleteMobilizerStretch', self.mmbID)

MobilizerStretch_ptr = POINTER(MobilizerStretch_wrapper)
# Mobilizers                                           
MMB.updateMobilizerStretch.argtypes     = [c_int, c_char_p, c_int, c_int, \
                                                  c_char_p, c_char_p]
MMB.getMobilizerStretches.argtypes      = [c_void_p, c_char_p]
MMB.getMobilizerStretches.restype      = c_int
MMB.deleteMobilizerStretch.argtypes     = [c_int, c_char_p]

## Return a list of Mobilizer stretch struct
#  @return a list of MobilizerStretch_wrapper
def getMobilizerStretches():
    mobilizers = POINTER(MobilizerStretch_wrapper)()
    nbMob = call('getMobilizerStretches', byref(mobilizers))
    return mobilizers[0:nbMob]

####################################################################
class MobilizerWithin_wrapper(Structure):
    _fields_ = [('mmbID', c_int),
                ('chainID', c_char_p),
                ('mobility', c_char_p),
                ('resID', c_int),
                ('radius', c_double)]

    ## Return the MMB command as a string
    def __str__(self):
        return "applyMobilizersWithin %s %.1f %s %i" % (self.mobility, 
                                                        self.radius, 
                                                        self.chainID, 
                                                        self.resID)

    def mmbAdd(self):
        cmd(str(self))

    ## Update a MobilizerWithin
    def mmbUpdate(self):
        return call("updateMobilizerWithin", self.mmbID, self.chainID, 
                                              self.resID, self.radius, 
                                              self.mobility)

    ## Delete a Mobilizer stretch identified by its mmbID
    #  @param mmbID the id in MMB's vector
    def mmbDelete(self):
        return call('deleteMobilizerWithin', self.mmbID)

MobilizerWithin_ptr = POINTER(MobilizerWithin_wrapper)
MMB.updateMobilizerWithin.argtypes      = [c_int, c_char_p, c_int, c_double, \
                                                  c_char_p, c_char_p]
MMB.getMobilizersWithin.argtypes        = [c_void_p, c_char_p]
MMB.getMobilizersWithin.restype         = c_int
MMB.deleteMobilizerWithin.argtypes      = [c_int, c_char_p]

## Return a list of MobilizerWithin struct
#  @return a list of MobilizerWithin_wrapper
def getMobilizersWithin():
    mobilizers = POINTER(MobilizerWithin_wrapper)()
    nbMob = call('getMobilizersWithin', byref(mobilizers))
    return mobilizers[0:nbMob]

####################################################################
class Constraint_wrapper(Structure):
    _fields_ = [('chain1', c_char_p),
                ('chain2', c_char_p),
                ('atom1', c_char_p),
                ('atom2', c_char_p),
                ('mmbID', c_int),
                ('res1', c_int),
                ('res2', c_int)]

    def mmbAdd(self):
        if self.mmbID > 0:
            raise MMBError("Constraint with id %i already exists in MMB." % self.mmbID)
        # if self.chain1=="":
        #     return cmd("rootMobilizer Weld")
        # if self.res1=="":
        #     return cmd("rootMobilizer " + self.chain1 + " Weld")
        com = "constraint " + self.chain1 + " " + str(self.res1) + " " + self.atom1
        if self.chain2 == "Ground":
            return cmd("constrainToGround " + self.chain1 + " " + str(self.res1))
        return cmd(com + " Weld " + self.chain2 + " " + str(self.res2) + " " + self.atom2)


    ## Update a Constraint
    def mmbUpdate(self):
        # pyMMB.updateConstraint(self.mmbID, 
        #                        self.chain1, res1, self.atom1, 
        #                        self.chain2, res2, self.atom2)
        call("updateConstraint", self.mmbID, 
                                 self.chain1, self.res1, self.atom1,
                                 self.chain2, self.res2, self.atom2)

    ## Delete a Constraint identified by its mmbID
    #  @param mmbID the id in MMB's vector
    def mmbDelete(self):
        return call('deleteConstraint', self.mmbID)


Constraint_ptr = POINTER(Constraint_wrapper)
# Constraints
MMB.getConstraints.argtypes             = [c_void_p, c_char_p]
MMB.deleteConstraint.argtypes           = [c_int, c_char_p]
MMB.updateConstraint.argtypes           = [c_int, c_char_p, c_int, c_char_p,\
                                           c_char_p, c_int, c_char_p, c_char_p]
MMB.getConstraints.restype              = c_int
## Return a list of Constraints_wrapper struct
#  @return a list of Constraints_wrapper
def getConstraints():
    constraints = POINTER(Constraint_wrapper)()
    nbConst = call('getConstraints', byref(constraints))
    return constraints[0:nbConst]

####################################################################
class ContactStretch_wrapper(Structure):
    _fields_ = [('mmbID', c_int),
                ('chainID', c_char_p),
                ('contactScheme', c_char_p),
                ('resStart', c_int),
                ('resEnd', c_int)]

    def __str__(self):
        resStart = str(self.resStart)
        resEnd = str(self.resEnd)

        if self.chainID == "All":
            self.chainID = ""
            resStart = ""
            resEnd = ""
        elif self.resStart == 0:
            resStart = ""
            resEnd = ""

        return "contact %s %s %s %s" % (self.contactScheme,
                                          self.chainID,
                                          resStart, resEnd)

    def mmbAdd(self):
        cmd(str(self))

    ## Update a Contact
    def mmbUpdate(self):
        call("updateContactStretch", self.mmbID, self.chainID, 
                                              self.resStart, self.resEnd, 
                                              self.contactScheme)

    ## Delete a Contact stretch identified by its mmbID
    #  @param mmbID the id in MMB's vector
    def mmbDelete(self):
        return call('deleteContactStretch', self.mmbID)

ContactStretch_ptr = POINTER(ContactStretch_wrapper)
# Contacts                                           
MMB.updateContactStretch.argtypes     = [c_int, c_char_p, c_int, c_int, \
                                                  c_char_p, c_char_p]
MMB.getContactStretches.argtypes      = [c_void_p, c_char_p]
MMB.getContactStretches.restype      = c_int
MMB.deleteContactStretch.argtypes     = [c_int, c_char_p]

## Return a list of Contact stretch struct
#  @return a list of ContactStretch_wrapper
def getContactStretches():
    contacts = POINTER(ContactStretch_wrapper)()
    nbMob = call('getContactStretches', byref(contacts))
    return contacts[0:nbMob]

####################################################################
class ContactWithin_wrapper(Structure):
    _fields_ = [('mmbID', c_int),
                ('chainID', c_char_p),
                ('contactScheme', c_char_p),
                ('resID', c_int),
                ('radius', c_double)]

    ## Return the MMB command as a string
    def __str__(self):
        return "applyContactsWithin %.1f %s %s %i" % (  self.radius, 
                                                        self.contactScheme,
                                                        self.chainID, 
                                                        self.resID)

    def mmbAdd(self):
        cmd(str(self))

    ## Update a ContactWithin
    def mmbUpdate(self):
        return call("updateContactWithin", self.mmbID, self.chainID, 
                                              self.resID, self.radius, 
                                              self.contactScheme)

    ## Delete a Contact stretch identified by its mmbID
    #  @param mmbID the id in MMB's vector
    def mmbDelete(self):
        return call('deleteContactWithin', self.mmbID)

ContactWithin_ptr = POINTER(ContactWithin_wrapper)
MMB.updateContactWithin.argtypes      = [c_int, c_char_p, c_int, c_double, \
                                                  c_char_p, c_char_p]
MMB.getContactsWithin.argtypes        = [c_void_p, c_char_p]
MMB.getContactsWithin.restype         = c_int
MMB.deleteContactWithin.argtypes      = [c_int, c_char_p]

## Return a list of ContactWithin struct
#  @return a list of ContactWithin_wrapper
def getContactsWithin():
    contacts = POINTER(ContactWithin_wrapper)()
    nbMob = call('getContactsWithin', byref(contacts))
    return contacts[0:nbMob]

####################################################################
class AllResiduesWithin_wrapper(Structure):
    _fields_ = [ 
            ('mmbID', c_int),
            ('chain', c_char_p),
            ('residue', c_int),
            ('radius', c_double)
            ]

    def mmbAdd(self):
        if self.mmbID > 0:
            raise MMBError("AllResiduesWithin with id %i already exists in MMB." % self.mmbID)
        cmd("includeAllResiduesWithin %.1f %s %i" % (self.radius, self.chain, self.residue))

    def mmbUpdate(self):
        call("updateAllResiduesWithin",self.mmbID, self.chain, self.residue, self.radius)

    def mmbDelete(self):
        call("deleteAllResiduesWithin", self.mmbID)

AllResiduesWithin_ptr = POINTER(AllResiduesWithin_wrapper)
MMB.deleteAllResiduesWithin.argtypes = [c_int, c_char_p]
MMB.updateAllResiduesWithin.argtypes = [c_int, c_char_p, c_int, c_double, c_char_p]

MMB.getAllResiduesWithins.argtypes = [c_void_p, c_char_p]
MMB.getAllResiduesWithins.restype = c_int
def getAllResiduesWithins():
    objs = AllResiduesWithin_ptr()
    nbObjs = call('getAllResiduesWithins', byref(objs))
    return objs[0:nbObjs]


####################################################################
class IncludeAllNonBondAtomsInResidue_wrapper(Structure):
    _fields_ = [ 
            ('mmbID', c_int),
            ('chain', c_char_p),
            ('residue', c_int)
            ]

    def mmbAdd(self):
        if self.mmbID > 0:
            raise MMBError("IncludeAllNonBondAtomsInResidue with id %i already exists in MMB." % self.mmbID)
        cmd("includeAllNonBondAtomsInResidues %s %i %i" % (self.chain, self.residue, self.resEnd))

    def mmbUpdate(self):
        call("updateIncludeAllNonBondAtomsInResidue",self.mmbID, self.chain, self.residue)

    def mmbDelete(self):
        print self.mmbID
        call("deleteIncludeAllNonBondAtomsInResidue", self.mmbID)


IncludeAllNonBondAtomsInResidue_ptr = POINTER(IncludeAllNonBondAtomsInResidue_wrapper)
MMB.deleteIncludeAllNonBondAtomsInResidue.argtypes = [c_int, c_char_p]
MMB.updateIncludeAllNonBondAtomsInResidue.argtypes = [c_int, c_char_p, c_int, c_char_p]

MMB.getIncludeAllNonBondAtomsInResidues.argtypes = [c_void_p, c_char_p]
MMB.getIncludeAllNonBondAtomsInResidues.restype = c_int
def getIncludeAllNonBondAtomsInResidues():
    objs = IncludeAllNonBondAtomsInResidue_ptr()
    nbObjs = call('getIncludeAllNonBondAtomsInResidues', byref(objs))
    return objs[0:nbObjs]

####################################################################
# Base interactions

class BaseInteraction_wrapper(Structure):
    _fields_ = [ 
            ('mmbID', c_int),
            ('edge1', c_char_p),
            ('edge2', c_char_p),
            ('poly1', c_char_p),
            ('poly2', c_char_p),
            ('res1', c_int),
            ('res2', c_int),
            ('bondOrient', c_char_p),
            #('BasePairIsTwoTransformForce', ////String: unknown ctype),
            #('BasePairPriority', ////int: unknown ctype),
            #('BasePairTemporary', ////int: unknown ctype),
            #('rotationCorrection1', //Rotation: unknown ctype),
            #('rotationCorrection2', //Rotation: unknown ctype),
            #('translationCorrection1', //Vec3: unknown ctype),
            #('translationCorrection2', //Vec3: unknown ctype),
            ('basePairSatisfied', c_bool),
            ('leontisWesthofBondRowIndex', c_int)
            ]

    ## Return an MMB baseInteraction command in a string
    #  @return a string like: baseInteraction A 1 WatsonCrick A 10 WatsonCrick Trans 
    def __str__(self):
        return "baseInteraction %s %i %s %s %i %s %s" % \
                    (self.poly1, self.res1, self.edge1, 
                     self.poly2, self.res2, self.edge2, self.bondOrient
                    )

    def mmbAdd(self):
        return cmd(str(self))

    ## Update the base interaction identified by the MMBid with the content of the base pair in argument.
    def mmbUpdate(self):
        return call('updateBasePair', self.mmbID, 
                                    self.poly1, self.res1, self.edge1, 
                                    self.poly2, self.res2, self.edge2, self.bondOrient)

    ## Delete a base pair identified by its mmbID
    def mmbDelete(self):
        return call('deleteBasePair', self.mmbID)

BaseInteraction_ptr = POINTER(BaseInteraction_wrapper)
MMB.getBaseInteractions.argtypes = [c_void_p, c_char_p]
MMB.getBaseInteractions.restype = c_int
def getBaseInteractions():
    objs = BaseInteraction_ptr()
    nbObjs = call('getBaseInteractions', byref(objs))
    return objs[0:nbObjs]


MMB.getBaseInteractionsStrings.argtypes = [c_char_p]
MMB.getBaseInteractionsStrings.restype = c_char_p
MMB.deleteBasePair.argtypes             = [c_int, c_char_p]
MMB.updateBasePair.argtypes             = [c_int, c_char_p, c_int, c_char_p, \
                                                  c_char_p, c_int, c_char_p, \
                                           c_char_p, c_char_p]
## Return all the base interactions currently in MMB as an array of strings.
#  MMBid ch1 resN1 edge1 ch2 resN2 edge2 orientation
#  @return list of MMB input strings
def getBaseInteractionsStrings():
    """
    String format:
    MMBid ch1 resN1 edge1 ch2 resN2 edge2 orientation 
    """
    interactions = call('getBaseInteractionsStrings');
    return interactions

## Update the base pair MMB_basePair_ID with
#  information from basePairString
#  @param MMB_basePair_ID MMB id of the base pair
#  @param basePairString baseInteraction ch1(str) res1(int) Edge1(str) ch2(str) res2(int) Edge2(str) Orientation(str)
def updateBasePair(MMB_basePair_ID, basePairString):
    p = basePairString.split()
    return call('updateBasePair', MMB_basePair_ID, 
                                   p[1], int(p[2]), p[3], 
                                   p[4], int(p[5]), p[6], 
                                   p[7])

## Remove the base interaction 'MMB_basePair_ID' from MMB.
#  Warning: shift the MMB id of the subsequent interactions from -1.
#  @param id MMB id of the interaction to delete
def deleteBasePair(MMB_basePair_ID):
    return call('deleteBasePair', MMB_basePair_ID)

MMB.getNumSatisfiedBasePairs.argtypes = [c_char_p]
MMB.getNumSatisfiedBasePairs.restype = c_int
def getNumSatisfiedBasePairs():
    return call('getNumSatisfiedBasePairs')

MMB.getNumUnSatisfiedBasePairs.argtypes = [c_char_p]
MMB.getNumUnSatisfiedBasePairs.restype = c_int
def getNumUnSatisfiedBasePairs():
    return call('getNumUnSatisfiedBasePairs')

## Remove all polymers from MMB.
def clearPolymers():
    return call('clearSequences')
MMB.clearSequences.argtypes              = [c_char_p]

## Remove all forces and constraints from MMB.
def clearForcesAndConstraints():
    return call('clearForcesAndConstraints')
MMB.clearForcesAndConstraints.argtypes  = [c_char_p]

## Return the total number of atoms of MMB's current polymers.
#  @return total number of atoms.
def getSystemNumAtoms():
    return call('getSystemNumAtoms')
MMB.getSystemNumAtoms.argtypes          = [c_char_p]
MMB.getSystemNumAtoms.restype          = c_int

## Return a linear float array of the coordinates of all polymers.
#  @return a numpy array [x0,y0,z0, x1,y1,z1, ..., xn,yn,zn]
def getSystemCoordinates():
    """
    Return a numpy array containing coordinates of all atoms.
    [[x0,y0,z0], [x1,y1,z1], ..., [xn,yn,zn]]
    """
    coordArray = numpy.zeros((getSystemNumAtoms(),3), dtype=numpy.float32)
    call('getSystemCoordinates', coordArray)
    return coordArray
MMB.getSystemCoordinates.argtypes       = [ndpointer(c_float), c_char_p]

def matchCoordinatesFromContent(chainID, pdbFileStream):
    return call('matchCoordinatesFromContent', chainID, pdbFileStream.getvalue())
MMB.matchCoordinatesFromContent.argtypes = [c_char_p, c_char_p, c_char_p]

###################################################################
class ThreadingStruct_wrapper(Structure):
    _fields_ = [ 
            ('mmbID', c_int),
            ('chainID1', c_char_p),
            ('residueStart1', c_int),
            ('residueEnd1', c_int),
            ('chainID2', c_char_p),
            ('residueStart2', c_int),
            ('residueEnd2', c_int),
            ('forceConstant', c_double),
            ('backboneOnly', c_bool)
            ]

    def mmbAdd(self):
        if self.mmbID > 0:
            raise MMBError("Threading with id %i already exists in MMB." % self.mmbID)
        call("addThreading", byref(self))

    def mmbUpdate(self):
        call("updateThreading",self.mmbID, byref(self))

    def mmbDelete(self):
        call("deleteThreading", self.mmbID)

    def __str__(self):
        return "threading %s %i %i %s %i %i %f" % (self.chainID1, self.residueStart1, self.residueEnd1,
                                                   self.chainID2, self.residueStart2, self.residueEnd2,
                                                   self.forceConstant)

class GappedThreadingStruct_wrapper(ThreadingStruct_wrapper):
    def mmbAdd(self):
        if self.mmbID > 0:
            raise MMBError("GappedThreading with id %i already exists in MMB." % self.mmbID)
        call("addGappedThreading", byref(self))

    def mmbUpdate(self):
        call("updateGappedThreading",self.mmbID, byref(self))

    def mmbDelete(self):
        call("deleteGappedThreading", self.mmbID)

    def __str__(self):
        return "gappedThreading %s %s %f" % (self.chainID1, self.chainID2, self.forceConstant)


MMB.addThreading.argtypes = [c_void_p, c_char_p]
MMB.deleteThreading.argtypes = [c_int, c_char_p]
MMB.updateThreading.argtypes = [c_int, c_void_p, c_char_p]

ThreadingStruct_ptr = POINTER(ThreadingStruct_wrapper)
MMB.getThreadingStructs.argtypes = [c_void_p, c_char_p]
MMB.getThreadingStructs.restype = c_int
def getThreadingStructs():
    objs = ThreadingStruct_ptr()
    nbObjs = call('getThreadingStructs', byref(objs))
    return objs[0:nbObjs]

MMB.addGappedThreading.argtypes = [c_void_p, c_char_p]
MMB.deleteGappedThreading.argtypes = [c_int, c_char_p]
MMB.updateGappedThreading.argtypes = [c_int, c_void_p, c_char_p]

GappedThreadingStruct_ptr = POINTER(GappedThreadingStruct_wrapper)
MMB.getGappedThreadingStructs.argtypes = [c_void_p, c_char_p]
MMB.getGappedThreadingStructs.restype = c_int
def getGappedThreadingStructs():
    objs = GappedThreadingStruct_ptr()
    nbObjs = call('getGappedThreadingStructs', byref(objs))
    return objs[0:nbObjs]

###################################################################
class AtomSpring_wrapper(Structure):
    _fields_ = [ 
            ('mmbID', c_int),
            ('atom1', c_char_p),
            ('atom2', c_char_p),
            ('res1', c_int),
            ('res2', c_int),
            ('chain1', c_char_p),
            ('chain2', c_char_p),
            ('toGround', c_bool),
            ('tether', c_bool),
            ('forceConstant', c_double),
            ('deadLength', c_double)
            ]

    def mmbAdd(self):
        if self.mmbID > 0:
            raise MMBError("AtomSpring with id %i already exists in MMB." % self.mmbID)
        call("addAtomSpring", byref(self))

    def mmbUpdate(self):
        call("updateAtomSpring",self.mmbID, byref(self))

    def mmbDelete(self):
        call("deleteAtomSpring", self.mmbID)

    def __str__(self):
        springType = "atomSpring" 
        if self.tether: springType = "atomTether"
        return "%s %s %i %s %s %i %s %f %f" % (springType,
                                               self.chain1, self.res1, self.atom1,
                                               self.chain2, self.res2, self.atom2,
                                               self.deadLength, self.forceConstant)

MMB.addAtomSpring.argtypes = [c_void_p, c_char_p]
MMB.deleteAtomSpring.argtypes = [c_int, c_char_p]
MMB.updateAtomSpring.argtypes = [c_int, c_void_p, c_char_p]

AtomSpring_ptr = POINTER(AtomSpring_wrapper)
MMB.getAtomSprings.argtypes = [c_void_p, c_char_p]
MMB.getAtomSprings.restype = c_int
def getAtomSprings():
    objs = AtomSpring_ptr()
    nbObjs = call('getAtomSprings', byref(objs))
    return objs[0:nbObjs]

###################################################################
class DensityStretch_wrapper(Structure):
    _fields_ = [ 
            ('mmbID', c_int),
            ('chainID', c_char_p),
            ('resStart', c_int),
            ('resEnd', c_int)
            ]

    def mmbAdd(self):
        if self.mmbID > 0:
            raise MMBError("Threading with id %i already exists in MMB." % self.mmbID)
        cmd(str(self))

    def mmbUpdate(self):
        call("updateDensityStretch",self.mmbID, byref(self))

    def mmbDelete(self):
        call("deleteDensityStretch", self.mmbID)

    def __str__(self):
        resStart = str(self.resStart)
        resEnd = str(self.resEnd)
        if self.chainID == "All":
            self.chainID = ""
            resStart = ""
            resEnd = ""
        elif self.resStart == 0:
            resStart = ""
            resEnd = ""
        return "fitToDensity %s %s %s" % (self.chainID, resStart, resEnd)

MMB.addDensityStretch.argtypes = [c_void_p, c_char_p]
MMB.deleteDensityStretch.argtypes = [c_int, c_char_p]
MMB.updateDensityStretch.argtypes = [c_int, c_void_p, c_char_p]

DensityStretch_ptr = POINTER(DensityStretch_wrapper)
MMB.getDensityStretches.argtypes = [c_void_p, c_char_p]
MMB.getDensityStretches.restype = c_int
def getDensityStretches():
    objs = DensityStretch_ptr()
    nbObjs = call('getDensityStretches', byref(objs))
    return objs[0:nbObjs]

###################################################################
MMB.getRootMobilizer.argtypes = [c_char_p, c_char_p, c_char_p]
def getRootMobilizer(chainID):
    rootmobilizer = create_string_buffer("",5)
    call('getRootMobilizer', chainID, rootmobilizer)
    return rootmobilizer.value

MMB.setRootMobilizer.argtypes = [c_char_p, c_char_p, c_char_p]
def setRootMobilizer(chainID, rootmobilizer):
    call('setRootMobilizer', chainID, rootmobilizer)



