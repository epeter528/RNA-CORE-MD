## @package docstring
# MMB UI Chimera plugin.
# 2013
# Author: Alex Tek - ICM, Uppsala University

import chimera
import Midas
from chimera import misc, OpenModels, Bond, MolResId
from chimera import colorTable

from Movie.gui import MovieDialog

import threading
import time
import sys
import os
import platform
from collections import OrderedDict

##################################################################################################
# pyMMB initialization
mmbUIPath = os.path.dirname(__file__)
workDir = mmbUIPath
print "MMB INIT ****************** "
print mmbUIPath

if platform.system() == "Darwin":
    os.environ['DYLD_FALLBACK_LIBRARY_PATH'] += ":"+mmbUIPath
elif platform.system() == "Windows":
    os.environ['PATH'] += ";"+mmbUIPath

import pyMMB
leontisWesthofFileName = os.path.join(mmbUIPath,"parameters.csv")
MMBparameters = pyMMB.MMBparameters

## MMB sessions number
MMB_SESSION_NUMBER = 0

## Flag to know if the polymers are initialized or not
polymersInitialized = False

## Chimera model containing the current biopolymers 
currentModel = None

## Chimera's Pseudobond group containing base pairs bonds to display 
basePairBonds = None

## Chimera's Pseudobond group containing constraints bonds to display
constraintsBonds = None

## Chimera's Pseudobond group containing atomSprings bonds to display
atomSpringsBonds = None

## Chimera MDMovie dialog used to navigate among simulation frames
movieDialog = None

## Ensemble object containing frames information
ensemble = None

## Current MMB simulation
currentSimulation = None

## Simulation flags
initFlag = False
runFlag = False

## Current BioPolymers in MMB. Key is the chainID
polymers = {}

def setWorkDir( newWD ):
    global workDir
    workDir = newWD
    os.chdir(workDir)

## Synchronise polymers sequences with the ones in MMB
def refreshSequences():
    polymers.clear()
    seqs = pyMMB.getSequences()
    for s in seqs:
        if s.chainID in polymers:
            raise pyMMB.MMBError("Invalid sequence: " + s.chainID + " already exists.")
        BioPolymer.mutateBioPolymerSeq_wrapper(s, currentModel)
        polymers[s.chainID] = s

## Use pyMMB to extract sequences from a pdb files
#  Populate the polymers global dict
#  @param the path of the pdb file
def loadSequencesFromPdb(pdbFileName):
    pyMMB.loadSequencesFromPdb(pdbFileName)
    refreshSequences();

def importSequencesFromChimera(models=[],selectedAtoms=[]):
    ficName = '/tmp/chimeraModels.pdb'
    if platform.system() == "Windows":
        ficName = os.path.join(os.environ["TEMP"],"chimeraModels.pdb")
    selOnly = bool(selectedAtoms)
    molecules = [m for m in chimera.openModels.list() if isinstance(m,chimera.Molecule)]
    models = [m for m in models if m]
    print models
    if models:
        molecules = [m for m in molecules if m in models]
    print molecules
    # for m in molecules:
        # print m.coordSets
        # print m.activeCoordSet
    if molecules:
        chimera.pdbWrite(molecules, 
                     chimera.Xform.identity(),
                     # molecules[0].openState.xform, 
                     ficName,
                     selectedOnly=selOnly,
                     selectionSet=set(selectedAtoms))
        pyMMB.loadSequencesFromPdb(ficName)
        refreshSequences();

def initMMB():
    global MMBparameters
    global polymersInitialized, currentModel, MMB_SESSION_NUMBER 
    global basePairBonds, constraintsBonds
    global initFlag, runFlag, currentSimulation
    global ensemble, movieDialog

    pyMMB.initMMB(leontisWesthofFileName)
    MMBparameters = pyMMB.MMBparameters
    pyMMB.cmd("lastStage 1")
    polymersInitialized = False

    if currentSimulation:
        currentSimulation.stop()
        currentSimulation = None
        # print [currentModel]
        # importSequencesFromChimera(models=[currentModel])

    if ensemble:
        ensemble = None
    if movieDialog:
        # movieDialog.destroy()
        movieDialog = None
    initFlag = False
    runFlag = False

    if currentModel:
        currentModel.name = "mmb_"+str(MMB_SESSION_NUMBER)
        currentModel = None
    if basePairBonds:
        basePairBonds.destroy()
        # chimera.openModels.close(basePairBonds)
        basePairBonds = None
    if constraintsBonds:
        constraintsBonds.destroy()
        # chimera.openModels.close(constraintsBonds)
        constraintsBonds = None
    MMB_SESSION_NUMBER += 1


## Change one MMB parameter from its name
def setMMBParameter(name, value):
    # print name, value
    # paramType = eval("type(MMBparameters.%s)" % name)
    MMBparameters.__setattr__(name, value)

## Get MMB parameter from its name
def getMMBParameter(name):
    return MMBparameters.__getattribute__(name)

initMMB()

##################################################################################################
# Trajectory

## Defines a Molecular model containing several conformations
# 
#  It is used by the MDMovie Chimera plugin 
class MolEnsemble:
    ## Constructor
    #  @param mol: a Chimera molecule
    def __init__(self, mol):
        self.molecule = mol
        self.name = "mmb_"+str(MMB_SESSION_NUMBER)
        self.startFrame = 1
        self.endFrame = len(mol.coordSets)

    ## The length of the ensemble is the number of conformations
    #  @return the number of conformations in this ensemble
    def __len__(self):
        return len(self.molecule.coordSets)

## Subclass of Chimera MovieDialog implementing frame management methods
# 
class MyMovieDialog(MovieDialog):
    ## Constructor
    #  @param ensemble a MolEnsemble object
    #  @kw a list of optional named parameters for MovieDialog instanciation
    def __init__(self, ensemble, **kw):
        MovieDialog.__init__(self, ensemble, **kw)

    ## Synchronize the dialog display with the ensemble content.
    def updateFrames(self):
        ensembleSize = len(self.ensemble)
        self.endFrame = ensembleSize
        self.triggers.activateTrigger(self.MORE_FRAMES, ensembleSize)

    ## Add a set of coordinates to the ensemble object and update the dialog display
    #  @param coordArray numpy array containing coordinates to add to the model
    def addFrame(self, coordArray):
        cs = self.ensemble.molecule.newCoordSet(len(self.ensemble)+1)
        chimera.fillCoordSet(cs, self.ensemble.molecule.atoms, coordArray)
        self.ensemble.molecule.activeCoordSet = cs
        self.updateFrames()
        # self.LoadFrame(len(self.ensemble))
        # print "wait"
        # Midas.wait(1)


##################################################################################################
# TODO: Move the headers to the respective locations in gui.py
## Headers names for mobilizers
mobilizersHeaders=("Id","Mobilizer","Mobility","Chain","Start Res","End Res","Valid","Select")

## Headers names for mobilizers within
mobilizersWithinHeaders=("Id","Mobilizer","Mobility","Chain","Res","Radius","Valid","Select")

## Headers names for constraints
constraintsHeaders=("Id","Chain1","Res1","Atom1"," ", "Chain2", "Res2", "Atom2", "Valid", "Select")


## Colors 
# validationColors = {True:"green", False:"red"}
# validationStates = {True:"disable", False:"normal"}

## Types of polymers. Indices match MMB types enum
PolyTypes = ['RNA','protein','DNA','Unassigned']

## Indicate how to initialize coordinates. 
# New: Use MMB coordinates generator.
# PDB: Extract coordinates from a PDB file.
Topologies = ['New','PDB']

## Types of pair interaction according to Leontis & Westhof.
PairTypes = ["WatsonCrick", "Hoogsteen", "SugarEdge", "HelicalStackingA3", "HelicalStackingA5"]
## Possible bond orientations according to Leontis & Westhof.
BondOrient = ["Cis", "Trans"]

## Default parameters for a Sequence.
Defaults = {
    "chainID": 'X',
    "sequence": '',
    "firstResNum": 1,
    "polyType":0,
    "topo":'New'
}

## Unicode representation of Leontis & Westhof interactions
LeontisWesthofSymbols = {
    ("WatsonCrick", "Cis"):  u"\u25CF",
    ("WatsonCrick", "Trans"):u"\u25CB",
    ("Hoogsteen",   "Cis"):  u"\u25A0",
    ("Hoogsteen",   "Trans"):u"\u25A1",
    ("SugarEdge",   "Cis"):  u"\u25B6",
    ("SugarEdge",   "Trans"):u"\u25B7"
}

allConstraintsRestraints = OrderedDict()

## Mobility Types
MobilityTypes = ["Default", "Rigid", "Free", "Torsion"]

## List of current base pairs, synchronized with MMB
basePairs = []
allConstraintsRestraints["basePairs"] = basePairs

## List of current mobilizers, synchronized with MMB
mobilizers = []
allConstraintsRestraints["mobilizers"] = mobilizers

## List of current mobilizersWithin, synchronized with MMB
mobilizersWithin = []
allConstraintsRestraints["mobilizersWithin"] = mobilizersWithin

## List of current rootMobilizers, synchronized with MMB
rootMobilizers = []
allConstraintsRestraints["rootMobilizers"] = rootMobilizers

## List of current contacts, synchronized with MMB
contacts = []
allConstraintsRestraints["contacts"] = contacts

## List of current contactsWithin, synchronized with MMB
contactsWithin = []
allConstraintsRestraints["contactsWithin"] = contactsWithin

## List of constraints, synchronized with MMB
constraints = []
allConstraintsRestraints["constraints"] = constraints

## List of atomSprings, synchronized with MMB
atomSprings = []
allConstraintsRestraints["atomSprings"] = atomSprings

## List of current AllResiduesWithin, synchronized with MMB
allResiduesWithins = []
allConstraintsRestraints["allResiduesWithins"] = allResiduesWithins

## List of current includeAllNonBondAtomsInResidues, synchronized with MMB
includeAllNonBondAtomsInResidues = []
allConstraintsRestraints["includeAllNonBondAtomsInResidues"] = includeAllNonBondAtomsInResidues

## List of current Threadings, synchronized with MMB
threadings = []
allConstraintsRestraints["threadings"] = threadings

## List of current GappedThreadings, synchronized with MMB
gappedThreadings = []
allConstraintsRestraints["gappedThreadings"] = gappedThreadings

## List of current density stretches, synchronized with MMB
densities = []
allConstraintsRestraints["densities"] = densities

## Current density map model
densityMapModel = None

## Name of the pdb file chosen by the user.
currentPdb = ""


##################################################################################################
## Defines a generic MMB object (mother class for interactions, constraints etc.)
#
class MMB_UI_Object:
    def __init__(self, model=currentModel, mmbID=-1, new=False):
        self.model = model
        self.mmbID = mmbID
        self.new   = new
        self.show  = True

    def setAttribute(self, attrName, value):
        if getattr(self, attrName) == value:
            return
        if attrName.startswith("res"):
            value = int(value.split()[1])
        setattr(self, attrName, value)

    def mmbAdd(self):
        if self.new:
            pyMMB.cmd(str(self))
        else:
            pass

    def mmbUpdate(self):
        pass

    def mmbDelete(self):
        pass

    def updateRepresentation(self):
        pass

    def updateRepresentationColor(self):
        pass

    def getChimeraSelection(self):
        return [None]

    def removeRepresentation(self):
        res = self.getChimeraSelection()
        if res:
            for r in res:
                if r:
                    r.label = ""
                    r.labelColor = None
                    r.ribbonColor = None
                    for a in r.atoms:
                        a.color = None

    ## Add a constraint to current chimera selection
    #  @param constraint to select
    def addToSelection(self):
        sel = self.getChimeraSelection()
        if sel:
            chimera.selection.addCurrent([s for s in sel if s])

    ## Remove a constraint from current chimera selection
    #  @param constraint to deselect
    def removeFromSelection(self):
        sel = self.getChimeraSelection()
        if sel:
            chimera.selection.removeCurrent([s for s in sel if s])


##################################################################################################
## Represents a biopolymer as in MMB
# 
class BioPolymer(pyMMB.BioPolymerSeq, MMB_UI_Object):
    representativeAtom = {"RNA": "C4'","protein": "CA","DNA": "C4'"}    
    
    ## Constructor
    def __init__(self,chainID='',sequence='', firstResNum=1, polyType=3,pdbFileName='',pdbChainID='', loadFromPdb=False, activePhysics=True, new=False, valid=True):
        self.chainID     = chainID
        self.sequence    = sequence
        self.firstResNum = firstResNum
        self.polyType    = polyType
        self.pdbFileName = pdbFileName
        self.pdbChainID  = pdbChainID
        self.loadFromPdb = loadFromPdb
        self.activePhysics = activePhysics
        self.initialization(currentModel, new, valid)

    def initialization(self, model, new, valid):
        self.topo        = "New"
        self.show  = True

        if self.loadFromPdb:
            self.topo = "PDB"

        self.new = new
        if self.new:
            self.topo = "New"
            self.loadFromPdb = False
        self.valid = valid
        self.select = False

        self.idsList = None

        # self.lastResNum = self.firstResNum + len(self.sequence) - 1

    @classmethod
    ## Mutate a MobilizerWithin_wrapper object's type to MobilizerWithin
    def mutateBioPolymerSeq_wrapper(cls, BioPolymerSeq_wrapper_object, model=currentModel, new=False, valid=True):
        BioPolymerSeq_wrapper_object.__class__ = cls 
        BioPolymerSeq_wrapper_object.initialization(model, new, valid)

    ## Length of the polymer == length of the sequence 
    #  @return the sequence's length
    def __len__(self):
        return len(sequence)

    def setAttribute(self, attrName, value):
        if getattr(self, attrName) == value:
            return
        if attrName == "polyType":
            setattr(self, attrName, PolyTypes.index(value))
            return
        if attrName in ["topo","sequence"]:
            setattr(self, attrName, value)
            if self.topo == "New":
                self.loadFromPdb = False
            else:
                self.loadFromPdb = True
            return
        setattr(self, attrName, value)

    # def mmbUpdate(self):
    #     pyMMB.updatePolymer(self.chainID, self.sequence, self.pdbFileName, self.loadFromPdb)

    def mmbDelete(self):
        pyMMB.deletePolymer(self.chainID)

    ## Add the polymer to MMB
    def mmbAdd(self):
        pyMMB.cmd(PolyTypes[self.polyType]+" "+self.chainID+" "+str(self.firstResNum)+" "+self.sequence)

    ## Returns the representative atom name
    def getRepresentativeAtom(self):
        return BioPolymer.representativeAtom[PolyTypes[self.polyType]]

    ## Builds a list of ids
    #  @return a list of ids
    def getIdsList(self, start=-1, end=sys.maxint):
        # if not self.idsList:
        self.idsList = []
        for res in currentModel.residues:
            if res.id.chainId == self.chainID and res.id.position >= start and res.id.position < end:
                self.idsList.append(res.id.position)
        self.idsList = sorted(self.idsList)
        return self.idsList

    ## Builds a string composed of the residue's 1-letter code and id 
    #  @param id id of the wanted residue
    #  @return a string like this: A 23
    def getResNameId(self, id):
        return self.getResName(id) + " " + str(id)

    ## Returns the residue id 1-letter name or " " if the id is outside the sequence's ids range
    #  @param id id of the wanted residue
    #  @return  1-letter name or " "
    def getResName(self, id):
        idsList = self.getIdsList()
        if id in idsList:
            return self.sequence[idsList.index(id)]
        else:
            return ""

    ## Build a list of resname + id strings
    #  @return a list containing "resname id" strings
    def getResNamesIDsList(self, start=-1):
        z = zip(self.sequence, self.getIdsList(start))
        return [rName+" "+str(rId) for rName, rId in z] 

    ## Return a list of Chimera Residue objects from start to end (included)
    #  @param start id of first residue
    #  @param end id of last residue (included)
    #  @return a list of Chimera Residue objects
    def getChimeraResidues(self, start=None, end=None):
        if not currentModel:
            return []

        if not end: end = sys.maxint

        selAdd = []    
        for res in currentModel.residues:
            if res.id.chainId == self.chainID:
                if res.id.position >= start and res.id.position <= end:
                    selAdd.append(res)
        return selAdd

    def updateRepresentation(self):
        pass

    def getChimeraSelection(self):
        return self.getChimeraResidues()

    def getRootMobilizer(self):
        return pyMMB.getRootMobilizer(self.chainID)

    def matchCoordinatesFromChimera(self):
        import cStringIO
        out = cStringIO.StringIO()
        chimera.pdbWrite([currentModel], currentModel.openState.xform, out)
        pyMMB.matchCoordinatesFromContent(self.chainID, out) 

    def show_hide(self):
        modelID = currentModel.id
        # print ":.%s#%i" % (self.chainID, modelID)
        if self.show:
            self.show = False
            # chimera.runCommand("~display #%i:.%s" % (modelID, self.chainID))
            chimera.runCommand("~ribbon #%i:.%s" % (modelID, self.chainID))
        else:
            self.show = True
            # chimera.runCommand("display #%i:.%s" % (modelID, self.chainID))
            chimera.runCommand("ribbon #%i:.%s" % (modelID, self.chainID))

## Check if all polymers have been validated in MMB
def checkSequences():
    for p in polymers.values():
        if not p.valid:
            return False
    return True

## Validate all polymers that needs it
def validateSequences():
    for p in polymers.values():
        if not p.valid:
            p.mmbUpdate()
            p.valid = True

## Initialize the polymers with the coordinates contained in the pdb file.
#  Open a new model in Chimera containing all the biopolymers.
#  @param polymers dictionary of Biopolymer objects with chainIDs as keys
#  @param the path of the pdb file 
def initPolymers(polymers, pdbFileName):
    global polymersInitialized
    #pyMMB.clearPolymers()
    for p in polymers.values():
        pdbfile = p.pdbFileName
        if not p.loadFromPdb:
            pdbfile = ""
        pyMMB.initializePolymer(p.chainID, pdbfile)
    polymersInitialized = True

def loadPolymers():
    global currentModel
    if os.getcwd() == "/"    :
        os.chdir (os.environ['HOME']);
    ficName = str(MMBparameters.workingDirectory) + '/mmb.pdb'
    #ficName = 'mmb.pdb'
    pyMMB.writeDefaultPdb(ficName)
    currentModel = chimera.openModels.open(ficName)[0]
    # print "CurrentModel coordSets: "
    # print currentModel.coordSets
    for m in chimera.openModels.list():
        if m.id != currentModel.id:
            m.display = False

## Module method returning a list of the polymers sorted by chain Id
#  @return a list of the biopolymers sorted by chain ID
def sortedPolymers():
    return [polymers[k] for k in sorted(polymers.keys())]

## Module method returning a sorted list of the biopolymer's chain ids.
def sortedChainsIDs():
    return sorted(polymers.keys())

## Module method returning a sorted list of the nucleic acids's chain ids.
def sortedNucleicAcidsIDs():
    return sorted([k for k in polymers.keys() if PolyTypes[polymers[k].polyType] in ("RNA","DNA") ])

## Return a list of atoms ids and names for a residue
def getAtomsList(chainID, resID, model=currentModel):
    res = model.findResidue(MolResId(chainID, resID))
    if not res:
        return []
    return [a.name for a in res.atoms]

## Return Chimera residues within a radius around the CA or C4 of :resID.chaindID
def getResiduesWithin(chainID, resID, radius):
    repAtom = polymers[chainID].getRepresentativeAtom()
    chimera.runCommand("sel :%i.%s@%s z<%.1f & @%s" % (resID, chainID, repAtom, radius, repAtom))
    return chimera.selection.currentResidues()

##################################################################################################
## Defines a base pair interaction as in MMB
# 
#  Some methods are used as callbacks for Tk widgets.
#  Others are used to get information from Chimera
class BasePairInteraction(pyMMB.BaseInteraction_wrapper, MMB_UI_Object):
    ## Constructor
    def __init__(self, model, mmbID, poly1, res1, edge1, poly2, res2, edge2, bondOrient, 
                 select=False, valid=True, newPair=False, chimBond=None):
        ## the Chimera model containing the residues of the interaction.
        MMB_UI_Object.__init__(self, model, mmbID)
        self.poly1      = poly1
        self.res1       = res1
        self.edge1      = edge1
        self.poly2      = poly2
        self.res2       = res2
        self.edge2      = edge2
        self.bondOrient = bondOrient
        self.initialization(model, newPair, valid, select, chimBond)
        
    def initialization(self, model=currentModel, new=False, valid=True, select=False, chimBond=None):
        self.model      = model
        self.select     = select
        self.valid      = valid
        self.new        = new
        ## The Chimera pseudobond used to display the interaction
        self.chimeraBond= chimBond


    @classmethod
    ## Mutate a MobilizerWithin_wrapper object's type to MobilizerWithin
    def mutateBaseInteraction_wrapper(cls, BaseInteraction_wrapper_object, model=currentModel, new=False, valid=True):
        BaseInteraction_wrapper_object.__class__ = cls 
        BaseInteraction_wrapper_object.initialization(model, new, valid)

    ## Build an interaction with default values
    #  @return a BasePairInteraction object with default values
    @classmethod
    def emptyPair(cls):
        ch1 = ""
        ch2 = ""
        try:
            ch1 = sortedNucleicAcidsIDs()[0]
            ch2 = sortedNucleicAcidsIDs()[0]
        except IndexError:
            pass

        return cls(currentModel,
                   -1, 
                   ch1, 0, PairTypes[0], 
                   ch2, 0, PairTypes[0], 
                   BondOrient[0], newPair=True)

    ## Return Chimera residue objects of the two residues
    #  @return a tuple containing the two residues. An element can be a NoneType.
    def getChimeraResidues(self):
        res1 = self.model.findResidue(MolResId(self.poly1, self.res1))
        res2 = self.model.findResidue(MolResId(self.poly2, self.res2))
        return res1, res2

    ## Return Chimera atom objects for the C5' atoms of the residues
    #  @return a tuple containing the two atomss. An element can be a NoneType.
    def getChimeraAtoms(self):
        res1,res2 = self.getChimeraResidues()
        atom1 = None
        atom2 = None
        if res1:
            atom1 = res1.findAtom(polymers[self.poly1].getRepresentativeAtom())
        if res2:
            atom2 = res2.findAtom(polymers[self.poly2].getRepresentativeAtom())
        return atom1, atom2

    ## Return a list of chimera selectables for the residues and the bond
    #  @return a list of chimera selectables. Can contain None values
    def getChimeraSelection(self):
        return list(self.getChimeraResidues()) + [self.chimeraBond]

    ## Build a string composed of the residues 1-letter name and the Leontis&Westhof symbol.
    #  @return a string like: A O U
    def getLabel(self):
        res1 = polymers[self.poly1].getResName(self.res1)
        res2 = polymers[self.poly2].getResName(self.res2)

        if self.edge1 == self.edge2:
            symbol = ""
            if (self.edge1, self.bondOrient) in LeontisWesthofSymbols:
                symbol = LeontisWesthofSymbols[(self.edge1, self.bondOrient)]
            return res1 + " - " + symbol + " - " + res2

        symbol1 = ""
        symbol2 = ""
        if (self.edge1, self.bondOrient) in LeontisWesthofSymbols:
            symbol1 = LeontisWesthofSymbols[(self.edge1, self.bondOrient)]
        if (self.edge2, self.bondOrient) in LeontisWesthofSymbols:
            symbol2 = LeontisWesthofSymbols[(self.edge2, self.bondOrient)]
        return res1 + " - " + symbol1 + " " + symbol2 + " - " + res2

    ## Update the bond display in Chimera. Create a new bond if necessary.
    def updateRepresentation(self):
        global basePairBonds

        b = self.chimeraBond
        if b:
            b.reuse(*(self.getChimeraAtoms()))
            b.label = self.getLabel()
            self.updateRepresentationColor()
            return
        # if new bond, check atoms
        a1, a2 = self.getChimeraAtoms()
        if a1 and a2:
            b = basePairBonds.newPseudoBond(a1,a2)
            b.display = 1
            b.drawMode = Bond.Wire
            b.label = self.getLabel()
            self.chimeraBond = b
            self.updateRepresentationColor()

    ## Update the bond's color according to state of the interaction in the GUI.
    #  New interaction: blue
    #  Non validated interaction: orange
    #  Normal state: black
    def updateRepresentationColor(self):
        if self.new:
            self.chimeraBond.color = colorTable.getColorByName("blue")
            return
        if not self.valid:
            self.chimeraBond.color = colorTable.getColorByName("orange")
            return
        
        self.chimeraBond.color = colorTable.getColorByName("black")

## Synchronize the base pair interactions with MMB
#  Update bond display as well
def refreshBasePairInteractions():
    """
    String format:
    MMBid ch1 resN1 edge1 ch2 resN2 edge2 orientation 
    """
    global basePairBonds

    interactions = pyMMB.getBaseInteractions()
    basePairs[:] = []
    for i in interactions:
        BasePairInteraction.mutateBaseInteraction_wrapper(i, currentModel)
        basePairs.append(i)

    refreshBasePairBonds()

## Update base pair bonds display
def refreshBasePairBonds():
    global basePairBonds

    if basePairBonds:
        basePairBonds.deleteAll()
    else:
        basePairBonds = misc.getPseudoBondGroup('BasePairs', modelID=OpenModels.Default, hidden=False)
        basePairBonds.name = "MMB Base Pairs"
        basePairBonds.lineType = chimera.Dash
        basePairBonds.lineWidth = 3
        basePairBonds.color = colorTable.getColorByName("black")
    
    [p.updateRepresentation() for p in basePairs]

##################################################################################################
## Defines a Nucleic Acid Duplex. Only for the gui as no such data structure exists in MMB.
#  Some methods are used as callbacks for Tk widgets.
#  Others are used to get information from Chimera
class NucleicAcidDuplex(MMB_UI_Object):
    ## Constructor
    def __init__(self, model=currentModel, mmbID=-1, chain1="", resStart1=0, resEnd1=0, 
                                                     chain2="", resStart2=0, resEnd2=0, new=True):
        MMB_UI_Object.__init__(self, model, mmbID)
        self.chainID1      = chain1
        self.residueStart1 = resStart1
        self.residueEnd1   = resEnd1
        self.chainID2      = chain2
        self.residueStart2 = resStart2
        self.residueEnd2   = resEnd2
        self.initialization(model, new)

    def initialization(self, model=currentModel, new=True, selection=None, valid=True):
        self.model      = model
        self.new        = new
        self.selection  = selection
        self.valid      = valid

    def mmbAdd(self):
        cmdTxt = "nucleicAcidDuplex %s %i %i %s %i %i" % ( self.chainID1,
                                                        self.residueStart1,
                                                        self.residueEnd1,
                                                        self.chainID2,
                                                        self.residueStart2,
                                                        self.residueEnd2
                                                       )
        pyMMB.cmd(cmdTxt)

    ## Initialize the residue when changing chain
    def setAttribute(self, attrName, value):
        if attrName == "chain1":
            setattr(self, "residueStart1", 0)
            setattr(self, "residueStart2", 0)
        elif attrName == "chain2":
            setattr(self, "residueStart1", 0)
            setattr(self, "residueStart2", 0)
        MMB_UI_Object.setAttribute(self, attrName, value)

    def getResStart1NameId(self):
        if self.residueStart1 == 0:
            return ""
        poly = polymers[self.chainID1]
        return poly.getResNameId(self.residueStart1)
    def getResStart2NameId(self):
        if self.residueStart2 == 0:
            return ""
        poly = polymers[self.chainID2]
        return poly.getResNameId(self.residueStart2)

    def getResEnd1NameId(self):
        if self.residueEnd1 == 0:
            return ""
        poly = polymers[self.chainID1]
        return poly.getResNameId(self.residueEnd1)
    def getResEnd2NameId(self):
        if self.residueEnd2 == 0:
            return ""
        poly = polymers[self.chainID2]
        return poly.getResNameId(self.residueEnd2)

    def getChimeraSelection(self):
        if self.chainID1 == ""  and self.chainID2 == "":
            return []
        sel = []
        poly = polymers[self.chainID1]
        sel += poly.getChimeraResidues(self.residueStart1, self.residueEnd1)

        if self.chainID2:
            poly = polymers[self.chainID2]
            sel += poly.getChimeraResidues(self.residueStart2, self.residueEnd2)
        return sel


##################################################################################################
mobiColors = {"Default":"orange", "Rigid":"dim gray", "Free":"yellow", "Torsion":"green"}
## Defines a MMB's mobilizer
#  Some methods are used as callbacks for Tk widgets.
#  Others are used to get information from Chimera
class Mobilizer(pyMMB.MobilizerStretch_wrapper, MMB_UI_Object):
    ## Constructor
    def __init__(self, model=currentModel, mmbID=-1, mobility="Default", chain="All", resStart=0, resEnd=-1, new=True):
        MMB_UI_Object.__init__(self, model, mmbID)
        self.mobility   = mobility
        self.chainID    = chain
        self.resStart   = resStart
        self.resEnd     = resEnd
        self.initialization(model, new)

    def initialization(self, model=currentModel, new=False, valid=True):
        self.model = model
        self.valid = valid
        self.new   = new

    @classmethod
    ## Mutate a MobilizerWithin_wrapper object's type to MobilizerWithin
    def mutateMobilizerStretch_wrapper(cls, MobilizerStretch_wrapper_object, model=currentModel, new=False, valid=True):
        MobilizerStretch_wrapper_object.__class__ = cls 
        MobilizerStretch_wrapper_object.initialization(model, new, valid)

    def setAttribute(self, attrName, value):
        if getattr(self, attrName) == value:
            return

        if attrName == "resStart" and value == "All":
            self.resStart = 0
            self.resEnd = -1
            return

        if attrName == "resStart":
            value = int(value.split()[1])
            self.resStart = value
            if value > self.resEnd:
                self.resEnd = value
            return

        if attrName.startswith("res"):
            value = int(value.split()[1])
        setattr(self, attrName, value)

    def getResStartNameId(self):
        if self.resStart == 0:
            return "All"
        poly = polymers[self.chainID]
        return poly.getResNameId(self.resStart)

    def getResEndNameId(self):
        if self.resEnd == -1:
            return ""
        poly = polymers[self.chainID]
        return poly.getResNameId(self.resEnd)

    ## Return a list of Chimera residue objects
    #  @return a list containing the residues.
    def getChimeraSelection(self):
        if self.chainID == "All":
            return currentModel.residues;
        poly = polymers[self.chainID]
        if self.resStart == "All":
            return poly.getChimeraResidues(start=poly.firstResNum)
        return poly.getChimeraResidues(self.resStart, self.resEnd)

    def updateRepresentation(self):
        resSel = self.getChimeraSelection()
        # for r in resSel:
        #     r.ribbonDisplay = True
        self.updateRepresentationColor(resSel)

    ## Update the color according to the mobility
    def updateRepresentationColor(self, resSel):
        color = colorTable.getColorByName(mobiColors[self.mobility])
        for r in resSel:
            r.ribbonColor = color
            for a in r.atoms:
                a.color = color

## Synchronize the mobilizers list with MMB
def refreshMobilizerStretches():
    global mobilizers

    mmbMobilizers = pyMMB.getMobilizerStretches()
    mobilizers[:] = []
    for m in mmbMobilizers:
        Mobilizer.mutateMobilizerStretch_wrapper(m, currentModel, False)
        mobilizers.append(m)

##################################################################################################
## Defines a MMB's mobilizerWithin
#  Some methods are used as callbacks for Tk widgets.
#  Others are used to get information from Chimera
class MobilizerWithin(pyMMB.MobilizerWithin_wrapper, MMB_UI_Object):
    ## Constructor
    def __init__(self, model=currentModel, mmbID=-1, mobility="Default", chain="", res=0, radius=0.4, new=True):
        MMB_UI_Object.__init__(self, model, mmbID)
        self.mobility   = mobility
        self.radius     = radius
        self.chainID    = chain
        self.resID      = res
        self.initialization(model, new)

    def initialization(self, model=currentModel, new=True, selection=None, valid=True):
        self.model      = model
        self.new        = new
        self.selection  = selection
        self.valid      = valid

    @classmethod
    ## Mutate a MobilizerWithin_wrapper object's type to MobilizerWithin
    def mutateMobilizerWithin_wrapper(cls, MobilizerWithin_wrapper_object, model=currentModel, new=False, selection=None, valid=True):
        MobilizerWithin_wrapper_object.__class__ = cls 
        MobilizerWithin_wrapper_object.initialization(model, new, selection, valid)

    ## Initialize the residue when changing chain
    def setAttribute(self, attrName, value):
        if attrName.startswith("chain"):
            setattr(self, "resID", 0)
        MMB_UI_Object.setAttribute(self, attrName, value)

    def getResNameId(self):
        if self.resID == 0:
            return ""
        poly = polymers[self.chainID]
        return poly.getResNameId(self.resID)

    def getChimeraSelection(self):
        if self.resID != 0 and self.radius > 0:
            # resWithin = pyMMB.getResiduesWithin(self.chainID, self.resID, self.radius)
            # return [self.model.findResidue( MolResId(x[0], x[1])) for x in resWithin]
            # chimera.selection.addCurrent(self.selection)
            return getResiduesWithin(self.chainID, self.resID, self.radius*10.0)

    def updateRepresentation(self):
        resSel = self.getChimeraSelection()
        self.updateRepresentationColor(resSel)

    ## Update the color according to the mobility
    def updateRepresentationColor(self, resSel):
        if resSel == None:
            return
        color = colorTable.getColorByName(mobiColors[self.mobility])
        for r in resSel:
            r.ribbonColor = color
            for a in r.atoms:
                a.color = color

## Synchronize the mobilizersWithin list with MMB
def refreshMobilizerWithin():
    global mobilizersWithin

    mmbMobilizers = pyMMB.getMobilizersWithin()
    mobilizersWithin[:] = []
    for m in mmbMobilizers:
        MobilizerWithin.mutateMobilizerWithin_wrapper(m, currentModel)
        mobilizersWithin.append(m)


##################################################################################################
## Defines a constraint as in MMB
#
class Constraint(pyMMB.Constraint_wrapper, MMB_UI_Object):
    ##Constructor
    def __init__(self, model=currentModel, mmbID=-1,
                 chain1="", res1=0, atom1="", 
                 chain2="Ground", res2=0, atom2="", 
                 new=True, valid=True):
        MMB_UI_Object.__init__(self, model, mmbID)
        self.chain1     = chain1
        self.res1       = res1
        self.atom1      = atom1.replace("*","'")
        self.chain2     = chain2
        self.res2       = res2
        self.atom2      = atom2.replace("*","'")
        self.initialization(model, new, valid)

    def initialization(self, model=currentModel, new=False, valid=True):
        self.model       = model
        self.atom1       = self.atom1.replace("*","'")
        self.atom2       = self.atom2.replace("*","'")
        self.valid       = valid
        self.new         = new
        self.chimeraBond = None

        if not self.chain2.isalpha():
            self.chain2 = "Ground"

        # if self.res1 == 0:
        #     self.res1 = ""
        # if self.res2 == 0:
        #     self.res2 = ""

    @classmethod
    ## Mutate a MobilizerWithin_wrapper object's type to MobilizerWithin
    def mutateConstraint_wrapper(cls, Constraint_wrapper_object, model=currentModel, new=False, valid=True):
        Constraint_wrapper_object.__class__ = cls 
        Constraint_wrapper_object.initialization(model, new, valid)

    ## Return an MMB command for this constraint
    def __str__(self):
        # if self.chain1=="":
        #     return "rootMobilizer Weld"
        # if self.res1=="":
        #     return "rootMobilizer " + self.chain1 + " Weld"
        com = "constraint " + self.chain1 + " " + str(self.res1) + " " + self.atom1
        if self.chain2 == "Ground":
            return "constrainToGround " + self.chain1 + " " + str(self.res1)
        return com + " Weld " + self.chain2 + " " + str(self.res2) + " " + self.atom2

    def setAttribute(self, attrName, value):
        if getattr(self, attrName) == value:
            return
        # if the attribute is res1 or res2 we set the corresponding atom to ""
        if attrName.startswith("res"):
            nAttr = attrName.split("s")[1]
            atomAttr = "atom"+nAttr
            setattr(self, atomAttr, "")
            if value:
                value = int(value.split()[1])
            else:
                value = 0
        # if the attribute is chain1 or chain2 we set the corresponding res,atom to 0,""
        elif attrName.startswith("chain"):
            nAttr = attrName.split("n")[1]
            setattr(self, "res"+nAttr, 0)
            setattr(self, "atom"+nAttr, "")
        setattr(self, attrName, value)

    def getRes1AtomsList(self):
        return getAtomsList(self.chain1, self.res1, self.model)

    def getRes2AtomsList(self):
        return getAtomsList(self.chain2, self.res2, self.model)

    def getRes1NameId(self):
        if self.res1 == 0:
            return ""
        return polymers[self.chain1].getResNameId(self.res1)

    def getRes2NameId(self):
        if self.chain2 == "Ground" or self.res2 == 0:
            return ""
        return polymers[self.chain2].getResNameId(self.res2)
        
    def getChimeraSelection(self):
        if self.chain1 == 0 or self.res1 == 0:
            return

        res1 = self.model.findResidue(MolResId(self.chain1, self.res1))
        selection = [self.chimeraBond]
        if self.atom1 != "":
            selection.append(res1.findAtom(self.atom1))
        else:
            selection.append(res1)
        if self.res2 == 0:
            return selection
        
        res2 = self.model.findResidue(MolResId(self.chain2, self.res2))
        if self.atom2 != "":
            selection.append(res2.findAtom(self.atom2))
        else:
            selection.append(res2)
        return selection

    ## Update the bond display in Chimera. Create a new bond if necessary.
    def updateRepresentation(self):
        global constraintsBonds

        # Contraint to ground
        if self.chain2 == "Ground" and self.res1:
            poly = polymers[self.chain1]
            resSel = poly.getChimeraResidues(self.res1, self.res1)
            for r in resSel:
                r.label = u"\u22A5"
            self.updateRepresentationColor()
            return

        # No bond
        if self.res2 == 0:
            if self.chimeraBond: 
                constraintsBonds.deletePseudoBond(self.chimeraBond)
            self.chimeraBond = None
            return

        res1 = self.model.findResidue(MolResId(self.chain1, self.res1))
        res2 = self.model.findResidue(MolResId(self.chain2, self.res2))
        a1 = self.atom1
        a2 = self.atom2
        if self.atom1 == "":
            a1 = polymers[self.chain1].getRepresentativeAtom()
        if self.atom2 == "":
            a2 = polymers[self.chain2].getRepresentativeAtom()

        b = self.chimeraBond
        if b:
            b.reuse(res1.findAtom(a1), res2.findAtom(a2))
            b.label = "||"
            self.updateRepresentationColor()
            return

        b = constraintsBonds.newPseudoBond(res1.findAtom(a1), res2.findAtom(a2))
        b.display = 1
        b.drawMode = Bond.Wire
        b.label = "||"
        self.chimeraBond = b
        self.updateRepresentationColor()

    ## Update the bond's color according to state of the interaction in the GUI.
    #  New constraint: blue
    #  Non validated constraint: orange
    #  Normal state: black
    def updateRepresentationColor(self):
        if self.chain2 == "Ground":
            poly = polymers[self.chain1]
            resSel = poly.getChimeraResidues(self.res1, self.res1)
            for r in resSel:
                r.labelColor = colorTable.getColorByName("black")
                r.ribbonColor = colorTable.getColorByName("black")
                for a in r.atoms:
                    a.color = colorTable.getColorByName("black")
            return

        if self.new:
            self.chimeraBond.color = colorTable.getColorByName("blue")
            return
        if not self.valid:
            self.chimeraBond.color = colorTable.getColorByName("orange")
            return
        
        self.chimeraBond.color = colorTable.getColorByName("black")

## Synchronize the constraints list with MMB
def refreshConstraints():
    global constraints
    mmbConstraints = pyMMB.getConstraints()
    constraints[:] = []
    for c in mmbConstraints:
        Constraint.mutateConstraint_wrapper(c, currentModel)
        constraints.append(c)

    refreshConstraintsBonds()

## Update base pair bonds display
def refreshConstraintsBonds():
    global constraintsBonds

    if constraintsBonds:
        constraintsBonds.deleteAll()
    else:
        constraintsBonds = misc.getPseudoBondGroup('constraints', modelID=OpenModels.Default, hidden=False)
        constraintsBonds.name = "MMB Constraints"
        constraintsBonds.lineWidth = 2
        constraintsBonds.color = colorTable.getColorByName("black")
    
    [p.updateRepresentation() for p in constraints]


##################################################################################################
## Defines an includeAllResiduesWithin command as in MMB
#
class AllResiduesWithin(pyMMB.AllResiduesWithin_wrapper, MMB_UI_Object):
    def __init__(self, model=None, mmbID=-1, chain="", residue=0, radius=1.0, valid=True, new=False,select=False):
        MMB_UI_Object.__init__(self, model, mmbID)
        self.mmbID = mmbID
        self.chain = chain
        self.residue = residue
        self.radius = radius
        self.initialization(model, new, valid, select)

    def initialization(self, model=currentModel, new=False, valid=True, select=False):
        self.model = model
        self.valid = valid
        self.new = new
        self.select = select

    @classmethod
    ## Mutate a MobilizerWithin_wrapper object's type to MobilizerWithin
    def mutateAllResiduesWithin_wrapper(cls, AllResiduesWithin_wrapper_object, model=currentModel, new=False, valid=True):
        AllResiduesWithin_wrapper_object.__class__ = cls 
        AllResiduesWithin_wrapper_object.initialization(model, new, valid)

    def getResNameId(self):
        if self.residue == 0:
            return ""
        poly = polymers[self.chain]
        return poly.getResNameId(self.residue)

    def getChimeraSelection(self):
        if self.residue != 0 and self.radius > 0:
            # resWithin = pyMMB.getResiduesWithin(self.chain, self.residue, self.radius)
            # return [self.model.findResidue( MolResId(x[0], x[1])) for x in resWithin]
            # chimera.selection.addCurrent(self.selection)
            return getResiduesWithin(self.chain, self.residue, self.radius*10.0)

    def updateRepresentation(self):
        resSel = self.getChimeraSelection()
        if resSel:
            for r in resSel:
                for a in r.atoms:
                    a.display = True
        
    def __str__(self):
        return "includeAllResiduesWithin %.1f %s %i" % (self.radius, self.chain, self.residue)

def refreshAllResiduesWithins():
    global allResiduesWithins

    allResiduesWithins[:] = []
    mmbAllResiduesWithins = pyMMB.getAllResiduesWithins()
    for a in mmbAllResiduesWithins:
        AllResiduesWithin.mutateAllResiduesWithin_wrapper(a, currentModel)
        allResiduesWithins.append(a)


##################################################################################################
## Defines an includeAllNonBondAtomsInResidue command as in MMB
#
class IncludeAllNonBondAtomsInResidue(pyMMB.IncludeAllNonBondAtomsInResidue_wrapper, MMB_UI_Object):
    def __init__(self, model=None, mmbID=-1, chain="", residue=0, valid=True, new=False,select=False):
        self.model = model
        self.mmbID = mmbID
        self.chain = chain
        self.residue = residue
        self.initialization(model, new, valid, select)

    def initialization(self, model=currentModel, new=False, valid=True, select=False):
        self.model = model
        self.valid = valid
        self.new = new
        self.select = select
        self.resEnd = self.residue

    @classmethod
    ## Mutate a MobilizerWithin_wrapper object's type to MobilizerWithin
    def mutateIncludeAllNonBondAtomsInResidue_wrapper(cls, IncludeAllNonBondAtomsInResidue_wrapper_object, model=currentModel, new=False, valid=True):
        IncludeAllNonBondAtomsInResidue_wrapper_object.__class__ = cls 
        IncludeAllNonBondAtomsInResidue_wrapper_object.initialization(model, new, valid)

    def setAttribute(self, attrName, value):
        MMB_UI_Object.setAttribute(self, attrName, value)
        if attrName == "residue" and self.residue > self.resEnd:
            self.resEnd = self.residue

    def getResStartNameId(self):
        if self.residue == 0:
            return ""
        poly = polymers[self.chain]
        return poly.getResNameId(self.residue)

    def getResEndNameId(self):
        if self.resEnd == 0:
            return ""
        poly = polymers[self.chain]
        return poly.getResNameId(self.resEnd)

    def getChimeraSelection(self):
        if self.chain != "":
            poly = polymers[self.chain]
            return poly.getChimeraResidues(self.residue, self.resEnd)

    def updateRepresentation(self):
        resSel = self.getChimeraSelection()
        if resSel:
            for r in resSel:
                for a in r.atoms:
                    a.display = True

    def __str__(self):
        return "includeAllNonBondAtomsInResidues %s %i %i" % (self.chain, self.residue, self.resEnd)

def refreshIncludeAllNonBondAtomsInResidues():
    global includeAllNonBondAtomsInResidues

    includeAllNonBondAtomsInResidues[:] = []
    mmbIncludeAllNonBondAtomsInResidues = pyMMB.getIncludeAllNonBondAtomsInResidues()
    for a in mmbIncludeAllNonBondAtomsInResidues:
        IncludeAllNonBondAtomsInResidue.mutateIncludeAllNonBondAtomsInResidue_wrapper(a, currentModel)
        includeAllNonBondAtomsInResidues.append(a)


##################################################################################################
## Defines a MMB's Contact stretch
#  Some methods are used as callbacks for Tk widgets.
#  Others are used to get information from Chimera
class Contact(pyMMB.ContactStretch_wrapper, MMB_UI_Object):

    typesList = ["AllAtomSterics", "AllHeavyAtomSterics", "SelectedAtoms"]

    ## Constructor
    def __init__(self, model=currentModel, mmbID=-1, contactScheme="AllAtomSterics", chain="", resStart=0, resEnd=-1, new=True):
        MMB_UI_Object.__init__(self, model, mmbID)
        self.contactScheme   = contactScheme
        self.chainID    = chain
        self.resStart   = resStart
        self.resEnd     = resEnd
        self.initialization(model, new)

    def initialization(self, model=currentModel, new=False, valid=True):
        self.model = model
        self.valid = valid
        self.new   = new

    @classmethod
    ## Mutate a ContactStretch_wrapper object's type to Contact
    def mutateContactStretch_wrapper(cls, ContactStretch_wrapper_object, model=currentModel, new=False, valid=True):
        ContactStretch_wrapper_object.__class__ = cls 
        ContactStretch_wrapper_object.initialization(model, new, valid)

    def setAttribute(self, attrName, value):
        if getattr(self, attrName) == value:
            return

        if attrName == "resStart" and value == "All":
            self.resStart = 0
            self.resEnd = -1
            return

        if attrName == "resStart":
            value = int(value.split()[1])
            self.resStart = value
            if value > self.resEnd:
                self.resEnd = value
            return

        if attrName.startswith("res"):
            value = int(value.split()[1])
        setattr(self, attrName, value)

    def getResStartNameId(self):
        if self.resStart == 0:
            return "All"
        poly = polymers[self.chainID]
        return poly.getResNameId(self.resStart)

    def getResEndNameId(self):
        if self.resEnd == -1:
            return ""
        poly = polymers[self.chainID]
        return poly.getResNameId(self.resEnd)

    ## Return a list of Chimera residue objects
    #  @return a list containing the residues.
    def getChimeraSelection(self):
        if self.chainID == "":
            return []
        poly = polymers[self.chainID]
        if self.resStart == "All":
            return poly.getChimeraResidues(start=poly.firstResNum)
        return poly.getChimeraResidues(self.resStart, self.resEnd)

# ## Synchronize the contacts list with MMB
def refreshContactStretches():
    global contacts

    mmbContacts = pyMMB.getContactStretches()
    contacts[:] = []
    for c in mmbContacts:
        Contact.mutateContactStretch_wrapper(c, currentModel, False)
        contacts.append(c)

##################################################################################################
## Defines a MMB's contactWithin
#  Some methods are used as callbacks for Tk widgets.
#  Others are used to get information from Chimera
class ContactWithin(pyMMB.ContactWithin_wrapper, MMB_UI_Object):
    ## Constructor
    def __init__(self, model=currentModel, mmbID=-1, contactScheme="AllAtomSterics", chain="", res=0, radius=0.4, new=True):
        MMB_UI_Object.__init__(self, model, mmbID)
        self.contactScheme   = contactScheme
        self.radius     = radius
        self.chainID    = chain
        self.resID      = res
        self.initialization(model, new)

    def initialization(self, model=currentModel, new=True, selection=None, valid=True):
        self.model      = model
        self.new        = new
        self.selection  = selection
        self.valid      = valid

    @classmethod
    ## Mutate a ContactWithin_wrapper object's type to ContactWithin
    def mutateContactWithin_wrapper(cls, ContactWithin_wrapper_object, model=currentModel, new=False, selection=None, valid=True):
        ContactWithin_wrapper_object.__class__ = cls 
        ContactWithin_wrapper_object.initialization(model, new, selection, valid)

    ## Initialize the residue when changing chain
    def setAttribute(self, attrName, value):
        if attrName.startswith("chain"):
            setattr(self, "resID", 0)
        MMB_UI_Object.setAttribute(self, attrName, value)

    def getResNameId(self):
        if self.resID == 0:
            return ""
        poly = polymers[self.chainID]
        return poly.getResNameId(self.resID)

    def getChimeraSelection(self):
        if self.resID != 0 and self.radius > 0:
            # resWithin = pyMMB.getResiduesWithin(self.chainID, self.resID, self.radius)
            # return [self.model.findResidue( MolResId(x[0], x[1])) for x in resWithin]
            return getResiduesWithin(self.chainID, self.resID, self.radius*10.0)

## Synchronize the contactsWithin list with MMB
def refreshContactWithin():
    global contactsWithin

    mmbContacts = pyMMB.getContactsWithin()
    contactsWithin[:] = []
    for m in mmbContacts:
        ContactWithin.mutateContactWithin_wrapper(m, currentModel)
        contactsWithin.append(m)


def selectionIntersection(sel1, sel2):
    sel1_IDs = set([(a.residue.id.position,a.name) for a in sel1])
    sel2_IDs = set([(a.residue.id.position,a.name) for a in sel2])
    intersect = sel1_IDs.intersection(sel2_IDs)
    sel1_final = [a for a in sel1 if (a.residue.id.position,a.name) in intersect]
    sel2_final = [a for a in sel2 if (a.residue.id.position,a.name) in intersect]
    return sel1_final, sel2_final

##################################################################################################
## Defines a MMB's threading
#  Some methods are used as callbacks for Tk widgets.
#  Others are used to get information from Chimera
class Threading(pyMMB.ThreadingStruct_wrapper, MMB_UI_Object):
    ## Constructor
    def __init__(self, model=currentModel, mmbID=-1, chain1="", resStart1=0, resEnd1=0, 
                                                     chain2="", resStart2=0, resEnd2=0, 
                                                     forceConstant=3.0, backboneOnly=False, new=True):
        MMB_UI_Object.__init__(self, model, mmbID)
        self.chainID1      = chain1
        self.residueStart1 = resStart1
        self.residueEnd1   = resEnd1
        self.chainID2      = chain2
        self.residueStart2 = resStart2
        self.residueEnd2   = resEnd2
        self.forceConstant = forceConstant
        self.backboneOnly  = backboneOnly
        self.initialization(model, new)

    def initialization(self, model=currentModel, new=True, selection=None, valid=True):
        self.model      = model
        self.new        = new
        self.selection  = selection
        self.valid      = valid

    @classmethod
    ## Mutate a ThreadingStruct_wrapper object's type to Threading
    def mutateThreadingStruct_wrapper(cls, ThreadingStruct_wrapper_object, model=currentModel, new=False, selection=None, valid=True):
        ThreadingStruct_wrapper_object.__class__ = cls 
        ThreadingStruct_wrapper_object.initialization(model, new, selection, valid)

    ## Initialize the residue when changing chain
    def setAttribute(self, attrName, value):
        if attrName == "chain1":
            setattr(self, "residueStart1", 0)
            setattr(self, "residueStart2", 0)
        elif attrName == "chain2":
            setattr(self, "residueStart1", 0)
            setattr(self, "residueStart2", 0)
        MMB_UI_Object.setAttribute(self, attrName, value)

    def getResStart1NameId(self):
        if self.residueStart1 == 0:
            return ""
        poly = polymers[self.chainID1]
        return poly.getResNameId(self.residueStart1)
    def getResStart2NameId(self):
        if self.residueStart2 == 0:
            return ""
        poly = polymers[self.chainID2]
        return poly.getResNameId(self.residueStart2)

    def getResEnd1NameId(self):
        if self.residueEnd1 == 0:
            return ""
        poly = polymers[self.chainID1]
        return poly.getResNameId(self.residueEnd1)
    def getResEnd2NameId(self):
        if self.residueEnd2 == 0:
            return ""
        poly = polymers[self.chainID2]
        return poly.getResNameId(self.residueEnd2)

    def getChimeraSelection1(self):
        sel = []
        if self.chainID1:
            poly = polymers[self.chainID1]
            sel += poly.getChimeraResidues(self.residueStart1, self.residueEnd1)
        return sel
    def getChimeraSelection2(self):
        sel = []
        if self.chainID2:
            poly = polymers[self.chainID2]
            sel += poly.getChimeraResidues(self.residueStart2, self.residueEnd2)
        return sel

    def getChimeraSelection(self):
        if self.chainID1 == ""  and self.chainID2 == "":
            return []
        sel = []
        sel += self.getChimeraSelection1()
        sel += self.getChimeraSelection2()
        return sel

    def matchSelections(self):
        selRes1 = self.getChimeraSelection1()
        selRes2 = self.getChimeraSelection2()
        sel1 = []
        [sel1.extend(r.atoms) for r in selRes1]
        sel2 = []
        [sel2.extend(r.atoms) for r in selRes2]
        inter1, inter2 = selectionIntersection(sel1, sel2)
        Midas.match(inter1, inter2, move="chains")

## Synchronize the mobilizersWithin list with MMB
def refreshThreadings():
    global threadings

    mmbThreadings = pyMMB.getThreadingStructs()
    threadings[:] = []
    for m in mmbThreadings:
        Threading.mutateThreadingStruct_wrapper(m, currentModel)
        threadings.append(m)

##################################################################################################
## Defines a MMB's GappedThreading
#  Some methods are used as callbacks for Tk widgets.
#  Others are used to get information from Chimera
class GappedThreading(pyMMB.GappedThreadingStruct_wrapper, MMB_UI_Object):
    ## Constructor
    def __init__(self, model=currentModel, mmbID=-1, chain1="",
                                                     chain2="", 
                                                     forceConstant=3.0, backboneOnly=False, new=True):
        MMB_UI_Object.__init__(self, model, mmbID)
        self.chainID1      = chain1
        self.chainID2      = chain2
        self.forceConstant = forceConstant
        self.backboneOnly  = backboneOnly
        self.initialization(model, new)

    def initialization(self, model=currentModel, new=True, selection=None, valid=True):
        self.model      = model
        self.new        = new
        self.selection  = selection
        self.valid      = valid

    @classmethod
    ## Mutate a GappedThreadingStruct_wrapper object's type to Threading
    def mutateGappedThreadingStruct_wrapper(cls, GappedThreadingStruct_wrapper_object, model=currentModel, new=False, selection=None, valid=True):
        GappedThreadingStruct_wrapper_object.__class__ = cls 
        GappedThreadingStruct_wrapper_object.initialization(model, new, selection, valid)

    ## Initialize the residue when changing chain
    def setAttribute(self, attrName, value):
        MMB_UI_Object.setAttribute(self, attrName, value)

    def getChimeraSelection1(self):
        sel = []
        if self.chainID1:
            poly = polymers[self.chainID1]
            sel += poly.getChimeraResidues()
        return sel
    def getChimeraSelection2(self):
        sel = []
        if self.chainID2:
            poly = polymers[self.chainID2]
            sel += poly.getChimeraResidues()
        return sel

    def getChimeraSelection(self):
        if self.chainID1 == ""  and self.chainID2 == "":
            return []
        sel = []
        sel += self.getChimeraSelection1()
        sel += self.getChimeraSelection2()
        return sel

    def matchSelections(self):
        selRes1 = self.getChimeraSelection1()
        selRes2 = self.getChimeraSelection2()
        sel1 = []
        [sel1.extend(r.atoms) for r in selRes1]
        sel2 = []
        [sel2.extend(r.atoms) for r in selRes2]
        inter1, inter2 = selectionIntersection(sel1, sel2)
        Midas.match(inter1, inter2, move="chains")

## Synchronize the mobilizersWithin list with MMB
def refreshGappedThreadings():
    global gappedThreadings

    mmbThreadings = pyMMB.getGappedThreadingStructs()
    gappedThreadings[:] = []
    for m in mmbThreadings:
        GappedThreading.mutateGappedThreadingStruct_wrapper(m, currentModel)
        gappedThreadings.append(m)

##################################################################################################
## Defines a MMB's AtomSpring
#  Some methods are used as callbacks for Tk widgets.
#  Others are used to get information from Chimera
class AtomSpring(pyMMB.AtomSpring_wrapper, MMB_UI_Object):
    ## Constructor
    def __init__(self, model=currentModel, mmbID=-1, chain1="",
                                                     res1 = 0,
                                                     atom1 = "",
                                                     chain2="",
                                                     res2 = 0,
                                                     atom2 = "",
                                                     tether = False, 
                                                     toGround = False, 
                                                     deadLength = 0.0,
                                                     forceConstant=3.0, 
                                                     new=True):

        MMB_UI_Object.__init__(self, model, mmbID)
        self.atom1          = atom1
        self.atom2          = atom2
        self.res1           = res1
        self.res2           = res2
        self.chain1         = chain1
        self.chain2         = chain2
        self.toGround     = toGround
        self.tether       = tether
        self.forceConstant= forceConstant
        self.deadLength   = deadLength
        self.initialization(model, new)

    def initialization(self, model=currentModel, new=True, selection=None, valid=True):
        self.model      = model
        self.new        = new
        self.selection  = selection
        self.valid      = valid
        self.chimeraBond = None

    @classmethod
    ## Mutate a AtomSpring_wrapper object's type to AtomSpring
    def mutateAtomSpring_wrapper(cls, AtomSpring_wrapper_object, model=currentModel, new=False, selection=None, valid=True):
        AtomSpring_wrapper_object.__class__ = cls 
        AtomSpring_wrapper_object.initialization(model, new, selection, valid)


    def setAttribute(self, attrName, value):
        if getattr(self, attrName) == value:
            return
        # if the attribute is res1 or res2 we set the corresponding atom to ""
        if "res" in attrName:
            nAttr = attrName[-1]
            atomAttr = "atom"+nAttr
            setattr(self, atomAttr, "")
            if value:
                value = int(value.split()[1])
            else:
                value = 0
        # if the attribute is chain1 or chain2 we set the corresponding res,atom to 0,""
        elif "Chain" in attrName:
            nAttr = attrName.split("n")[1]
            setattr(self, "res"+nAttr, 0)
            setattr(self, "atom"+nAttr, "")
        setattr(self, attrName, value)

    def getRes1AtomsList(self):
        return getAtomsList(self.chain1, self.res1, self.model)

    def getRes2AtomsList(self):
        return getAtomsList(self.chain2, self.res2, self.model)

    def getRes1NameId(self):
        if self.res1 == 0:
            return ""
        return polymers[self.chain1].getResNameId(self.res1)

    def getRes2NameId(self):
        if self.res2 == 0:
            return ""
        return polymers[self.chain2].getResNameId(self.res2)
        
    def getChimeraSelection(self):
        if self.chain1 == 0 or self.res1 == 0:
            return

        res1 = self.model.findResidue(MolResId(self.chain1, self.res1))
        selection = [self.chimeraBond]
        if self.atom1 != "":
            selection.append(res1.findAtom(self.atom1))
        else:
            selection.append(res1)
        if self.res2 == 0:
            return selection
        
        res2 = self.model.findResidue(MolResId(self.chain2, self.res2))
        if self.atom2 != "":
            selection.append(res2.findAtom(self.atom2))
        else:
            selection.append(res2)
        return selection

    ## Update the bond display in Chimera. Create a new bond if necessary.
    def updateRepresentation(self):
        global atomSpringsBonds
        # No bond
        if self.res2 == 0:
            if self.chimeraBond: 
                atomSpringsBonds.deletePseudoBond(self.chimeraBond)
            self.chimeraBond = None
            return

        res1 = self.model.findResidue(MolResId(self.chain1, self.res1))
        res2 = self.model.findResidue(MolResId(self.chain2, self.res2))
        a1 = self.atom1
        a2 = self.atom2
        if self.atom1 == "":
            a1 = polymers[self.chain1].getRepresentativeAtom()
        if self.atom2 == "":
            a2 = polymers[self.chain2].getRepresentativeAtom()

        b = self.chimeraBond
        if b:
            b.reuse(res1.findAtom(a1), res2.findAtom(a2))
            # b.label = "||"
            self.updateRepresentationColor()
            return

        b = atomSpringsBonds.newPseudoBond(res1.findAtom(a1), res2.findAtom(a2))
        b.display = 1
        b.drawMode = Bond.Spring
        # b.label = "||"
        self.chimeraBond = b
        self.updateRepresentationColor()

    ## Update the bond's color according to state of the interaction in the GUI.
    #  New atomspring: blue
    #  Non validated atomspring: orange
    #  Normal state: black
    def updateRepresentationColor(self):
        if self.new:
            self.chimeraBond.color = colorTable.getColorByName("blue")
            return
        if not self.valid:
            self.chimeraBond.color = colorTable.getColorByName("orange")
            return
        
        self.chimeraBond.color = colorTable.getColorByName("black")

## Synchronize the atomsprings list with MMB
def refreshAtomSprings():
    global atomSprings
    mmbAtomsprings = pyMMB.getAtomSprings()
    atomSprings[:] = []
    for c in mmbAtomsprings:
        AtomSpring.mutateAtomSpring_wrapper(c, currentModel)
        atomSprings.append(c)

    refreshAtomSpringsBonds()

## Update base pair bonds display
def refreshAtomSpringsBonds():
    global atomSpringsBonds

    if atomSpringsBonds:
        atomSpringsBonds.deleteAll()
    else:
        atomSpringsBonds = misc.getPseudoBondGroup('atomsprings', modelID=OpenModels.Default, hidden=False)
        atomSpringsBonds.name = "MMB AtomSprings"
        atomSpringsBonds.lineWidth = 2
        atomSpringsBonds.color = colorTable.getColorByName("black")
    
    [p.updateRepresentation() for p in atomSprings]


##################################################################################################
class RootMobilizer(MMB_UI_Object):

    rootMobilizerTypes = ["Free", "Weld"]

    ## Constructor
    def __init__(self, model=currentModel, mmbID=-1, chainID="", rootMobilizer=""):
        MMB_UI_Object.__init__(self, model, mmbID)
        self.chainID = chainID
        self.rootMobilizer = rootMobilizer
        self.new = False
        self.valid = True

    def __str__(self):
        return "rootMobilizer " + self.chainID + " " + self.rootMobilizer

    def mmbUpdate(self):
        pyMMB.setRootMobilizer(self.chainID, self.rootMobilizer)

    def getChimeraSelection(self):
        if self.chainID == "":
            return []
        sel = []
        poly = polymers[self.chainID]
        sel += poly.getChimeraResidues()
        return sel

    def updateRepresentation(self):
        if self.chainID:
            r = self.getChimeraSelection()[0]
            r.label = u"\u22A5"
            self.updateRepresentationColor()

    def updateRepresentationColor(self):
        if self.chainID:
            r = self.getChimeraSelection()[0]
            r.labelColor = colorTable.getColorByName("black")
            r.ribbonColor = colorTable.getColorByName("black")
            for a in r.atoms:
                a.color = colorTable.getColorByName("black")

def refreshRootMobilizers():
    global rootMobilizers

    rootMobilizers[:] = []
    for p in polymers.values():
        mobi = p.getRootMobilizer()
        rootMobilizers.append(RootMobilizer(chainID=p.chainID, rootMobilizer=mobi))


##################################################################################################
## Defines a MMB's Density stretch
#  Some methods are used as callbacks for Tk widgets.
#  Others are used to get information from Chimera
class Density(pyMMB.DensityStretch_wrapper, MMB_UI_Object):

    ## Constructor
    def __init__(self, model=currentModel, mmbID=-1, chain="", resStart=0, resEnd=-1, new=True):
        MMB_UI_Object.__init__(self, model, mmbID)
        self.chainID    = chain
        self.resStart   = resStart
        self.resEnd     = resEnd
        self.initialization(model, new)

    def initialization(self, model=currentModel, new=False, valid=True):
        self.model = model
        self.valid = valid
        self.new   = new

    @classmethod
    ## Mutate a DensityStretch_wrapper object's type to Density
    def mutateDensityStretch_wrapper(cls, DensityStretch_wrapper_object, model=currentModel, new=False, valid=True):
        DensityStretch_wrapper_object.__class__ = cls 
        DensityStretch_wrapper_object.initialization(model, new, valid)

    def setAttribute(self, attrName, value):
        if getattr(self, attrName) == value:
            return

        if attrName == "resStart" and value == "All":
            self.resStart = 0
            self.resEnd = -1
            return

        if attrName == "resStart":
            value = int(value.split()[1])
            self.resStart = value
            if value > self.resEnd:
                self.resEnd = value
            return

        if attrName.startswith("res"):
            value = int(value.split()[1])
        setattr(self, attrName, value)

    def getResStartNameId(self):
        if self.resStart == 0:
            return "All"
        poly = polymers[self.chainID]
        return poly.getResNameId(self.resStart)

    def getResEndNameId(self):
        if self.resEnd == -1:
            return ""
        poly = polymers[self.chainID]
        return poly.getResNameId(self.resEnd)

    ## Return a list of Chimera residue objects
    #  @return a list containing the residues.
    def getChimeraSelection(self):
        if self.chainID == "":
            return []
        poly = polymers[self.chainID]
        if self.resStart == "All":
            return poly.getChimeraResidues(start=poly.firstResNum)
        return poly.getChimeraResidues(self.resStart, self.resEnd)

# ## Synchronize the densitys list with MMB
def refreshDensityStretches():
    global densities

    mmbDensities = pyMMB.getDensityStretches()
    densities[:] = []
    for c in mmbDensities:
        Density.mutateDensityStretch_wrapper(c, currentModel, False)
        densities.append(c)

##################################################################################################
## Open a density map (Xplor format only) and set MMB's parameters accordingly
def openDensityMap(filepath):
    global densityMapModel

    if densityMapModel and filepath == MMBparameters.densityFileName:
        return
    try:
        newModel = chimera.openModels.open(filepath)[0]
    except Exception as e:
        return

    if densityMapModel:
        densityMapModel.close()
    densityMapModel = newModel
    MMBparameters.densityFileName = filepath


##################################################################################################
## Display molecules with MMB information
def applyMMBVisualization():
    from chimera import runCommand as rc

    rc("~display")
    rc("color none")
    rc("~rlabel")
    rc("ribspline bspline")
    rc("ribbackbone")
    rc("ribscale 'Chimera default'")
    for array in allConstraintsRestraints.values():
        for command in array:
            command.updateRepresentation()
    rc("window")


##################################################################################################
## Commands to string
def commandsToStr():
    buff = ""

    buff += "loadSequencesFromPdb mmb.pdb"

    if (allResiduesWithins + includeAllNonBondAtomsInResidues) or getMMBParameter("physicsRadius")>0:
        buff += ("setDefaultMDParameters\n\n")

    for p in polymers.values():
        if p.activePhysics == False:
            buff += "deactivatePhysics "+p.chainID+"\n"
    buff += "\n"

    for array in allConstraintsRestraints.values():
        for command in array:
            buff += (str(command)+"\n")
        buff += ("\n")

    # print MMBparameters.nonDefaultParameters
    for p in MMBparameters.nonDefaultParameters:
        if p != "converged":
            buff += p + " " + str(getMMBParameter(p)) + "\n"

    return buff

## Dump commands
def dumpCommands(fileName):
    buff = commandsToStr()

    fileOut = open(fileName, "w")
    fileOut.write(buff)
    fileOut.close()

    return buff


# def resetMMB(keepCurrentFrame=False):
#     if keepCurrentFrame:
        
        
##################################################################################################
# Bindings with pyMMB



## Execute a set of MMB commands.
#  @param a list of MMB commands lines
def sendMMBCmds(cmdsLines):
    ignore = ["RNA","DNA","protein"]
    ignored = []
    notSupported = ["read"]
    for l in cmdsLines:
        l = l.strip()
        if not l:
            continue
        tmp = [l.startswith(s) for s in ignore]
        if True in tmp:
            print "Ignored the following command:"
            print l
            ignored.append(l)
            continue
        tmp = [l.startswith(s) for s in notSupported]
        if True in tmp:
            return ">>> "+l.split()[0]
        if l.startswith("loadSequencesFromPdb"):
            if polymersInitialized:
                print "Ignored the following command:"
                print l
                ignored.append(l)
                continue 
            cmd = l.split()
            if ( len(cmd) == 1 ):
                fileName = "last.0.pdb"
            else:
                fileName = cmd[1]
            fileName = workDir+'/'+fileName
            l = 'loadSequencesFromPdb '+fileName
        pyMMB.cmd(l)
    return ignored

## Add an MMB object to MMB using its MMB command returned by str(mmbObj)
#  @param mmbObj the object to add
def addMMBObject(mmbObj):
    pyMMB.cmd(str(mmbObj))

## Clear MMB forces and constraints
def clearForcesAndConstraints():
    pyMMB.clearForcesAndConstraints()


def setDefaultMDParameters():
    pyMMB.cmd("setDefaultMDParameters")

##################################################################################################
# MMB simulation control
import timeit
## Asynchronous control of MMB simulation
class MMBSimulationThread(threading.Thread):
    test = 0
    ## Constructor
    #  @param callback a function called after of each frame
    #  @param maxFrame maximum number of frame to compute
    def __init__(self, maxFrame = 0, endCallback = None):
        threading.Thread.__init__(self)
        self.running = threading.Event()
        self.running.set()

        self.callbacks = []

        if endCallback and not hasattr(endCallback, '__call__'):
            raise MMBError("MMBSimulationThread::addCallback", "Argument is not callable")
        self.endCallback = endCallback

        self.maxFrame = maxFrame
        self.currentFrame = 0
        pyMMB.MMBparameters.converged = False
        MMBSimulationThread.test += 1
        print "Ready to run"

    def run(self):
        global runFlag
        runFlag = True
        t0 = time.clock()
        print "HELLO:", MMBSimulationThread.test
        while self.currentFrame < self.maxFrame and not pyMMB.MMBparameters.converged:
            print "Running"
            pyMMB.runOneStep()
            self.currentFrame += 1
            for cb in self.callbacks:
                cb()          
            self.running.wait()
        runFlag = False
        print time.clock() - t0
        self.endCallback()

    def addCallback(self, callback):
        if callback in self.callbacks:
            return
        if not hasattr(callback, '__call__'):
            raise MMBError("MMBSimulationThread::addCallback", "Argument is not callable")
        self.callbacks.append(callback)

    def pause(self):
        self.running.clear()

    def restart(self):
        self.running.set()

    def stop(self):
        self.currentFrame = self.maxFrame

def addLastFrameFromMMB(movieDialog):
    global currentModel
    pos = pyMMB.getSystemCoordinates()
    # print len(pos)
    movieDialog.addFrame(pos)
    # cs = currentModel.newCoordSet(len(currentModel.coordSets))
    # chimera.fillCoordSet(cs, currentModel.atoms, pos)
    # currentModel.activeCoordSet = cs
    # movieDialog.updateFrames()

## Init MMB simulation and display the MDMovie dialog
def initSimulation():
    global polymers
    for p in polymers.values():
        p.matchCoordinatesFromChimera()
    pyMMB.cmd("numReportingIntervals 1000")
    pyMMB.postInitialize()
    pyMMB.initializeBodies()
    pyMMB.initDynamics()

    global currentModel
    cs = currentModel.activeCoordSet
    chimera.fillCoordSet(cs, currentModel.atoms, pyMMB.getSystemCoordinates())
    chimera.runCommand("window #"+str(currentModel.id))

    global ensemble
    if not ensemble:
        ensemble = MolEnsemble(currentModel)

    global movieDialog
    if not movieDialog:
        movieDialog = MyMovieDialog(ensemble, externalEnsemble=True)

    global initFlag 
    initFlag = True

## Run n intervals of the simulation.
#  Each new frame is added to the MDMovie dialog
#  The callback function is called after each frame
def runIntervals(n=1, callback=None, endCallback=None):
    """
    n is the number of intervals to compute
    """
    global currentSimulation
    global movieDialog

    if currentSimulation and currentSimulation.isAlive():
        raise pyMMB.MMBError("MMB_UI.runIntervals", "A simulation is currently running.")
    currentSimulation = MMBSimulationThread(maxFrame=n, endCallback=endCallback)
    currentSimulation.addCallback(lambda md=movieDialog: addLastFrameFromMMB(md))
    currentSimulation.addCallback(callback)
    currentSimulation.start()

## Pause a simulation
def pauseSimulation():
    global currentSimulation
    if currentSimulation and currentSimulation.isAlive():
        currentSimulation.pause()
    else:
        raise pyMMB.MMBError("MMB_UI.pauseSimulation", "No simulation is currently running.")

## Continue the current paused simulation
def continueSimulation():
    global currentSimulation
    if currentSimulation and not currentSimulation.isAlive():
        raise pyMMB.MMBError("MMB_UI.continueSimulation", "No simulation is currently running.")
    currentSimulation.restart()

## Stop a simulation
def stopSimulation():
    global currentSimulation
    if currentSimulation and currentSimulation.isAlive():
        currentSimulation.stop()
    else:
        raise pyMMB.MMBError("MMB_UI.stopSimulation", "No simulation is currently running.")

## Return the number of satisfied and unsatisfied base pairs in a tuple
#  @return a tuple with the numbers
def getSatisfiedBasePairs():
    oks = pyMMB.getNumSatisfiedBasePairs()
    notOks = pyMMB.getNumUnSatisfiedBasePairs()
    return (oks, notOks)

def  getCurrentFrameNumber():
    if currentSimulation:
        return currentSimulation.currentFrame
    return 0




