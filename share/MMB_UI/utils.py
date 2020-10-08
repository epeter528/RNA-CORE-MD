import Tkinter as Tk
import MMB_UI
import pyMMB
import tkMessageBox
import Tix
import ttk
import string

# import pdb

#############################################################################
class ComboBox(Tk.OptionMenu):
    """
    A combobox with a self managed traced variable
    """
    def __init__(self, master, options, default, data=None, command=None, *args, **kw):
        self.var = Tk.StringVar(master)
        self.var.set(default)
        self.data = data
        if not command:
            self.var.trace("w", lambda name, index, mode: self.defaultCommand())
        else:
            self.var.trace("w", lambda name, index, mode, x=self: command(x))
        params = list(options) + list(args)
        Tk.OptionMenu.__init__(self, master, self.var, *params, **kw)

    def get(self):
        return self.var.get()

    def set(self, val):
        self.var.set(val)

    def defaultCommand(self):
        pass

#############################################################################
class EmptyLabel(Tk.Label):
    def __init__(self, parent, *args, **kw):
        Tk.Label.__init__(self, parent, text="", *args, **kw)

#############################################################################
class ValidatingEntry(Tk.Entry):
    def __init__(self, parent, value=None, command=None, data=None, extVariable=None, *args, **kw):
        # extVariable and command,value are exclusive. If you want to use an external variable, set it up yourself

        # valid percent substitutions (from the Tk entry man page)
        # %d = Type of action (1=insert, 0=delete, -1 for others)
        # %i = index of char string to be inserted/deleted, or -1
        # %P = value of the entry if the edit is allowed
        # %s = value of entry prior to editing
        # %S = the text string being inserted or deleted, if any
        # %v = the type of validation that is currently set
        # %V = the type of validation that triggered the callback
        #      (key, focusin, focusout, forced)
        # %W = the tk name of the widget
        vcmd = (parent.register(self.OnValidate), 
                        '%d', '%i', '%P', '%s', '%S', '%v', '%V', '%W')
        self.command = command or self.dummyCmd
        self.data = data
        if extVariable:
            self.textVar = extVariable
        else:
            self.textVar = Tk.StringVar()
            if value: self.textVar.set(str(value))
            self.textVar.trace("w", lambda name, index, mode, x=self: self.command(x))
        Tk.Entry.__init__(self, parent, validate="key", 
                          validatecommand=vcmd, 
                          textvariable=self.textVar, 
                          *args, **kw)
        self.balloon = Tix.Balloon()

    def dummyCmd(self, *args):
        pass

    def get(self):
        return Tk.Entry.get(self)

    def set(self, val):
        self.textVar.set(str(val))

    def validate(self):
        self.configure(background="white")

    def invalidate(self, message):
        self.configure(background="red")
        # self.balloon.bind_widget(self, msg=message)

    def OnValidate(self, action, index, value_if_allowed,
                   prior_value, text, validation_type, trigger_type, widget_name):
        return True


#############################################################################
class NumericEntry(ValidatingEntry):
    def get(self):
        try:
            return float(Tk.Entry.get(self))
        except ValueError:
            return 0.0

    def OnValidate(self, action, index, value_if_allowed,
                   prior_value, text, validation_type, trigger_type, widget_name):
        if value_if_allowed == "":
            return True
        try:
            float(value_if_allowed)
            return True
        except ValueError:
            return False
        return False

#############################################################################
class IntegerEntry(ValidatingEntry):

    def get(self):
        try:
            return int(Tk.Entry.get(self))
        except ValueError:
            return 0

    def OnValidate(self, action, index, value_if_allowed,
                   prior_value, text, validation_type, trigger_type, widget_name):
        if value_if_allowed == "":
            return True
        try:
            int(value_if_allowed)
            return True
        except ValueError:
            return False
        return False

#############################################################################
class SequenceEntry(ValidatingEntry):
    def get(self):
        return Tk.Entry.get(self).upper()

    def OnValidate(self, action, index, value_if_allowed,
                   prior_value, text, validation_type, trigger_type, widget_name):
        if action == "0":
            return True
        if not text.isalpha():
            return False
        # pdb.set_trace()
        return True

#############################################################################
class ChainIDEntry(ValidatingEntry):
    def get(self):
        return Tk.Entry.get(self).upper()
    def OnValidate(self, action, index, value_if_allowed,
                   prior_value, text, validation_type, trigger_type, widget_name):
        if action == "0":
            return True
        if text not in string.ascii_letters+string.digits:
            return False
        if len(value_if_allowed) > 1:
            return False
        return True

#############################################################################
class TracedCheckButton(Tk.Checkbutton):
    def __init__(self, parent, text="", value=False, command=None, data=None, extVariable=None, *args, **kw):
        # extVariable and command,value are exclusive. If you want to use an external variable, set it up yourself
        cmd = None
        if extVariable:
            self.varBool = extVariable
        else:
            self.varBool = Tk.BooleanVar()
            self.varBool.set(value)
            cmd = self.dummyCmd
            if command:
                cmd = lambda wid=self: command(wid)
        self.data = data
        Tk.Checkbutton.__init__(self, parent, text=text, 
                                variable=self.varBool, 
                                command=cmd)
        self.bind("<Button-1>", self.giveMeFocus)

    def dummyCmd(self, *args):
        pass

    def get(self):
        return self.varBool.get()

    def set(self, value):
        if value:
            self.select()
        else:
            self.deselect()

    def validate(self):
        pass

    def invalidate(self):
        self.set(self.varBool.get())

    def giveMeFocus(self, event):
        self.focus_set()


#############################################################################
class LabeledWidget(Tk.Frame):
    def __init__(self, parent, entryclass=ValidatingEntry, labelText="", value=None, command=None, data=None, *args, **kw):
        Tk.Frame.__init__(self,parent)
        self.label = Tk.Label(self, text=labelText)
        self.label.pack(side="left")
        self.entry = entryclass(self, command=command,value=value,data=data, *args, **kw)
        self.entry.pack(side="left")
        if hasattr(parent, "tkTags"):
            self.bindtags(self.bindtags()+parent.tkTags)
            self.label.bindtags(self.label.bindtags()+parent.tkTags)
            self.entry.bindtags(self.entry.bindtags()+parent.tkTags)

    def get(self):
        return self.entry.get()
    def set(self, val):
        return self.entry.set(val)
    def config(self, cnf=None, **kw):
        try:
            return Tk.Frame.config(self, cnf, **kw)
        except Tk.TclError:
            pass
        return self.entry.config(cnf, **kw)
    def configure(self, cnf=None, **kw):
        self.config(cnf, **kw)


#############################################################################
## Management of MMB parameters variables by Tk Variables
class ParametersTraces():
    entryTypes = {  int:  Tk.IntVar,
                    float:Tk.DoubleVar,
                    bool: Tk.BooleanVar,
                    str:  Tk.StringVar
                 }
    traces = {}

    @classmethod
    def getVariable(cls, parameterName):
        paramType = type(MMB_UI.getMMBParameter(parameterName))
        var = ParametersTraces.traces.get(parameterName)
        if not var:
            var = ParametersTraces.traces.setdefault(parameterName, ParametersTraces.entryTypes[paramType]())
        # Make sure the trace is set to the freshest value from MMB
        value = MMB_UI.getMMBParameter(parameterName)
        var.set(value)
        return var

    @classmethod
    def get(cls, parameterName):
        var = ParametersTraces.getVariable(parameterName)
        return var.get()

    @classmethod
    def set(cls, parameterName, value):
        var = ParametersTraces.traces.get(parameterName)
        var.set(value)

#############################################################################
## Generation of LabeledWidget widgets for MMB parameters
class ParameterEntry(LabeledWidget):
    entryTypes = {  int:  IntegerEntry,
                    float:NumericEntry,
                    bool: TracedCheckButton,
                    str:  ValidatingEntry
                 }

    def __init__(self, parent, parameterName, guitext="", command=None, width=5):
        if not guitext:
            guitext = parameterName
        value = MMB_UI.getMMBParameter(parameterName)
        paramType = type(value)
        var = ParametersTraces.getVariable(parameterName)
        LabeledWidget.__init__(self, parent, 
                               ParameterEntry.entryTypes[paramType], 
                               labelText=guitext,
                               extVariable=var,
                               data=parameterName, 
                               width=width
                               )
        self.traceName = var.trace_variable("w", lambda name, index, mode, x=self.entry: ParameterEntry.changeParameter(x))

    ## We have to delete the callback when the widget is destroyed
    def destroy(self):
        var = ParametersTraces.getVariable(self.entry.data)
        var.trace_vdelete("w", self.traceName)
        LabeledWidget.destroy(self)


    @classmethod
    def changeParameter(cls, widget):
        try:
            if widget.master.focus_get() == widget:
                MMB_UI.setMMBParameter(widget.data, widget.get())
                widget.validate()
        except pyMMB.MMBError as e:
            widget.invalidate()
            tkMessageBox.showerror("MMB Error", e.msg)

#############################################################################
from chimera.baseDialog import ModelessDialog
## Dialog containing widgets for all MMB parameters
class AllParametersDialog(ModelessDialog):
    name = "MMB - Parameters"
    buttons = ("Close")

    title = "MMB - Parameters"

    ## GUI code
    def fillInUI(self, parent):
        paramNames = sorted( [x[0] for x in MMB_UI.MMBparameters._fields_] )
        
        self.fr = ScrolledFrame(parent, tkTags=("AllParametersFrameTag",))
        # self.fr.bindtags(self.fr.bindtags()+self.fr.tkTags)
        for i,p in enumerate(paramNames):
            w = ParameterEntry(self.fr, p, width=15)
            w.grid(row=i, column=0, sticky=Tk.S+Tk.N+Tk.E)



#############################################################################
## Defines a composition of widgets displayed on the same grid line
#  Widgets variables names must begin by 'wid'
class WidgetsLine:
    def __init__(self, parent, guiID, data, select=False):
        self.listFrame = parent
        self.guiID = guiID
        self.data = data

        self.varSelect = Tk.BooleanVar()
        self.varSelect.set(select)

        self.widValid = None
        self.widNew = None
        self.widMatch = None

        self.balloon = Tix.Balloon()

    ## Generator yielding +1 from 0. Useful for grid column option.
    def iter1(self):
        n = 0
        while 1:
            yield n
            n += 1

    def changeAttribute(self, x, redisplay=True):
        """
        x.data is the name of the attribute to update
        x.get() return the new value of this attribute
        """
        # self.data.removeRepresentation()
        self.data.removeFromSelection()
        attr_key = x.data
        attr = x.get()
        self.data.setAttribute(attr_key, attr)
        self.data.valid = False
        self.select(True)
        self.display()
        self.data.updateRepresentation()

    def applyTags(self):
        if hasattr(self.listFrame,"tkTags"):
            [w.bindtags(w.bindtags()+self.listFrame.tkTags) for w in [self.__dict__[k] for k in self.__dict__ if k.startswith("wid")] if w]

    def select(self, val):
        self.varSelect.set(val)
        self.selectLine()

    def selectLine(self):
        WidgetsLine.updateBgColor(self)
        if self.varSelect.get():
            self.data.addToSelection()
        else:
            self.data.removeFromSelection()

    def updateBgColor(self):
        # determine background color
        col = "snow"
        if self.varSelect.get():
            col = "pale green"
        if not self.data.valid:
            col = "light salmon"
        if self.data.new:
            col = "light blue"

        [w.config(bg=col) for w in [self.__dict__[k] for k in self.__dict__ if k.startswith("wid")] if w]

        if self.widNew:
            self.widNew.config(highlightbackground=col)
        
        if self.widMatch:
            self.widMatch.config(highlightbackground=col)

        if self.widValid:
            self.widValid.config(highlightbackground=col)

    def drawValidNewButtons(self, colID):
        if self.data.new:
            self.widNew = Tk.Button(self.listFrame, text="Add", command=self.addCommand)
            self.balloon.bind_widget(self.widNew, msg="Validate and add this mobilizer to MMB")
            self.widValid = EmptyLabel(self.listFrame)
            self.widNew.grid(row=self.guiID, column=colID,sticky=Tk.W+Tk.E+Tk.S+Tk.N) 
        elif not self.data.valid:
            self.widValid = Tk.Button(self.listFrame, text=u"\u2713", command=self.updateMMB)
            self.balloon.bind_widget(self.widValid, msg="Validate and commit changes to MMB")
            self.widNew = EmptyLabel(self.listFrame)
            self.widValid.grid(row=self.guiID, column=colID,sticky=Tk.W+Tk.E+Tk.S+Tk.N)
        else:
            self.widNew = EmptyLabel(self.listFrame)
            self.widValid = EmptyLabel(self.listFrame)
            self.widNew.grid(row=self.guiID, column=colID,sticky=Tk.W+Tk.E+Tk.S+Tk.N) 

    def updateMMB(self):
        try:
            self.data.mmbUpdate()
            # pyMMB.postInitialize()
            self.data.valid = True
            self.data.updateRepresentation()
            self.display()
        except pyMMB.MMBError as e:
            tkMessageBox.showerror("MMB Error", e.msg)


#############################################################################
class PolymerWidgetsLine(WidgetsLine):
    def __init__(self, listFrame, guiID, polymer, select=False, addCommand=None):
        WidgetsLine.__init__(self, listFrame, guiID, polymer)
        self.addCommand = addCommand


    ## Redefinition of base class changeAttribute because the behavior is different
    def changeAttribute(self, x, redisplay=True):
        """
        x.data is the name of the attribute to update
        x.get() return the new value of this attribute
        """
        self.data.removeFromSelection()
        attr_key = x.data
        attr = x.get()
        self.data.setAttribute(attr_key, attr)
        self.data.addToSelection()
        self.data.valid = False
        if redisplay:
            self.display()
        else:
            self.drawValidNewButtons(5)
            self.updateBgColor()
        # self.data.updateRepresentation()

    def display(self):
        # clear my widgets
        [w.destroy() for w in [self.__dict__[k] for k in self.__dict__ if k.startswith("wid")] if w]

        colNum = self.iter1()

        # self.widChainLabel = Tk.Label(self.listFrame, text=str(self.data.chainID))
        # self.widChainLabel.grid(row=self.guiID, column=colNum.next(),sticky=Tk.W+Tk.E+Tk.S+Tk.N)

        # Polymer's type
        polyType = MMB_UI.PolyTypes[self.data.polyType]
        self.widPolyTypeOption = ComboBox(self.listFrame, MMB_UI.PolyTypes, polyType, command=self.changeAttribute, data="polyType")
        self.widPolyTypeOption.grid(row=self.guiID, column=colNum.next(),sticky=Tk.W+Tk.E)
        if not self.data.new:
            self.widPolyTypeOption.config(state="disabled")

        self.widChainEntry = ChainIDEntry(self.listFrame, value= self.data.chainID, width=1, command=lambda x, redisplay=False: self.changeAttribute(x,redisplay), data="chainID")
        # self.widChainEntry.insert(0, self.data.chainID)
        self.widChainEntry.grid(row=self.guiID, column=colNum.next(),sticky=Tk.W+Tk.E+Tk.S+Tk.N)
        if not self.data.new:
            self.widChainEntry.config(state="disabled")
        else:
            self.widChainEntry.config(selectbackground="dodger blue")

        # First Res # textfield
        self.widFirstRes = IntegerEntry(self.listFrame, value=self.data.firstResNum, width=3, command=lambda x, redisplay=False: self.changeAttribute(x,redisplay), data="firstResNum")
        # self.widFirstRes.insert(0, self.data.firstResNum)
        self.widFirstRes.grid(row=self.guiID, column=colNum.next(),sticky=Tk.W+Tk.E+Tk.S+Tk.N)
        if not self.data.new:
            self.widFirstRes.config(state="disabled")
        else:
            self.widFirstRes.config(selectbackground="dodger blue")
        self.balloon.bind_widget(self.widFirstRes, msg="Number of the first res. in the sequence")

        # Sequence
        self.widSeqEntry = SequenceEntry(self.listFrame, value=self.data.sequence, width=50, command=lambda x, redisplay=False: self.changeAttribute(x,redisplay), data="sequence")
        # self.widSeqEntry.insert(0, self.data.sequence)
        self.widSeqEntry.grid(row=self.guiID, column=colNum.next(),sticky=Tk.W+Tk.E+Tk.S+Tk.N)
        self.balloon.bind_widget(self.widSeqEntry, msg="Polymer's sequence")
        if MMB_UI.polymersInitialized:
            self.widSeqEntry.config(state="disabled")
        if self.data.new:
            self.widSeqEntry.config(selectbackground="dodger blue")


        # Topolgy (new or from a pdb file)
        self.widTopoOption = ComboBox(self.listFrame, MMB_UI.Topologies, self.data.topo, command=self.changeAttribute, data="topo")
        self.widTopoOption.grid(row=self.guiID, column=colNum.next(),sticky=Tk.W+Tk.E)
        if self.data.new or MMB_UI.polymersInitialized:
            self.widTopoOption.config(state="disabled")

        self.drawValidNewButtons(colNum.next())
        
        self.varSelect.set(self.data.select)
        self.widSelect = Tk.Checkbutton(self.listFrame, text="", variable=self.varSelect, command=self.selectLine)
        self.widSelect.grid(row=self.guiID, column=colNum.next(),sticky=Tk.W+Tk.E+Tk.S+Tk.N)

        self.widShowHide = Tk.Button(self.listFrame, text="Show/Hide", command=self.data.show_hide)
        self.widShowHide.grid(row=self.guiID, column=colNum.next(),sticky=Tk.W+Tk.E+Tk.S+Tk.N)

        self.updateBgColor()
        self.applyTags()

#############################################################################
class BasePairWidgetsLine(WidgetsLine):
    def __init__(self, listFrame, guiID, pair, select=False, addCommand=None):
        WidgetsLine.__init__(self, listFrame, guiID, pair)
        self.addCommand = addCommand

    def display(self):
        # clear my widgets
        [w.destroy() for w in [self.__dict__[k] for k in self.__dict__ if k.startswith("wid")] if w]



        self.widLabel = Tk.Label(self.listFrame, text=str(self.guiID))
        self.widLabel.grid(row=self.guiID, column=0,sticky=Tk.W+Tk.E+Tk.S+Tk.N)

        p1 = MMB_UI.polymers[self.data.poly1]
        p2 = MMB_UI.polymers[self.data.poly2]

        self.widPoly1 = ComboBox(self.listFrame, MMB_UI.sortedNucleicAcidsIDs(), self.data.poly1, data="poly1", command=self.changeAttribute)
        self.widPoly1.grid(row=self.guiID, column=1,sticky=Tk.W+Tk.E)

        self.widRes1 = ComboBox(self.listFrame, p1.getResNamesIDsList(), p1.getResNameId(self.data.res1), data="res1", command=self.changeAttribute)
        self.widRes1.grid(row=self.guiID, column=2,sticky=Tk.W+Tk.E)

        self.widEdge1 = ComboBox(self.listFrame, MMB_UI.PairTypes, self.data.edge1, data="edge1", command=self.changeAttribute)
        self.widEdge1.grid(row=self.guiID, column=3,sticky=Tk.W+Tk.E)

        self.widPoly2 = ComboBox(self.listFrame, MMB_UI.sortedNucleicAcidsIDs(), self.data.poly2, data="poly2", command=self.changeAttribute)
        self.widPoly2.grid(row=self.guiID, column=4,sticky=Tk.W+Tk.E)

        self.widRes2 = ComboBox(self.listFrame, p2.getResNamesIDsList(), p2.getResNameId(self.data.res2), data="res2", command=self.changeAttribute)
        self.widRes2.grid(row=self.guiID, column=5,sticky=Tk.W+Tk.E)

        self.widEdge2 = ComboBox(self.listFrame, MMB_UI.PairTypes, self.data.edge2, data="edge2", command=self.changeAttribute)
        self.widEdge2.grid(row=self.guiID, column=6,sticky=Tk.W+Tk.E)

        self.widOrient = ComboBox(self.listFrame, MMB_UI.BondOrient, self.data.bondOrient, data="bondOrient", command=self.changeAttribute)
        self.widOrient.grid(row=self.guiID, column=7,sticky=Tk.W+Tk.E)

        self.drawValidNewButtons(8)

        self.varSelect.set(self.data.select)
        self.widSelect = Tk.Checkbutton(self.listFrame, text="", variable=self.varSelect, command=self.selectLine)
        self.widSelect.grid(row=self.guiID, column=9,sticky=Tk.W+Tk.E+Tk.S+Tk.N)

        # print "Select ", self.pair.MMBiD

        # Apply background color to all widgets
        self.updateBgColor()
        self.applyTags()

#############################################################################
class NucleicAcidDuplexWidgetsLine(WidgetsLine):
    def __init__(self, listFrame, guiID, duplex, select=False, addCommand=None):
        WidgetsLine.__init__(self, listFrame, guiID, duplex)
        self.addCommand = addCommand

    def display(self):
        # clear my widgets
        [w.destroy() for w in [self.__dict__[k] for k in self.__dict__ if k.startswith("wid")] if w]

        colNum = self.iter1()

        chainsList = MMB_UI.sortedChainsIDs()
        if self.data.new: chainsList.insert(0,"")
        self.widChain1 = ComboBox(self.listFrame, chainsList, self.data.chainID1, 
                                 data="chainID1", command=self.changeAttribute)
        self.widChain1.grid(row=self.guiID, column=colNum.next(),sticky=Tk.W+Tk.E)

        resNamesIdsList = []
        if self.widChain1.get() != "":
            poly = MMB_UI.polymers[self.data.chainID1]
            resNamesIdsList = poly.getResNamesIDsList()
            self.widStartRes1 = ComboBox(self.listFrame, resNamesIdsList,  self.data.getResStart1NameId(), 
                                     data="residueStart1", command=self.changeAttribute)

            resNamesIdsList = poly.getResNamesIDsList()
            self.widEndRes1 = ComboBox(self.listFrame, resNamesIdsList,  self.data.getResEnd1NameId(), 
                                     data="residueEnd1", command=self.changeAttribute)
        else:
            self.widStartRes1 = EmptyLabel(self.listFrame)
            self.widEndRes1   = EmptyLabel(self.listFrame)

        self.widStartRes1.grid(row=self.guiID, column=colNum.next(),sticky=Tk.W+Tk.E+Tk.S+Tk.N)
        self.widEndRes1.grid(row=self.guiID, column=colNum.next(),sticky=Tk.W+Tk.E+Tk.S+Tk.N)

        self.widChain2 = ComboBox(self.listFrame, chainsList, self.data.chainID2, 
                                 data="chainID2", command=self.changeAttribute)
        self.widChain2.grid(row=self.guiID, column=colNum.next(),sticky=Tk.W+Tk.E)

        resNamesIdsList = []
        if self.widChain2.get() != "":
            poly = MMB_UI.polymers[self.data.chainID2]
            resNamesIdsList = poly.getResNamesIDsList()
            self.widStartRes2 = ComboBox(self.listFrame, resNamesIdsList,  self.data.getResStart2NameId(), 
                                     data="residueStart2", command=self.changeAttribute)

            resNamesIdsList = poly.getResNamesIDsList()
            self.widEndRes2 = ComboBox(self.listFrame, resNamesIdsList,  self.data.getResEnd2NameId(), 
                                     data="residueEnd2", command=self.changeAttribute)
        else:
            self.widStartRes2 = EmptyLabel(self.listFrame)
            self.widEndRes2   = EmptyLabel(self.listFrame)

        self.widStartRes2.grid(row=self.guiID, column=colNum.next(),sticky=Tk.W+Tk.E+Tk.S+Tk.N)
        self.widEndRes2.grid(row=self.guiID, column=colNum.next(),sticky=Tk.W+Tk.E+Tk.S+Tk.N)

        self.drawValidNewButtons(colNum.next())


        self.widSelect = Tk.Checkbutton(self.listFrame, text="", variable=self.varSelect, command=self.selectLine)
        self.widSelect.grid(row=self.guiID, column=colNum.next(),sticky=Tk.W+Tk.E+Tk.S+Tk.N)

        self.updateBgColor()
        self.applyTags()


#############################################################################
class MobilizerWidgetsLine(WidgetsLine):
    def __init__(self, listFrame, guiID, mobilizer, select=False, addCommand=None):
        WidgetsLine.__init__(self, listFrame, guiID, mobilizer)
        self.addCommand = addCommand

    def display(self):
        # clear my widgets
        [w.destroy() for w in [self.__dict__[k] for k in self.__dict__ if k.startswith("wid")] if w]

        self.widLabel = Tk.Label(self.listFrame, text=str(self.guiID))
        self.widLabel.grid(row=self.guiID, column=0,sticky=Tk.W+Tk.E+Tk.S+Tk.N)

        self.widMobiLabel = Tk.Label(self.listFrame, text="Mobilizer")
        self.widMobiLabel.grid(row=self.guiID, column=1,sticky=Tk.W+Tk.E+Tk.S+Tk.N)

        self.widMobiType = ComboBox(self.listFrame, MMB_UI.MobilityTypes, self.data.mobility, 
                                 data="mobility", command=self.changeAttribute)
        self.widMobiType.grid(row=self.guiID, column=2,sticky=Tk.W+Tk.E)

        chainsList = MMB_UI.sortedChainsIDs()
        if self.data.new: chainsList.insert(0,"All")
        self.widChain = ComboBox(self.listFrame, chainsList, self.data.chainID, 
                                 data="chainID", command=self.changeAttribute)
        self.widChain.grid(row=self.guiID, column=3,sticky=Tk.W+Tk.E)

        if self.widChain.get() != "All":
            poly = MMB_UI.polymers[self.data.chainID]
            resNamesIdsList = poly.getResNamesIDsList()
            if self.data.new: resNamesIdsList.insert(0,"All")
            self.widStartRes = ComboBox(self.listFrame, resNamesIdsList,  self.data.getResStartNameId(), 
                                     data="resStart", command=self.changeAttribute)

            if self.widStartRes.get() != "All":
                resNamesIdsList = poly.getResNamesIDsList(self.data.resStart)
                self.widEndRes = ComboBox(self.listFrame, resNamesIdsList,  self.data.getResEndNameId(), 
                                         data="resEnd", command=self.changeAttribute)
            else:
                self.widEndRes = EmptyLabel(self.listFrame)
        else:
            self.widStartRes = EmptyLabel(self.listFrame)
            self.widEndRes   = EmptyLabel(self.listFrame)

        self.widStartRes.grid(row=self.guiID, column=4,sticky=Tk.W+Tk.E+Tk.S+Tk.N)
        self.widEndRes.grid(row=self.guiID, column=5,sticky=Tk.W+Tk.E+Tk.S+Tk.N)

        self.drawValidNewButtons(6)


        self.widSelect = Tk.Checkbutton(self.listFrame, text="", variable=self.varSelect, command=self.selectLine)
        self.widSelect.grid(row=self.guiID, column=7,sticky=Tk.W+Tk.E+Tk.S+Tk.N)

        self.updateBgColor()
        self.applyTags()

#############################################################################
class MobilizerWithinWidgetsLine(WidgetsLine):
    def __init__(self, listFrame, guiID, mobilizerWithin, select=False, addCommand=None):
        WidgetsLine.__init__(self, listFrame, guiID, mobilizerWithin)
        self.addCommand = addCommand
        self.widRadius = None

    def changeRadius(self, *args):
        self.data.removeFromSelection()
        self.data.radius = float(self.widRadius.get())/10.0
        self.data.valid = False
        self.select(True)
        self.display()

    def display(self):
        # clear my widgets
        [w.destroy() for w in [self.__dict__[k] for k in self.__dict__ if k.startswith("wid") and k!="widRadius"] if w]

        self.widLabel = EmptyLabel(self.listFrame)
        self.widLabel.grid(row=self.guiID, column=0,sticky=Tk.W+Tk.E+Tk.S+Tk.N)

        self.widMobiLabel = Tk.Label(self.listFrame, text="MobilizersWithin")
        self.widMobiLabel.grid(row=self.guiID, column=1,sticky=Tk.W+Tk.E+Tk.S+Tk.N)

        self.widMobiType = ComboBox(self.listFrame, MMB_UI.MobilityTypes, self.data.mobility, 
                                 data="mobility", command=self.changeAttribute)
        self.widMobiType.grid(row=self.guiID, column=2,sticky=Tk.W+Tk.E)

        chainsList = MMB_UI.sortedChainsIDs()
        self.widChain = ComboBox(self.listFrame, chainsList, self.data.chainID, 
                                 data="chainID", command=self.changeAttribute)
        self.widChain.grid(row=self.guiID, column=3,sticky=Tk.W+Tk.E)

        resID = ""
        resNamesIdsList = [""]
        if self.data.chainID != "":
            poly = MMB_UI.polymers[self.data.chainID]
            resNamesIdsList = poly.getResNamesIDsList()
            if self.data.resID != "":
                resID = self.data.getResNameId()
        
        self.widRes = ComboBox(self.listFrame, resNamesIdsList, resID, data="resID", command=self.changeAttribute)
        self.widRes.grid(row=self.guiID, column=4,sticky=Tk.W+Tk.E+Tk.S+Tk.N)

        if not self.widRadius:
            self.widRadius = NumericEntry(self.listFrame,
                                          value=self.data.radius*10.0,
                                          command=self.changeAttribute,
                                          data="radius", 
                                          width=5)
            if self.data.new: self.widRadius.config(selectbackground="dodger blue")
            self.widRadius.grid(row=self.guiID, column=5,sticky=Tk.W+Tk.E+Tk.S+Tk.N)

        self.drawValidNewButtons(6)

        self.widSelect = Tk.Checkbutton(self.listFrame, text="", variable=self.varSelect, command=self.selectLine)
        self.widSelect.grid(row=self.guiID, column=7,sticky=Tk.W+Tk.E+Tk.S+Tk.N)

        self.updateBgColor()
        self.applyTags()

#############################################################################
class RootMobilizerWidgetsLine(WidgetsLine):
    def __init__(self, listFrame, guiID, rootmobilizer, select=False, addCommand=None):
        WidgetsLine.__init__(self, listFrame, guiID, rootmobilizer)
        self.addCommand = addCommand

    def display(self):
        # clear my widgets
        [w.destroy() for w in [self.__dict__[k] for k in self.__dict__ if k.startswith("wid")] if w]

        colNum = self.iter1()

        self.widLabel = Tk.Label(self.listFrame, text="RootMobilizer")
        self.widLabel.grid(row=self.guiID, column=colNum.next(),sticky=Tk.W+Tk.E+Tk.S+Tk.N)

        self.widChainLabel = Tk.Label(self.listFrame, text=self.data.chainID)
        self.widChainLabel.grid(row=self.guiID, column=colNum.next(),sticky=Tk.W+Tk.E+Tk.S+Tk.N)

        mobiTypesList = MMB_UI.RootMobilizer.rootMobilizerTypes
        self.widMobi = ComboBox(self.listFrame, mobiTypesList, self.data.rootMobilizer, 
                                 data="rootMobilizer", command=self.changeAttribute)
        self.widMobi.grid(row=self.guiID, column=colNum.next(),sticky=Tk.W+Tk.E)

        self.drawValidNewButtons(colNum.next())

        self.widSelect = Tk.Checkbutton(self.listFrame, text="", variable=self.varSelect, command=self.selectLine)
        self.widSelect.grid(row=self.guiID, column=colNum.next(),sticky=Tk.W+Tk.E+Tk.S+Tk.N)

        self.updateBgColor()
        self.applyTags()

#############################################################################
class ConstraintWidgetsLine(WidgetsLine):
    def __init__(self, listFrame, guiID, constraint, select=False, addCommand=None):
        WidgetsLine.__init__(self, listFrame, guiID, constraint)
        self.addCommand = addCommand

    def display(self):
        # clear my widgets
        [w.destroy() for w in [self.__dict__[k] for k in self.__dict__ if k.startswith("wid")] if w]

        self.widLabel = Tk.Label(self.listFrame, text=str(self.guiID))
        self.widLabel.grid(row=self.guiID, column=0,sticky=Tk.W+Tk.E+Tk.S+Tk.N)

        chainsList = [""] + MMB_UI.sortedChainsIDs()
        
        # Chain 1
        self.widChain1 = ComboBox(self.listFrame, chainsList, self.data.chain1, 
                                 data="chain1", command=self.changeAttribute)
        self.widChain1.grid(row=self.guiID, column=1,sticky=Tk.W+Tk.E)

        # Residue 1
        res1List = [""]
        poly1 = MMB_UI.BioPolymer()
        res1State = "normal"
        if self.data.chain1 != "": 
            poly1 = MMB_UI.polymers[self.data.chain1]
            res1List = res1List + poly1.getResNamesIDsList()
        else:
            res1State = "disabled"
        self.widRes1 = ComboBox(self.listFrame, res1List, self.data.getRes1NameId(), 
                                 data="res1", command=self.changeAttribute)
        self.widRes1.config(state=res1State)
        self.widRes1.grid(row=self.guiID, column=2,sticky=Tk.W+Tk.E)

        # Atom 1
        atom1List = [""]
        atom1State = "normal"
        if self.data.chain1 != "" and self.data.res1 != "":
            atom1List = atom1List + self.data.getRes1AtomsList()
        elif self.data.res1 == "":
            atom1State = "disabled"
        self.widAtom1 = ComboBox(self.listFrame, atom1List, self.data.atom1, 
                                 data="atom1", command=self.changeAttribute)
        self.widAtom1.config(state=atom1State)
        self.widAtom1.grid(row=self.guiID, column=3,sticky=Tk.W+Tk.E)

        # Weld label
        self.widWeldLabel = Tk.Label(self.listFrame, text="Weld")
        self.widWeldLabel.grid(row=self.guiID, column=4,sticky=Tk.W+Tk.E+Tk.S+Tk.N)

        # Chain 2
        chainsList = ["Ground"] + chainsList[1:] # we skip the " " chain
        self.widChain2 = ComboBox(self.listFrame, chainsList, self.data.chain2, 
                                 data="chain2", command=self.changeAttribute)
        self.widChain2.grid(row=self.guiID, column=5,sticky=Tk.W+Tk.E)

        # Residue 2
        res2List = [""]
        poly2 = MMB_UI.BioPolymer()
        res2State = "normal"
        if self.data.chain2 != "Ground": 
            poly2 = MMB_UI.polymers[self.data.chain2]
            res2List = res2List + poly2.getResNamesIDsList()
        else:
            res2State = "disabled"
        self.widRes2 = ComboBox(self.listFrame, res2List, self.data.getRes2NameId(), 
                                 data="res2", command=self.changeAttribute)
        self.widRes2.config(state=res2State)
        self.widRes2.grid(row=self.guiID, column=6,sticky=Tk.W+Tk.E)

        # Atom 2
        atom2List = [""]
        atom2State = "normal"
        if self.data.chain2 != "Ground" and self.data.res2 != "":
            atom2List = atom2List + self.data.getRes2AtomsList()
        elif self.data.res2 == "":
            atom2State = "disabled"
        self.widAtom2 = ComboBox(self.listFrame, atom2List, self.data.atom2, 
                                 data="atom2", command=self.changeAttribute)
        self.widAtom2.config(state=atom2State)
        self.widAtom2.grid(row=self.guiID, column=7,sticky=Tk.W+Tk.E)

        self.drawValidNewButtons(8)

        self.widSelect = Tk.Checkbutton(self.listFrame, text="", variable=self.varSelect, command=self.selectLine)
        self.widSelect.grid(row=self.guiID, column=9,sticky=Tk.W+Tk.E+Tk.S+Tk.N)

        self.updateBgColor()
        self.applyTags()

#############################################################################
class AllResiduesWithinWidgetsLine(WidgetsLine):
    def __init__(self, listFrame, guiID, allResiduesWithin, select=False, addCommand=None):
        WidgetsLine.__init__(self, listFrame, guiID, allResiduesWithin)
        self.addCommand = addCommand
        self.widRadius = None

    def changeRadius(self, *args):
        self.data.removeFromSelection()
        self.data.radius = float(self.widRadius.get())/10.0
        self.data.valid = False
        self.select(True)
        self.display()

    def display(self):
        # clear my widgets
        [w.destroy() for w in [self.__dict__[k] for k in self.__dict__ if k.startswith("wid") and k!="widRadius"] if w]

        self.widLabel = Tk.Label(self.listFrame, text=str(self.guiID))
        self.widLabel.grid(row=self.guiID, column=0,sticky=Tk.W+Tk.E+Tk.S+Tk.N)

        self.widARWLabel = Tk.Label(self.listFrame, text="AllResiduesWithin")
        self.widARWLabel.grid(row=self.guiID, column=1,sticky=Tk.W+Tk.E+Tk.S+Tk.N)

        if not self.widRadius:
            self.widRadius = NumericEntry(self.listFrame,
                                          value=self.data.radius*10.0,
                                          command=self.changeRadius, 
                                          width=5)
            if self.data.new: self.widRadius.config(selectbackground="dodger blue")
            self.widRadius.grid(row=self.guiID, column=2,sticky=Tk.W+Tk.E+Tk.S+Tk.N)

        chainsList = MMB_UI.sortedChainsIDs()
        self.widChain = ComboBox(self.listFrame, chainsList, self.data.chain, 
                                 data="chain", command=self.changeAttribute)
        self.widChain.grid(row=self.guiID, column=3,sticky=Tk.W+Tk.E)

        resID = ""
        resNamesIdsList = [""]
        if self.data.chain != "":
            poly = MMB_UI.polymers[self.data.chain]
            resNamesIdsList = poly.getResNamesIDsList()
            if self.data.residue != "":
                resID = self.data.getResNameId()
        
        self.widRes = ComboBox(self.listFrame, resNamesIdsList, resID, data="residue", command=self.changeAttribute)
        self.widRes.grid(row=self.guiID, column=4,sticky=Tk.W+Tk.E+Tk.S+Tk.N)

        self.drawValidNewButtons(5)

        self.widSelect = Tk.Checkbutton(self.listFrame, text="", variable=self.varSelect, command=self.selectLine)
        self.widSelect.grid(row=self.guiID, column=6,sticky=Tk.W+Tk.E+Tk.S+Tk.N)

        self.updateBgColor()
        self.applyTags()

#############################################################################
class IncludeNonBondAllAtomsInResidueWidgetsLine(WidgetsLine):
    def __init__(self, listFrame, guiID, includeNonBondAllAtomsInResidue, select=False, addCommand=None):
        WidgetsLine.__init__(self, listFrame, guiID, includeNonBondAllAtomsInResidue)
        self.addCommand = addCommand
        self.widRadius = None

    def display(self):
        # clear my widgets
        [w.destroy() for w in [self.__dict__[k] for k in self.__dict__ if k.startswith("wid") and k!="widRadius"] if w]

        self.widLabel = Tk.Label(self.listFrame, text=str(self.guiID))
        self.widLabel.grid(row=self.guiID, column=0,sticky=Tk.W+Tk.E+Tk.S+Tk.N)

        self.widARWLabel = Tk.Label(self.listFrame, text="IncludeNonBondAllAtomsInResidue")
        self.widARWLabel.grid(row=self.guiID, column=1,sticky=Tk.W+Tk.E+Tk.S+Tk.N)

        chainsList = MMB_UI.sortedChainsIDs()
        self.widChain = ComboBox(self.listFrame, chainsList, self.data.chain, 
                                 data="chain", command=self.changeAttribute)
        self.widChain.grid(row=self.guiID, column=2,sticky=Tk.W+Tk.E)

        resID = ""
        resNamesIdsList = [""]
        if self.data.chain != "":
            poly = MMB_UI.polymers[self.data.chain]
            resNamesIdsList = poly.getResNamesIDsList()
            if self.data.residue != "":
                resID = self.data.getResStartNameId()
        
        self.widRes = ComboBox(self.listFrame, resNamesIdsList, resID, data="residue", command=self.changeAttribute)
        self.widRes.grid(row=self.guiID, column=3,sticky=Tk.W+Tk.E+Tk.S+Tk.N)
        self.balloon.bind_widget(self.widRes, msg="From")

        if self.data.new:
            resEnd = ""
            if self.data.chain != "" and self.data.resEnd != "":
                resNamesIdsList = poly.getResNamesIDsList(self.data.residue)
                resEnd = self.data.getResEndNameId()
            self.widResEnd = ComboBox(self.listFrame, resNamesIdsList, resEnd, data="resEnd", command=self.changeAttribute)
            self.widResEnd.grid(row=self.guiID, column=4,sticky=Tk.W+Tk.E+Tk.S+Tk.N)
            self.balloon.bind_widget(self.widResEnd, msg="To")
        else:
            self.widResEnd = EmptyLabel(self.listFrame)
            self.widResEnd.grid(row=self.guiID, column=4,sticky=Tk.W+Tk.E+Tk.S+Tk.N)

        self.drawValidNewButtons(5)

        self.widSelect = Tk.Checkbutton(self.listFrame, text="", variable=self.varSelect, command=self.selectLine)
        self.widSelect.grid(row=self.guiID, column=6,sticky=Tk.W+Tk.E+Tk.S+Tk.N)

        self.updateBgColor()
        self.applyTags()

#############################################################################
class DeactivatePhysicsWidgetsLine(WidgetsLine):
    def __init__(self, listFrame, guiID, data, select=False, addCommand=None):
        WidgetsLine.__init__(self, listFrame, guiID, data)
        self.addCommand = addCommand

        self.varActivate = Tk.BooleanVar()
        self.varActivate.set(self.data.activePhysics)

    def changeActivePhysics(self, *args):
        self.data.removeFromSelection()
        self.data.activePhysics = self.varActivate.get()
        self.data.valid = False
        self.select(True)
        self.display()

    def display(self):
        # clear my widgets
        [w.destroy() for w in [self.__dict__[k] for k in self.__dict__ if k.startswith("wid")] if w]

        colNum = self.iter1()

        self.widLabel = Tk.Label(self.listFrame, text="Physics")
        self.widLabel.grid(row=self.guiID, column=colNum.next(),sticky=Tk.W+Tk.E+Tk.S+Tk.N)

        self.widChainLabel = Tk.Label(self.listFrame, text=self.data.chainID)
        self.widChainLabel.grid(row=self.guiID, column=colNum.next(),sticky=Tk.W+Tk.E+Tk.S+Tk.N)

        self.widActivate = Tk.Checkbutton(self.listFrame, text="", variable=self.varActivate, command=self.changeActivePhysics)
        self.widActivate.grid(row=self.guiID, column=colNum.next(), sticky=Tk.W+Tk.E+Tk.S+Tk.N)

        self.drawValidNewButtons(colNum.next())

        self.widSelect = Tk.Checkbutton(self.listFrame, text="", variable=self.varSelect, command=self.selectLine)
        self.widSelect.grid(row=self.guiID, column=colNum.next(),sticky=Tk.W+Tk.E+Tk.S+Tk.N)

        self.updateBgColor()
        self.applyTags()


#############################################################################
class ContactWidgetsLine(WidgetsLine):
    def __init__(self, listFrame, guiID, contact, select=False, addCommand=None):
        WidgetsLine.__init__(self, listFrame, guiID, contact)
        self.addCommand = addCommand

    def display(self):
        # clear my widgets
        [w.destroy() for w in [self.__dict__[k] for k in self.__dict__ if k.startswith("wid")] if w]

        self.widLabel = Tk.Label(self.listFrame, text=str(self.guiID))
        self.widLabel.grid(row=self.guiID, column=0,sticky=Tk.W+Tk.E+Tk.S+Tk.N)

        self.widContactLabel = Tk.Label(self.listFrame, text="Contact")
        self.widContactLabel.grid(row=self.guiID, column=1,sticky=Tk.W+Tk.E+Tk.S+Tk.N)

        self.widContactType = ComboBox(self.listFrame, MMB_UI.Contact.typesList, self.data.contactScheme, 
                                 data="contactScheme", command=self.changeAttribute)
        self.widContactType.grid(row=self.guiID, column=2,sticky=Tk.W+Tk.E)

        chainsList = MMB_UI.sortedChainsIDs()
        if self.data.new: chainsList.insert(0,"")
        self.widChain = ComboBox(self.listFrame, chainsList, self.data.chainID, 
                                 data="chainID", command=self.changeAttribute)
        self.widChain.grid(row=self.guiID, column=3,sticky=Tk.W+Tk.E)

        resNamesIdsList = []
        if self.widChain.get() != "":
            poly = MMB_UI.polymers[self.data.chainID]
            resNamesIdsList = poly.getResNamesIDsList()
            resNamesIdsList.insert(0,"All")
            self.widStartRes = ComboBox(self.listFrame, resNamesIdsList,  self.data.getResStartNameId(), 
                                     data="resStart", command=self.changeAttribute)

            if self.widStartRes.get() != "All":
                resNamesIdsList = poly.getResNamesIDsList(self.data.resStart)
                self.widEndRes = ComboBox(self.listFrame, resNamesIdsList,  self.data.getResEndNameId(), 
                                         data="resEnd", command=self.changeAttribute)
            else:
                self.widEndRes = EmptyLabel(self.listFrame)
        else:
            self.widStartRes = EmptyLabel(self.listFrame)
            self.widEndRes   = EmptyLabel(self.listFrame)

        self.widStartRes.grid(row=self.guiID, column=4,sticky=Tk.W+Tk.E+Tk.S+Tk.N)
        self.widEndRes.grid(row=self.guiID, column=5,sticky=Tk.W+Tk.E+Tk.S+Tk.N)

        self.drawValidNewButtons(6)


        self.widSelect = Tk.Checkbutton(self.listFrame, text="", variable=self.varSelect, command=self.selectLine)
        self.widSelect.grid(row=self.guiID, column=7,sticky=Tk.W+Tk.E+Tk.S+Tk.N)

        self.updateBgColor()
        self.applyTags()

#############################################################################
class ContactWithinWidgetsLine(WidgetsLine):
    def __init__(self, listFrame, guiID, contactWithin, select=False, addCommand=None):
        WidgetsLine.__init__(self, listFrame, guiID, contactWithin)
        self.addCommand = addCommand
        self.widRadius = None

    def changeRadius(self, *args):
        self.data.removeFromSelection()
        self.data.radius = float(self.widRadius.get())/10.0
        self.data.valid = False
        self.select(True)
        self.display()

    def display(self):
        # clear my widgets
        [w.destroy() for w in [self.__dict__[k] for k in self.__dict__ if k.startswith("wid") and k!="widRadius"] if w]

        self.widLabel = EmptyLabel(self.listFrame)
        self.widLabel.grid(row=self.guiID, column=0,sticky=Tk.W+Tk.E+Tk.S+Tk.N)

        self.widContactLabel = Tk.Label(self.listFrame, text="ContactsWithin")
        self.widContactLabel.grid(row=self.guiID, column=1,sticky=Tk.W+Tk.E+Tk.S+Tk.N)

        self.widContactType = ComboBox(self.listFrame, MMB_UI.Contact.typesList, self.data.contactScheme, 
                                 data="contactScheme", command=self.changeAttribute)
        self.widContactType.grid(row=self.guiID, column=2,sticky=Tk.W+Tk.E)

        chainsList = MMB_UI.sortedChainsIDs()
        self.widChain = ComboBox(self.listFrame, chainsList, self.data.chainID, 
                                 data="chainID", command=self.changeAttribute)
        self.widChain.grid(row=self.guiID, column=3,sticky=Tk.W+Tk.E)

        resID = ""
        resNamesIdsList = [""]
        if self.data.chainID != "":
            poly = MMB_UI.polymers[self.data.chainID]
            resNamesIdsList = poly.getResNamesIDsList()
            if self.data.resID != "":
                resID = self.data.getResNameId()
        
        self.widRes = ComboBox(self.listFrame, resNamesIdsList, resID, data="resID", command=self.changeAttribute)
        self.widRes.grid(row=self.guiID, column=4,sticky=Tk.W+Tk.E+Tk.S+Tk.N)

        if not self.widRadius:
            self.widRadius = NumericEntry(self.listFrame,
                                          value=self.data.radius*10.0,
                                          command=self.changeRadius, 
                                          width=5)
            if self.data.new: self.widRadius.config(selectbackground="dodger blue")
            self.widRadius.grid(row=self.guiID, column=5,sticky=Tk.W+Tk.E+Tk.S+Tk.N)

        self.drawValidNewButtons(6)

        self.widSelect = Tk.Checkbutton(self.listFrame, text="", variable=self.varSelect, command=self.selectLine)
        self.widSelect.grid(row=self.guiID, column=7,sticky=Tk.W+Tk.E+Tk.S+Tk.N)

        self.updateBgColor()
        self.applyTags()

#############################################################################
class ThreadingWidgetsLine(WidgetsLine):
    def __init__(self, listFrame, guiID, threading, select=False, addCommand=None):
        WidgetsLine.__init__(self, listFrame, guiID, threading)
        self.addCommand = addCommand
        self.widForceConstant = None
        self.varBackbone = Tk.BooleanVar()
        self.varBackbone.set(self.data.backboneOnly)

    def changeBackbone(self, *args):
        self.data.removeFromSelection()
        self.data.backboneOnly = self.varBackbone.get()
        self.data.valid = False
        self.select(True)
        self.display()

    def changeConstant(self, *args):
        self.data.removeFromSelection()
        self.data.forceConstant = float(self.widForceConstant.get())
        self.data.valid = False
        self.select(True)
        self.display()

    def display(self):
        # clear my widgets
        [w.destroy() for w in [self.__dict__[k] for k in self.__dict__ if k.startswith("wid") and k !="widForceConstant"] if w]

        colNum = self.iter1()

        self.widLabel = Tk.Label(self.listFrame, text=str(self.guiID))
        self.widLabel.grid(row=self.guiID, column=colNum.next(),sticky=Tk.W+Tk.E+Tk.S+Tk.N)

        chainsList = MMB_UI.sortedChainsIDs()
        if self.data.new: chainsList.insert(0,"")
        self.widChain1 = ComboBox(self.listFrame, chainsList, self.data.chainID1, 
                                 data="chainID1", command=self.changeAttribute)
        self.widChain1.grid(row=self.guiID, column=colNum.next(),sticky=Tk.W+Tk.E)

        resNamesIdsList = [""]
        if self.widChain1.get() != "":
            poly = MMB_UI.polymers[self.data.chainID1]
            resNamesIdsList = poly.getResNamesIDsList()
        self.widStartRes1 = ComboBox(self.listFrame, resNamesIdsList, self.data.getResStart1NameId(), 
                                     data="residueStart1", command=self.changeAttribute)

        if self.data.getResStart1NameId():
            resNamesIdsList = poly.getResNamesIDsList(self.data.residueStart1)
        self.widEndRes1 = ComboBox(self.listFrame, resNamesIdsList, self.data.getResEnd1NameId(), 
                                   data="residueEnd1", command=self.changeAttribute)

        self.widStartRes1.grid(row=self.guiID, column=colNum.next(),sticky=Tk.W+Tk.E+Tk.S+Tk.N)
        self.widEndRes1.grid(row=self.guiID, column=colNum.next(),sticky=Tk.W+Tk.E+Tk.S+Tk.N)

        self.widChain2 = ComboBox(self.listFrame, chainsList, self.data.chainID2, 
                                 data="chainID2", command=self.changeAttribute)
        self.widChain2.grid(row=self.guiID, column=colNum.next(),sticky=Tk.W+Tk.E)

        resNamesIdsList = [""]
        if self.widChain2.get() != "":
            poly = MMB_UI.polymers[self.data.chainID2]
            resNamesIdsList = poly.getResNamesIDsList()
        self.widStartRes2 = ComboBox(self.listFrame, resNamesIdsList, self.data.getResStart2NameId(), 
                                     data="residueStart2", command=self.changeAttribute)

        if self.data.getResStart2NameId():
            resNamesIdsList = poly.getResNamesIDsList(self.data.residueStart2)
        self.widEndRes2 = ComboBox(self.listFrame, resNamesIdsList, self.data.getResEnd2NameId(), 
                                   data="residueEnd2", command=self.changeAttribute)

        self.widStartRes2.grid(row=self.guiID, column=colNum.next(),sticky=Tk.W+Tk.E+Tk.S+Tk.N)
        self.widEndRes2.grid(row=self.guiID, column=colNum.next(),sticky=Tk.W+Tk.E+Tk.S+Tk.N)

        if not self.widForceConstant:
            self.widForceConstant = NumericEntry(self.listFrame,
                                          value=self.data.forceConstant,
                                          command=self.changeConstant, 
                                          width=5)
            if self.data.new: self.widForceConstant.config(selectbackground="dodger blue")
            self.widForceConstant.grid(row=self.guiID, column=colNum.next(),sticky=Tk.W+Tk.E+Tk.S+Tk.N)
        else:
            colNum.next()

        self.widBackbone = Tk.Checkbutton(self.listFrame, text="", variable=self.varBackbone, command=self.changeBackbone)
        self.widBackbone.grid(row=self.guiID, column=colNum.next(), sticky=Tk.W+Tk.E+Tk.S+Tk.N)

        if not self.data.new:
            self.widMatch = Tk.Button(self.listFrame, text="Match", command=self.data.matchSelections)
            self.balloon.bind_widget(self.widMatch, msg="RMS fit the residues with Chimera's match command")
        else:
            self.widMatch = EmptyLabel(self.listFrame)
        self.widMatch.grid(row=self.guiID, column=colNum.next(),sticky=Tk.W+Tk.E+Tk.S+Tk.N) 

        self.drawValidNewButtons(colNum.next())


        self.widSelect = Tk.Checkbutton(self.listFrame, text="", variable=self.varSelect, command=self.selectLine)
        self.widSelect.grid(row=self.guiID, column=colNum.next(),sticky=Tk.W+Tk.E+Tk.S+Tk.N)

        self.updateBgColor()
        self.applyTags()


#############################################################################
class GappedThreadingWidgetsLine(WidgetsLine):
    def __init__(self, listFrame, guiID, gappedThreading, select=False, addCommand=None):
        WidgetsLine.__init__(self, listFrame, guiID, gappedThreading)
        self.addCommand = addCommand
        self.widForceConstant = None
        self.varBackbone = Tk.BooleanVar()
        self.varBackbone.set(self.data.backboneOnly)

    def changeBackbone(self, *args):
        self.data.removeFromSelection()
        self.data.backboneOnly = self.varBackbone.get()
        self.data.valid = False
        self.select(True)
        self.display()

    def changeConstant(self, *args):
        self.data.removeFromSelection()
        self.data.forceConstant = float(self.widForceConstant.get())
        self.data.valid = False
        self.select(True)
        self.display()

    def display(self):
        # clear my widgets
        [w.destroy() for w in [self.__dict__[k] for k in self.__dict__ if k.startswith("wid") and k !="widForceConstant"] if w]

        colNum = self.iter1()

        self.widLabel = Tk.Label(self.listFrame, text=str(self.guiID))
        self.widLabel.grid(row=self.guiID, column=colNum.next(),sticky=Tk.W+Tk.E+Tk.S+Tk.N)

        chainsList = MMB_UI.sortedChainsIDs()
        if self.data.new: chainsList.insert(0,"")
        self.widChain1 = ComboBox(self.listFrame, chainsList, self.data.chainID1, 
                                 data="chainID1", command=self.changeAttribute)
        self.widChain1.grid(row=self.guiID, column=colNum.next(),sticky=Tk.W+Tk.E)

        self.widChain2 = ComboBox(self.listFrame, chainsList, self.data.chainID2, 
                                 data="chainID2", command=self.changeAttribute)
        self.widChain2.grid(row=self.guiID, column=colNum.next(),sticky=Tk.W+Tk.E)

        if not self.widForceConstant:
            self.widForceConstant = NumericEntry(self.listFrame,
                                          value=self.data.forceConstant,
                                          command=self.changeConstant, 
                                          width=5)
            if self.data.new: self.widForceConstant.config(selectbackground="dodger blue")
            self.widForceConstant.grid(row=self.guiID, column=colNum.next(),sticky=Tk.W+Tk.E+Tk.S+Tk.N)
        else:
            colNum.next()

        self.widBackbone = Tk.Checkbutton(self.listFrame, text="", variable=self.varBackbone, command=self.changeBackbone)
        self.widBackbone.grid(row=self.guiID, column=colNum.next(), sticky=Tk.W+Tk.E+Tk.S+Tk.N)

        if not self.data.new:
            self.widMatch = Tk.Button(self.listFrame, text="Match", command=self.data.matchSelections)
            self.balloon.bind_widget(self.widMatch, msg="RMS fit the residues with Chimera's match command")
        else:
            self.widMatch = EmptyLabel(self.listFrame)
        self.widMatch.grid(row=self.guiID, column=colNum.next(),sticky=Tk.W+Tk.E+Tk.S+Tk.N) 

        self.drawValidNewButtons(colNum.next())


        self.widSelect = Tk.Checkbutton(self.listFrame, text="", variable=self.varSelect, command=self.selectLine)
        self.widSelect.grid(row=self.guiID, column=colNum.next(),sticky=Tk.W+Tk.E+Tk.S+Tk.N)

        self.updateBgColor()
        self.applyTags()


#############################################################################
class AtomSpringWidgetsLine(WidgetsLine):
    def __init__(self, listFrame, guiID, atomspring, select=False, addCommand=None):
        WidgetsLine.__init__(self, listFrame, guiID, atomspring)
        self.addCommand = addCommand
        self.widForceConstant = None
        self.widDeadLength = None
        self.varTether = Tk.BooleanVar()
        self.varTether.set(self.data.tether)

    def changeTether(self, *args):
        # self.data.removeFromSelection()
        self.data.tether = self.varTether.get()
        self.data.valid = False
        self.select(True)
        self.display()

    def changeConstant(self, *args):
        self.data.forceConstant = float(self.widForceConstant.get())
        self.data.valid = False
        self.select(True)
        self.display()

    def changeDeadLength(self, *args):
        self.data.deadLength = float(self.widDeadLength.get())/10.0
        self.data.valid = False
        self.select(True)
        self.display()

    def display(self):
        # clear my widgets
        [w.destroy() for w in [self.__dict__[k] for k in self.__dict__ if k.startswith("wid") and k not in ["widForceConstant", "widDeadLength"]] if w]

        colNum = self.iter1()

        self.widLabel = Tk.Label(self.listFrame, text=str(self.guiID))
        self.widLabel.grid(row=self.guiID, column=colNum.next(),sticky=Tk.W+Tk.E+Tk.S+Tk.N)

        chainsList = [""] + MMB_UI.sortedChainsIDs()
        
        # Chain 1
        self.widChain1 = ComboBox(self.listFrame, chainsList, self.data.chain1, 
                                 data="chain1", command=self.changeAttribute)
        self.widChain1.grid(row=self.guiID, column=colNum.next(),sticky=Tk.W+Tk.E)

        # Residue 1
        res1List = [""]
        poly1 = MMB_UI.BioPolymer()
        res1State = "normal"
        if self.data.chain1 != "": 
            poly1 = MMB_UI.polymers[self.data.chain1]
            res1List = res1List + poly1.getResNamesIDsList()
        else:
            res1State = "disabled"
        self.widRes1 = ComboBox(self.listFrame, res1List, self.data.getRes1NameId(), 
                                 data="res1", command=self.changeAttribute)
        self.widRes1.config(state=res1State)
        self.widRes1.grid(row=self.guiID, column=colNum.next(),sticky=Tk.W+Tk.E)

        # Atom 1
        atom1List = [""]
        atom1State = "normal"
        if self.data.chain1 != "" and self.data.res1 != "":
            atom1List = atom1List + self.data.getRes1AtomsList()
        elif self.data.res1 == "":
            atom1State = "disabled"
        self.widAtom1 = ComboBox(self.listFrame, atom1List, self.data.atom1, 
                                 data="atom1", command=self.changeAttribute)
        self.widAtom1.config(state=atom1State)
        self.widAtom1.grid(row=self.guiID, column=colNum.next(),sticky=Tk.W+Tk.E)

        # Weld label
        self.widWeldLabel = Tk.Label(self.listFrame, text="Attached to")
        self.widWeldLabel.grid(row=self.guiID, column=colNum.next(),sticky=Tk.W+Tk.E+Tk.S+Tk.N)

        # Chain 2
        self.widChain2 = ComboBox(self.listFrame, chainsList, self.data.chain2, 
                                 data="chain2", command=self.changeAttribute)
        self.widChain2.grid(row=self.guiID, column=colNum.next(),sticky=Tk.W+Tk.E)

        # Residue 2
        res2List = [""]
        poly2 = MMB_UI.BioPolymer()
        res2State = "normal"
        if self.data.chain2 != "": 
            poly2 = MMB_UI.polymers[self.data.chain2]
            res2List = res2List + poly2.getResNamesIDsList()
        else:
            res2State = "disabled"
        self.widRes2 = ComboBox(self.listFrame, res2List, self.data.getRes2NameId(), 
                                 data="res2", command=self.changeAttribute)
        self.widRes2.config(state=res2State)
        self.widRes2.grid(row=self.guiID, column=colNum.next(),sticky=Tk.W+Tk.E)

        # Atom 2
        atom2List = [""]
        atom2State = "normal"
        if self.data.chain2 != "" and self.data.res2 != "":
            atom2List = atom2List + self.data.getRes2AtomsList()
        elif self.data.res2 == "":
            atom2State = "disabled"
        self.widAtom2 = ComboBox(self.listFrame, atom2List, self.data.atom2, 
                                 data="atom2", command=self.changeAttribute)
        self.widAtom2.config(state=atom2State)
        self.widAtom2.grid(row=self.guiID, column=colNum.next(),sticky=Tk.W+Tk.E)

        self.widTether = Tk.Checkbutton(self.listFrame, text="", variable=self.varTether, command=self.changeTether)
        self.widTether.grid(row=self.guiID, column=colNum.next(), sticky=Tk.W+Tk.E+Tk.S+Tk.N)

        if not self.widForceConstant:
            self.widForceConstant = NumericEntry(self.listFrame,
                                          value=self.data.forceConstant,
                                          command=self.changeConstant, 
                                          width=5)
            if self.data.new: self.widForceConstant.config(selectbackground="dodger blue")
            self.widForceConstant.grid(row=self.guiID, column=colNum.next(),sticky=Tk.W+Tk.E+Tk.S+Tk.N)
        else:
            colNum.next()

        if not self.widDeadLength:
            self.widDeadLength = NumericEntry(self.listFrame,
                                          value=self.data.deadLength*10.0,
                                          command=self.changeDeadLength, 
                                          width=5)
            if self.data.new: self.widDeadLength.config(selectbackground="dodger blue")
            self.widDeadLength.grid(row=self.guiID, column=colNum.next(),sticky=Tk.W+Tk.E+Tk.S+Tk.N)
        else:
            colNum.next()

        self.drawValidNewButtons(colNum.next())

        self.widSelect = Tk.Checkbutton(self.listFrame, text="", variable=self.varSelect, command=self.selectLine)
        self.widSelect.grid(row=self.guiID, column=colNum.next(),sticky=Tk.W+Tk.E+Tk.S+Tk.N)

        self.updateBgColor()
        self.applyTags()


#############################################################################
class DensityWidgetsLine(WidgetsLine):
    def __init__(self, listFrame, guiID, density, select=False, addCommand=None):
        WidgetsLine.__init__(self, listFrame, guiID, density)
        self.addCommand = addCommand

    def display(self):
        # clear my widgets
        [w.destroy() for w in [self.__dict__[k] for k in self.__dict__ if k.startswith("wid")] if w]
        colNum = self.iter1()

        self.widLabel = Tk.Label(self.listFrame, text=str(self.guiID))
        self.widLabel.grid(row=self.guiID, column=colNum.next(),sticky=Tk.W+Tk.E+Tk.S+Tk.N)

        self.widDensityLabel = Tk.Label(self.listFrame, text="Density")
        self.widDensityLabel.grid(row=self.guiID, column=colNum.next(),sticky=Tk.W+Tk.E+Tk.S+Tk.N)

        chainsList = MMB_UI.sortedChainsIDs()
        if self.data.new: chainsList.insert(0,"")
        self.widChain = ComboBox(self.listFrame, chainsList, self.data.chainID, 
                                 data="chainID", command=self.changeAttribute)
        self.widChain.grid(row=self.guiID, column=colNum.next(),sticky=Tk.W+Tk.E)

        resNamesIdsList = []
        if self.widChain.get() != "":
            poly = MMB_UI.polymers[self.data.chainID]
            resNamesIdsList = poly.getResNamesIDsList()
            resNamesIdsList.insert(0,"All")
            self.widStartRes = ComboBox(self.listFrame, resNamesIdsList,  self.data.getResStartNameId(), 
                                     data="resStart", command=self.changeAttribute)

            if self.widStartRes.get() != "All":
                resNamesIdsList = poly.getResNamesIDsList(self.data.resStart)
                self.widEndRes = ComboBox(self.listFrame, resNamesIdsList,  self.data.getResEndNameId(), 
                                         data="resEnd", command=self.changeAttribute)
            else:
                self.widEndRes = EmptyLabel(self.listFrame)
        else:
            self.widStartRes = EmptyLabel(self.listFrame)
            self.widEndRes   = EmptyLabel(self.listFrame)

        self.widStartRes.grid(row=self.guiID, column=colNum.next(),sticky=Tk.W+Tk.E+Tk.S+Tk.N)
        self.widEndRes.grid(row=self.guiID, column=colNum.next(),sticky=Tk.W+Tk.E+Tk.S+Tk.N)

        self.drawValidNewButtons(colNum.next())


        self.widSelect = Tk.Checkbutton(self.listFrame, text="", variable=self.varSelect, command=self.selectLine)
        self.widSelect.grid(row=self.guiID, column=colNum.next(),sticky=Tk.W+Tk.E+Tk.S+Tk.N)

        self.updateBgColor()
        self.applyTags()


#############################################################################
class ScrolledFrame(Tk.Frame):
    def __init__(self, parent, tkTags=()):
        self.tkTags = tkTags

        self.canvas = Tk.Canvas(parent)
        self.canvas.bindtags(self.canvas.bindtags() + self.tkTags)

        Tk.Frame.__init__(self, self.canvas)
        self.bindtags(self.bindtags() + self.tkTags)

        self.canvas.create_window((0,0), window=self, anchor="nw", tags="frame")

        self.vsb = ttk.Scrollbar(parent, orient="vertical", command=self.canvas.yview)
        self.canvas.configure(yscrollcommand=self.vsb.set)
        self.vsb.pack(side="right", fill="y")
        self.canvas.pack(fill="both", expand=True)
        self.canvas.xview_moveto(0)
        self.canvas.yview_moveto(0)

        self.update_idletasks()
        self.bind("<Configure>", self.OnScrolledFrameConfigure)
        self.canvas.bind_class(self.tkTags, "<Button-4>",  self.OnScrolledFrameMouseWheel)
        self.canvas.bind_class(self.tkTags, "<Button-5>",  self.OnScrolledFrameMouseWheel)
        self.canvas.bind_class(self.tkTags, "<MouseWheel>",self.OnScrolledFrameMouseWheel)
        
        
    ## Callback when resizing a frame within a self.canvas - sets the scrollbar size
    def OnScrolledFrameConfigure(self, event):
        w,h = event.width, event.height
        wnatural = self.winfo_reqwidth()
        self.canvas.configure(width= w if w>wnatural else wnatural)
        hnatural = self.winfo_reqheight()
        # self.canvas.configure(height= h if h>hnatural else hnatural)
        self.canvas.configure(scrollregion=[0,0,wnatural,hnatural])

        # sbtop, sbbottom = self.vsb.get()
        # if sbbottom - sbtop < 1.0:
        #     self.canvas.bind_class(self.tkTags, "<Button-4>",  self.OnScrolledFrameMouseWheel)
        #     self.canvas.bind_class(self.tkTags, "<Button-5>",  self.OnScrolledFrameMouseWheel)
        #     self.canvas.bind_class(self.tkTags, "<MouseWheel>",self.OnScrolledFrameMouseWheel)
        # else:
        #     self.canvas.unbind_class(self.tkTags, "<Button-4>" )
        #     self.canvas.unbind_class(self.tkTags, "<Button-5>" )
        #     self.canvas.unbind_class(self.tkTags, "<MouseWheel>")


    ## Callback for mouse wheel events in ScrolledFrame
    def OnScrolledFrameMouseWheel(self, event):
        # print event.delta
        sbtop, sbbottom = self.vsb.get()
        if sbbottom - sbtop < 1.0:
            self.canvas.yview_scroll(-event.delta,"units")

#############################################################################
def createCustomTextEditors(parent, name):
    cmdGlobalLabel = Tk.Label(parent, text=name)
    cmdGlobalLabel.pack()
    cmdEditor = Tk.Text(parent)

    scroll = Tk.Scrollbar(cmdEditor)
    scroll.pack(side='right', fill='y')
    scroll.config(command=cmdEditor.yview)
    cmdEditor.config(yscrollcommand = scroll.set)

    return cmdEditor

#############################################################################
class CmdPrompt(Tk.Frame):
    def __init__(self, parent, history = [], value=None, command=None, label=None, *args, **kw):
        Tk.Frame.__init__(self,parent)

        self.command = command
        self.textVar = Tk.StringVar()
        if value: self.textVar.set(str(value))
        
        if label:
            self.label = Tk.Label(self, text=label)
            self.label.pack(side="left")

        self.entry = Tk.Entry(self, 
                          textvariable=self.textVar, 
                          *args, **kw)
        self.entry.pack(side="left",fill="x",expand=True)
        self.entry.bind("<Return>", self.onEnterKey)
        self.entry.bind("<KP_Enter>", self.onEnterKey)
        self.entry.bind("<Key-Up>", self.onUpKey)
        self.entry.bind("<Key-Down>", self.onDownKey)

        self.button = Tk.Button(self, text="Send", command=lambda event=None: self.onEnterKey(event))
        self.button.pack(side="right")

        self.history = history
        self.historyPos = len(self.history)

    def onEnterKey(self, event):
        if self.textVar.get() != "":
            if self.command:
                self.command(self.textVar.get())
            self.history.append(self.textVar.get())
            self.historyPos = len(self.history)
            self.entry.selection_range(0, Tk.END)

    def onUpKey(self, event):
        self.historyPos -= 1
        if self.historyPos <= 0: 
            self.historyPos = 1
            return
        if len(self.history) > 0:
            self.textVar.set(self.history[self.historyPos-1])
            self.entry.selection_range(0, Tk.END)

    def onDownKey(self, event):
        self.historyPos += 1
        if self.historyPos > len(self.history): 
            self.historyPos = len(self.history)+1
            self.textVar.set("")
            return 
        if len(self.history) > 0:
            self.textVar.set(self.history[self.historyPos-1])
            self.entry.selection_range(0, Tk.END)

#############################################################################
class LabeledButton(Tk.Frame):
    def __init__(self, parent, history = [], command=None, labelText=" ", buttonText=" " , *args, **kw):
        Tk.Frame.__init__(self,parent)
        self.command = command
        self.label = Tk.Label(self, text=labelText)
        self.label.pack(side="left")
        self.button = Tix.Button(self, text=buttonText,command=command)
        self.button.pack(side="left")

    def setLabelText(self, text):
        self.label.configure(text=text)

#############################################################################
## Enable/Disable a widget according to a boolean
def changeWidgetState(widget, boolean):
    if boolean:
        widget.configure(state="normal")
    else:
        widget.configure(state="disabled")

#############################################################################
from OpenSave import OpenModeless
#* Open file chimera Style
class InputPdbDialog(OpenModeless):
    default = 'Open Pdb File'
    title = 'MMB -- Polymers'
    def OpenPdbFile(self, *args):
        self.command(self.getPaths()[0])
        self.destroy()

class InputMMBFileDialog(OpenModeless):
    default = 'Open MMB File'
    title = 'MMB -- Command file'
    def OpenMMBFile(self, *args):
        self.command(self.getPaths()[0])
        self.destroy()

class OpenXplorMapDialog(OpenModeless):
    default = 'Open Map'
    def __init__(self, command=None, multiple=False):
        title = 'MMB -- Xplor Density Map'
        OpenModeless.__init__(self, title = title, 
                              filters=[("Xplor density maps","*.xplor")], 
                              defaultFilter="Xplor density maps",
                              command=command, multiple=multiple)

    def OpenMap(self, *args):
        self.command(self.getPaths()[0])
        self.destroy()




