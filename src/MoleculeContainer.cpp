

#include "MoleculeContainer.h"
#include <fstream>

void CompoundObjectMapContainer::loadCompoundMap() {
    cout <<__FILE__<<":"<<__LINE__<<" starting loadCompoundMap() "<<endl;
    compoundMap.insert (std::pair <const String , Compound> ("PurineBaseCore",PurineBaseCore() ));
    RibonucleotideResidue::Guanylate myGuanylate; myGuanylate.assignBiotypes();
    compoundMap.insert (std::pair <const String , Compound> ("Guanylate",myGuanylate ));
    //compoundMap.insert (std::pair <const String , Compound> ("Guanylate",RibonucleotideResidue::Guanylate() ));
    // from molmodel/include/molmodel/internal/NA.h:
    compoundMap.insert (std::pair <const String , Compound> ("NaPhosphodiesterLinkage",NaPhosphodiesterLinkage("NaPhosphodiesterLinkage") ));
    compoundMap.insert (std::pair <const String , Compound> ("FivePrimeNaPhosphateGroup",FivePrimeNaPhosphateGroup("FivePrimeNaPhosphateGroup") ));
    compoundMap.insert (std::pair <const String , Compound> ("ThreePrimeNaPhosphateGroup",ThreePrimeNaPhosphateGroup("ThreePrimeNaPhosphateGroup") ));
    compoundMap.insert (std::pair <const String , Compound> ("RibonucleosideResidue",RibonucleosideResidue("RibonucleosideResidue") ));
    compoundMap.insert (std::pair <const String , Compound> ("MethyleneGroup",MethyleneGroup() ));
    compoundMap.insert (std::pair <const String , Compound> ("methyl",MethylGroup() ));
    compoundMap.insert (std::pair <const String , Compound> ("MethylGroup",MethylGroup() ));
    compoundMap.insert (std::pair <const String , Compound> ("Methane",Methane()     ));
    compoundMap.insert (std::pair <const String , Compound> ("MagnesiumIon",MagnesiumIon()     ));
    compoundMap.insert (std::pair <const String , Compound> ("AromaticSixMemberedCHGroup",AromaticSixMemberedCHGroup() ));
    compoundMap.insert (std::pair <const String , Compound> ("AliphaticHydrogen",AliphaticHydrogen() ));

// From Compound.h:
// SingleAtom
    /*compoundMap.insert (std::pair <const String , Compound> ("UnivalentAtom",
    compoundMap.insert (std::pair <const String , Compound> ("BivalentAtom",
    compoundMap.insert (std::pair <const String , Compound> ("TrivalentAtom",
    compoundMap.insert (std::pair <const String , Compound> ("QuadrivalentAtom",*/
    compoundMap.insert (std::pair <const String , Compound> ("AlcoholOHGroup",AlcoholOHGroup()       ));
    compoundMap.insert (std::pair <const String , Compound> ("PrimaryAmineGroup",PrimaryAmineGroup() ));
    compoundMap.insert (std::pair <const String , Compound> ("CarboxylateGroup",CarboxylateGroup()   ));


    cout <<__FILE__<<":"<<__LINE__<<" done with loadCompoundMap() "<<endl;
}

void CompoundObjectMapContainer::printCompoundMap() {
    cout <<__FILE__<<":"<<__LINE__<<" Available compounds are: "<<endl;
    map <const String , Compound>::iterator compoundMapIterator = compoundMap.begin();
    for (compoundMapIterator = compoundMap.begin(); compoundMapIterator != compoundMap.end(); compoundMapIterator++) {
        std::cout<<compoundMapIterator->first<<std::endl;
    }
}

/*void CompoundObjectMapContainer::loadSingleAtomMap() {
    cout <<__FILE__<<":"<<__LINE__<<" starting loadSingleAtomMap() "<<endl;
    singleAtomMap.insert (std::pair <const String , Compound::SingleAtom> ("AliphaticHydrogen",AliphaticHydrogen() ));
    singleAtomMap.insert (std::pair <const String , Compound::SingleAtom> ("AliphaticCarbon",AliphaticCarbon() ));
    singleAtomMap.insert (std::pair <const String , Compound::SingleAtom> ("UnivalentAtom",UnivalentAtom() ));
    singleAtomMap.insert (std::pair <const String , Compound::SingleAtom> ("BivalentAtom",BivalentAtom() ));
    singleAtomMap.insert (std::pair <const String , Compound::SingleAtom> ("TrivalentAtom",TrivalentAtom() ));
    singleAtomMap.insert (std::pair <const String , Compound::SingleAtom> ("QuadrivalentAtom",QuadrivalentAtom() ));

    cout <<__FILE__<<":"<<__LINE__<<" done with loadSingleAtomMap() "<<endl;
}*/
void CompoundObjectMapContainer::loadBiotypeMap() {
    cout <<__FILE__<<":"<<__LINE__<<" starting loadBiotypeMap() "<<endl;
    biotypeMap.insert (std::pair <const String , Biotype> ("MethaneH",Biotype::MethaneH() ));
    biotypeMap.insert (std::pair <const String , Biotype> ("MethaneC",Biotype::MethaneC() ));
    biotypeMap.insert (std::pair <const String , Biotype> ("SerineN", Biotype::SerineN()  ));
    //biotypeMap.insert (std::pair <const String , Biotype> ("Magnesium Ion", Biotype::MagnesiumIon()  ));
    cout <<__FILE__<<":"<<__LINE__<<" done with loadBiotypeMap() "<<endl;
}

Compound CompoundObjectMapContainer::fetchCompound(const String compoundName)  {
    //printCompoundMap();
    if (compoundMap.find(compoundName) == compoundMap.end()) {
        ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" No Compound object found with name "<<compoundName <<endl; ErrorManager::instance.treatError();
    }
    else return compoundMap[compoundName];
}


Biotype CompoundObjectMapContainer::fetchBiotype   (const String compoundName)  {
    cout<<__FILE__<<":"<<__LINE__<<" You have requested biotype named : >"<<compoundName<<"< "<<endl;
    if (biotypeMap.find(compoundName) == biotypeMap.end()) {
        ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" No Biotype  object found with name >"<<compoundName <<"< "<<endl; 
        ErrorManager::instance.treatError();
    }
    else {
        Biotype myBiotype = biotypeMap[compoundName];
        BiotypeIndex myBiotypeIndex = myBiotype.getIndex();
        if (!(myBiotype.exists(myBiotypeIndex))) {
            ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" This biotype doesn't appear to exist! >"<<compoundName <<"<, >"<<myBiotypeIndex<<endl; 
            ErrorManager::instance.treatError();
        } 
        cout<<__FILE__<<":"<<__LINE__<<" Returning biotype with index : >"<<myBiotypeIndex<<"< "<<endl;
        return myBiotype;                
    }
}

const Element CompoundObjectMapContainer::fetchElement (String elementName) {
    cout<<__FILE__<<":"<<__LINE__<<" About to fetch element with name "<<elementName<<endl;
    if (elementName.compare("Hydrogen") == 0) {cout<<__FILE__<<":"<<__LINE__<< ( Element::Hydrogen()).getName()<<endl; return Element::         Hydrogen();}
    //cout<<__FILE__<<":"<<__LINE__<<endl;
    if (elementName.compare("Deuterium") == 0) return Element::        Deuterium();
    if (elementName.compare("Helium") == 0) return Element::           Helium();
    if (elementName.compare("Lithium") == 0) return Element::          Lithium();
    if (elementName.compare("Beryllium") == 0) return Element::        Beryllium();
    if (elementName.compare("Boron") == 0) return Element::            Boron();
    if (elementName.compare("Carbon") == 0) return Element::           Carbon();
    if (elementName.compare("Nitrogen") == 0) return Element::         Nitrogen();
    if (elementName.compare("Oxygen") == 0) return Element::           Oxygen();
    if (elementName.compare("Fluorine") == 0) return Element::         Fluorine();
    if (elementName.compare("Neon") == 0) return Element::             Neon();
    if (elementName.compare("Sodium") == 0) return Element::           Sodium();
    if (elementName.compare("Magnesium") == 0) return Element::        Magnesium();
    if (elementName.compare("Aluminum") == 0) return Element::         Aluminum();
    if (elementName.compare("Silicon") == 0) return Element::          Silicon();
    cout<<__FILE__<<":"<<__LINE__<<endl;
    if (elementName.compare("Phosphorus") == 0) return Element::       Phosphorus();
    cout<<__FILE__<<":"<<__LINE__<<endl;
    if (elementName.compare("Sulfur") == 0) return Element::           Sulfur();
    if (elementName.compare("Chlorine") == 0) return Element::         Chlorine();
    if (elementName.compare("Argon") == 0) return Element::            Argon();
    if (elementName.compare("Potassium") == 0) return Element::        Potassium();
    if (elementName.compare("Calcium") == 0) return Element::          Calcium();
    if (elementName.compare("Scandium") == 0) return Element::         Scandium();
    if (elementName.compare("Titanium") == 0) return Element::         Titanium();
    if (elementName.compare("Vanadium") == 0) return Element::         Vanadium();
    if (elementName.compare("Chromium") == 0) return Element::         Chromium();
    if (elementName.compare("Manganese") == 0) return Element::        Manganese();
    if (elementName.compare("Iron") == 0) return Element::             Iron();
    if (elementName.compare("Cobalt") == 0) return Element::           Cobalt();
    if (elementName.compare("Nickel") == 0) return Element::           Nickel();
    if (elementName.compare("Copper") == 0) return Element::           Copper();
    if (elementName.compare("Zinc") == 0) return Element::             Zinc();
    if (elementName.compare("Gallium") == 0) return Element::          Gallium();
    if (elementName.compare("Germanium") == 0) return Element::        Germanium();
    if (elementName.compare("Arsenic") == 0) return Element::          Arsenic();
    if (elementName.compare("Selenium") == 0) return Element::         Selenium();
    if (elementName.compare("Bromine") == 0) return Element::          Bromine();
    if (elementName.compare("Krypton") == 0) return Element::          Krypton();
    if (elementName.compare("Rubidium") == 0) return Element::         Rubidium();
    if (elementName.compare("Strontium") == 0) return Element::        Strontium();
    if (elementName.compare("Yttrium") == 0) return Element::          Yttrium();
    if (elementName.compare("Zirconium") == 0) return Element::        Zirconium();
    if (elementName.compare("Niobium") == 0) return Element::          Niobium();
    if (elementName.compare("Molybdenum") == 0) return Element::       Molybdenum();
    if (elementName.compare("Technetium") == 0) return Element::       Technetium();
    if (elementName.compare("Ruthenium") == 0) return Element::        Ruthenium();
    if (elementName.compare("Rhodium") == 0) return Element::          Rhodium();
    if (elementName.compare("Palladium") == 0) return Element::        Palladium();
    if (elementName.compare("Silver") == 0) return Element::           Silver();
    if (elementName.compare("Cadmium") == 0) return Element::          Cadmium();
    if (elementName.compare("Indium") == 0) return Element::           Indium();
    if (elementName.compare("Tin") == 0) return Element::              Tin();
    if (elementName.compare("Antimony") == 0) return Element::         Antimony();
    if (elementName.compare("Tellurium") == 0) return Element::        Tellurium();
    if (elementName.compare("Iodine") == 0) return Element::           Iodine();
    if (elementName.compare("Xenon") == 0) return Element::            Xenon();
    if (elementName.compare("Cesium") == 0) return Element::           Cesium();
    if (elementName.compare("Barium") == 0) return Element::           Barium();
    if (elementName.compare("Lanthanum") == 0) return Element::        Lanthanum();
    if (elementName.compare("Cerium") == 0) return Element::           Cerium();
    if (elementName.compare("Praseodymium") == 0) return Element::     Praseodymium();
    if (elementName.compare("Neodymium") == 0) return Element::        Neodymium();
    if (elementName.compare("Promethium") == 0) return Element::       Promethium();
    if (elementName.compare("Samarium") == 0) return Element::         Samarium();
    if (elementName.compare("Europium") == 0) return Element::         Europium();
    if (elementName.compare("Gadolinium") == 0) return Element::       Gadolinium();
    if (elementName.compare("Terbium") == 0) return Element::          Terbium();
    if (elementName.compare("Dysprosium") == 0) return Element::       Dysprosium();
    if (elementName.compare("Holmium") == 0) return Element::          Holmium();
    if (elementName.compare("Erbium") == 0) return Element::           Erbium();
    if (elementName.compare("Thulium") == 0) return Element::          Thulium();
    if (elementName.compare("Ytterbium") == 0) return Element::        Ytterbium();
    if (elementName.compare("Lutetium") == 0) return Element::         Lutetium();
    if (elementName.compare("Hafnium") == 0) return Element::          Hafnium();
    if (elementName.compare("Tantalum") == 0) return Element::         Tantalum();
    if (elementName.compare("Tungsten") == 0) return Element::         Tungsten();
    if (elementName.compare("Rhenium") == 0) return Element::          Rhenium();
    if (elementName.compare("Osmium") == 0) return Element::           Osmium();
    if (elementName.compare("Iridium") == 0) return Element::          Iridium();
    if (elementName.compare("Platinum") == 0) return Element::         Platinum();
    if (elementName.compare("Gold") == 0) return Element::             Gold();
    if (elementName.compare("Mercury") == 0) return Element::          Mercury();
    if (elementName.compare("Thallium") == 0) return Element::         Thallium();
    if (elementName.compare("Lead") == 0) return Element::             Lead();
    if (elementName.compare("Bismuth") == 0) return Element::          Bismuth();
    if (elementName.compare("Polonium") == 0) return Element::         Polonium();
    if (elementName.compare("Astatine") == 0) return Element::         Astatine();
    if (elementName.compare("Radon") == 0) return Element::            Radon();
    if (elementName.compare("Francium") == 0) return Element::         Francium();
    if (elementName.compare("Radium") == 0) return Element::           Radium();
    if (elementName.compare("Actinium") == 0) return Element::         Actinium();
    if (elementName.compare("Thorium") == 0) return Element::          Thorium();
    if (elementName.compare("Protactinium") == 0) return Element::     Protactinium();
    if (elementName.compare("Uranium") == 0) return Element::          Uranium();
    if (elementName.compare("Neptunium") == 0) return Element::        Neptunium();
    if (elementName.compare("Plutonium") == 0) return Element::        Plutonium();
    if (elementName.compare("Americium") == 0) return Element::        Americium();
    if (elementName.compare("Curium") == 0) return Element::           Curium();
    if (elementName.compare("Berkelium") == 0) return Element::        Berkelium();
    if (elementName.compare("Californium") == 0) return Element::      Californium();
    if (elementName.compare("Einsteinium") == 0) return Element::      Einsteinium();
    if (elementName.compare("Fermium") == 0) return Element::          Fermium();
    if (elementName.compare("Mendelevium") == 0) return Element::      Mendelevium();
    if (elementName.compare("Nobelium") == 0) return Element::         Nobelium();
    if (elementName.compare("Lawrencium") == 0) return Element::       Lawrencium();
    if (elementName.compare("Rutherfordium") == 0) return Element::    Rutherfordium();
    if (elementName.compare("Dubnium") == 0) return Element::          Dubnium();
    if (elementName.compare("Seaborgium") == 0) return Element::       Seaborgium();
    if (elementName.compare("Bohrium") == 0) return Element::          Bohrium();
    if (elementName.compare("Hassium") == 0) return Element::          Hassium();
    if (elementName.compare("Meitnerium") == 0) return Element::       Meitnerium();
    if (elementName.compare("Darmstadtium") == 0) return Element::     Darmstadtium();
    if (elementName.compare("Roentgenium") == 0) return Element::      Roentgenium();
    if (elementName.compare("Ununbium") == 0) return Element::         Ununbium();
    if (elementName.compare("Ununtrium") == 0) return Element::        Ununtrium();
    if (elementName.compare("Ununquadium") == 0) return Element::      Ununquadium();
    if (elementName.compare("Ununpentium") == 0) return Element::      Ununpentium();
    if (elementName.compare("Ununhexium") == 0) return Element::       Ununhexium();
    // If we didn't return after the above, something is wrong:
    ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" Element "<<elementName<<" not found! Did you use proper capitalization (e.g. Carbon)?"<<endl;
    ErrorManager::instance.treatError();       
    //return Element::Hydrogen();
}




Compound::SingleAtom  CompoundObjectMapContainer::fetchSingleAtom(const String className, Compound::AtomName& atomName , String elementName, Angle angle1 = 180*Deg2Rad )  {
    cout<<__FILE__<<":"<<__LINE__<<" Fetching single atom of class >"<<className<<"< "<<endl;
    const Element myElement = fetchElement(elementName);
    cout<<__FILE__<<":"<<__LINE__<<" using element : >"<<myElement.getName()<<"< "<<endl;
    cout<<__FILE__<<":"<<__LINE__<<endl;
    if (className.compare("AliphaticHydrogen") == 0) {return AliphaticHydrogen(atomName);}
    if (className.compare("AliphaticCarbon") == 0) {return AliphaticCarbon(atomName);}
    if (className.compare("UnivalentAtom") == 0) {
	cout<<__FILE__<<":"<<__LINE__<<" using element : >"<<myElement.getName()<<"< "<<endl;
        return UnivalentAtom(atomName,myElement);
        //return UnivalentAtom(atomName,   myElement);
	cout<<__FILE__<<":"<<__LINE__<<endl;
    }
    if (className.compare("BivalentAtom") == 0) {return BivalentAtom(atomName,fetchElement(elementName) , angle1);} // the angle should be adjusted later
    if (className.compare("TrivalentAtom") == 0) {return TrivalentAtom(atomName,fetchElement(elementName));}
    if (className.compare("QuadrivalentAtom") == 0) {return QuadrivalentAtom(atomName,fetchElement(elementName));}
    // If we didn't return after the above, something is wrong:
    ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" The SingleAtom subclass "<<className<<" was not found! Accepted subclasses are UnivalentAtom, BivalentAtom, BivalentAtom, TrivalentAtom, QuadrivalentAtom"<<endl;
    ErrorManager::instance.treatError();       
}

CustomMolecule::CustomMolecule(vector <vector <String> > moleculeBuildCommandVector, DuMMForceFieldSubsystem & dumm) {
    CompoundObjectMapContainer compoundObjectMapContainer;
    if (moleculeBuildCommandVector.size() < 1) {
        ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" setBaseCompound must be the first command you issue!"<<endl;
        ErrorManager::instance.treatError();       
    }
    for (int i = 0; i < moleculeBuildCommandVector.size(); i ++) {
        cout<<__FILE__<<":"<<__LINE__<<" Contents of moleculeBuildCommandVector[i] : ("<<moleculeBuildCommandVector[i].size()<<" elements) " <<endl;;
        for (int j = 0; j < moleculeBuildCommandVector[i].size(); j++) {
            cout<<">"<<moleculeBuildCommandVector[i][j]<<"< ";
        }
        cout<<"."<<endl;
	if ((moleculeBuildCommandVector[i])[0].compare("setBaseCompound") == 0) {  // element 0 of moleculeBuildCommand is always a command
            //cout<<__FILE__<<":"<<__LINE__<<" Contents of moleculeBuildCommandVector[i] : >"<<moleculeBuildCommandVector[i][0]<<"<, >"<<moleculeBuildCommandVector[i][1]<<"< ."<<endl; ;
            compoundObjectMapContainer.printCompoundMap();
	    if (i > 0) {ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" setBaseCompound must be the first command you issue!"<<endl; ErrorManager::instance.treatError();}
	    if ((moleculeBuildCommandVector[i]).size() != 2) {ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" Wrong number of parameters!"<<endl; ErrorManager::instance.treatError();}
	    //if ((moleculeBuildCommandVector[i]).size() > 2) {ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" Too many parameters!"<<endl; ErrorManager::instance.treatError();}
            cout<<__FILE__<<":"<<__LINE__<<" About to fetchCompound ("<<moleculeBuildCommandVector[i][1]<<") "<<endl;
	    Compound myCompound = compoundObjectMapContainer.fetchCompound(moleculeBuildCommandVector[i][1]);
	    setBaseCompound(moleculeBuildCommandVector[i][1],                                    // name of the Compound (e.g. "MethyleneGroup") is element [1] of moleculeBuildCommand
		myCompound               //followed by the corresponding Compound object
		); 
	    inheritAtomNames(moleculeBuildCommandVector[i][1]);                                  // for now I am assuming we will always want to inherit atom names of this compound.
	}
	else if ((moleculeBuildCommandVector[i])[0].compare("setBaseAtom") == 0) {  // element 0 of moleculeBuildCommand is always a command
            cout<<__FILE__<<":"<<__LINE__<<" Syntax: setBaseAtom <object type> <atom name> <element name> [<angle1>] "<<endl;
	    if (i > 0) {ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" If used, this must be the first command you issue!"<<endl; ErrorManager::instance.treatError();}
	    if ((moleculeBuildCommandVector[i]).size() > 5) {ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" Too many parameters!"<<endl; ErrorManager::instance.treatError();}
	    String objectType = moleculeBuildCommandVector[i][1]; 
	    String atomName = moleculeBuildCommandVector[i][2]; 
	    String elementname = moleculeBuildCommandVector[i][3]; 
            Angle myAngle1 = 180*Deg2Rad;
            if (moleculeBuildCommandVector[i].size() > 4)
                myAngle1 = atof(moleculeBuildCommandVector[i][4].c_str());

            cout<<__FILE__<<":"<<__LINE__<<endl;
	    setBaseAtom    (
		compoundObjectMapContainer.fetchSingleAtom(objectType, atomName, elementname , myAngle1)               //followed by the corresponding Compound object
		); 
            cout<<__FILE__<<":"<<__LINE__<<endl;
	    //inheritAtomNames(moleculeBuildCommandVector[i][1]);                                  // for now I am assuming we will always want to inherit atom names of this compound.
	}

	else if ((moleculeBuildCommandVector[i])[0].compare("convertInboardBondCenterToOutboard") == 0 ) {
            //cout<<__FILE__<<":"<<__LINE__<<" Contents of moleculeBuildCommandVector[i] : >"<<moleculeBuildCommandVector[i][0]<<"<"<<endl;
	    if (i == 0) {ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" This command cannot be the first to issue.  You first need to setBaseCompound!"<<endl; ErrorManager::instance.treatError();}
	    convertInboardBondCenterToOutboard();}

	else if ((moleculeBuildCommandVector[i])[0].compare("bondCompound") == 0 ) {
            cout<<__FILE__<<":"<<__LINE__<<" Syntax : bondCompound <name of added compound in parent> <compound to add> <name of parent bond at which to attach e.g. methyl/bond >  "<<endl;
            /*cout<<__FILE__<<":"<<__LINE__<<" Contents of moleculeBuildCommandVector[i] : "<<endl;;
            for (int j = 0; j < moleculeBuildCommandVector[i].size(); j++) {
                cout<<">"<<moleculeBuildCommandVector[i][j]<<"< ";
            }
            cout<<endl;*/
            if (moleculeBuildCommandVector[i].size  () >4){
	        ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" Too many parameters!"<<endl; ErrorManager::instance.treatError();
            }
            String subcompoundNameInParent = moleculeBuildCommandVector[i][1];
            String compoundToAdd = moleculeBuildCommandVector[i][2];
            String bondName = moleculeBuildCommandVector[i][3];
            Compound myCompound = compoundObjectMapContainer.fetchCompound(compoundToAdd);
	    bondCompound(
                subcompoundNameInParent,
                myCompound,
		bondName,  // name of bond at which to attach the atom    
                .14 //hack, fix later
		); 
	    // For example:
            // bondCompound("methyl2", MethylGroup(), "methyl1/bond");
            
        }
	else if ((moleculeBuildCommandVector[i])[0].compare("bondAtom") == 0 ) {
            cout<<__FILE__<<":"<<__LINE__<<" Syntax : bondAtom <molmodel class of added atom (AliphaticHydrogen, UnivalentAtom, DivalentAtom, etc> <name of atom to be added, e.g. H1> <name of element to be added e.g. Hydrogen > <name of bond at which to attach atom e.g. methyl/bond > <bond length> [<dihedral angle, degrees>] [bond mobility: Default, Free, Torsion, or Rigid]"<<endl;
            if (moleculeBuildCommandVector[i].size  () <6){
	        ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" Too few parameters!"<<endl; ErrorManager::instance.treatError();
            }
            if (moleculeBuildCommandVector[i].size  () >8){
	        ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" Too many parameters!"<<endl; ErrorManager::instance.treatError();
            }
            String singleAtomSubclass = moleculeBuildCommandVector[i][1];
            String addedAtomName = moleculeBuildCommandVector[i][2];
            String addedElementName = moleculeBuildCommandVector[i][3]; 
            String bondName = moleculeBuildCommandVector[i][4];
            double bondLength = atof(moleculeBuildCommandVector[i][5].c_str()) ;
	    
            Angle myAngle1 = 180*Deg2Rad; // this is a bond angle 
            //if (moleculeBuildCommandVector[i].size() > 6)
            //    myAngle1 = atof(moleculeBuildCommandVector[i][6].c_str());
            Compound::SingleAtom addedAtom = compoundObjectMapContainer.fetchSingleAtom(singleAtomSubclass,addedAtomName,addedElementName );

            Angle myDihedral = 180*Deg2Rad;
            if (moleculeBuildCommandVector[i].size() > 6)
                myDihedral = atof(moleculeBuildCommandVector[i][6].c_str())*Deg2Rad ;

            String myBondMobilityString = "Torsion";
            if (moleculeBuildCommandVector[i].size() > 7) {
                myBondMobilityString = moleculeBuildCommandVector[i][7];}
            BondMobility::Mobility myBondMobility = stringToBondMobility(myBondMobilityString);
            
	    bondAtom(
                addedAtom ,
		bondName,  // name of bond at which to attach the atom    
		bondLength,
                myDihedral,          // This is the default dihedral angle parameter, defaults to 180 * Deg2Rad (i.e. pi)
                myBondMobility
                ); // Default bond length
	    // For example:
	    // bondAtom(AliphaticHydrogen("H4"), "methyl/bond", 0.1112);
	}

	else if ((moleculeBuildCommandVector[i])[0].compare("setBiotypeIndex") == 0 ) {
            cout<<__FILE__<<":"<<__LINE__<<" Syntax: setBiotypeIndex <specific atom name e.g. H4> <Biotype e.g. MethaneH>"<<endl;
            cout<<__FILE__<<":"<<__LINE__<<"     or: setBiotypeIndex <specific atom name e.g. OP4> <Residue e.g. Phosphate,?RNA> <generic atom name e.g. OP> <ordinality : Initial, Any, or Final> "<<endl;
            String specificAtomName = moleculeBuildCommandVector[i][1];
            if ((moleculeBuildCommandVector[i]).size() == 3) {
                String biotypeName = moleculeBuildCommandVector[i][2];
                for (int k = 0; k < biotypeName.length(); k++) {
                    if (String(biotypeName[k]).compare("?") ==0) {
                        biotypeName[k] = (String(" "))[0];
                    }
                }
		setBiotypeIndex(specificAtomName,               //name
		    compoundObjectMapContainer.fetchBiotype(biotypeName).getIndex() //moleculeBuildCommandVector[i][2]).getIndex() // BiotypeIndex
		    );
		// For example:
		// setBiotypeIndex( "H4", Biotype::MethaneH().getIndex() );
            } else if ((moleculeBuildCommandVector[i]).size() == 5) {
		String specificAtomName = moleculeBuildCommandVector[i][1];
		String residueName      = moleculeBuildCommandVector[i][2];
                for (int k = 0; k < residueName.length(); k++) {
                    if (String(residueName[k]).compare("?") ==0) {
                        residueName[k] = (String(" "))[0];
                    }
                }
		String genericAtomName  = moleculeBuildCommandVector[i][3];
                cout<<__FILE__<<":"<<__LINE__<<" specificAtomName is now : >"<<specificAtomName<<"< "<<endl;
                cout<<__FILE__<<":"<<__LINE__<<" Residue name is now : >"<<residueName<<"< "<<endl;
                cout<<__FILE__<<":"<<__LINE__<<" genericAtomName is now : >"<<genericAtomName<<"< "<<endl;
                String ordinalityString = moleculeBuildCommandVector[i][4];
                enum  	ordinalityEnum { Any = 1, Initial = 2, Final = 3 };
                Ordinality::Residue myOrdinality;
                if (ordinalityString.compare("Initial") == 0) myOrdinality =  SimTK::Ordinality::Initial;
                else if (ordinalityString.compare("Final") == 0) myOrdinality =  SimTK::Ordinality::Final;
                else if (ordinalityString.compare("Any") == 0) myOrdinality =  SimTK::Ordinality::Any;
                else {
		    ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" Invalid ordinality : >"<<ordinalityString<<"< "<<endl; 
		    ErrorManager::instance.treatError();
                }
		BiotypeIndex myBiotypeIndex =  Biotype::get(residueName, genericAtomName, 
                        myOrdinality
                        ).getIndex();                
                cout<<__FILE__<<":"<<__LINE__<<" BiotypeIndex = "<<myBiotypeIndex<<endl;
                setBiotypeIndex(specificAtomName,myBiotypeIndex);
            } else {
	        ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" Wrong number of parameters ("<<(moleculeBuildCommandVector[i]).size()<<")!"<<endl; 
                ErrorManager::instance.treatError();
            }
	}
	else if ((moleculeBuildCommandVector[i])[0].compare("defineBiotype") == 0 ) {
            cout<<__FILE__<<":"<<__LINE__<<" Syntax: defineBiotype <element symbol (e.g. O, H, C)> <valence (integer)> <residue name> <generic atom name, e.g. Oxygen>"<<endl;
            if ((moleculeBuildCommandVector[i]).size() != 5) {
		    ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" Wrong number of parameters! "<<endl; 
		    ErrorManager::instance.treatError();
                }
            Biotype::defineBiotype (
                Element::getBySymbol(moleculeBuildCommandVector[i][1]),
		atoi(moleculeBuildCommandVector[i][2].c_str()),
                moleculeBuildCommandVector[i][3] ,
                moleculeBuildCommandVector[i][4]              
            );

        }

	else if ((moleculeBuildCommandVector[i])[0].compare("setCompoundName") == 0 ) {
            //cout<<__FILE__<<":"<<__LINE__<<" Contents of moleculeBuildCommandVector[i] : >"<<moleculeBuildCommandVector[i][0]<<"<, >"<<moleculeBuildCommandVector[i][1]<<"<, >"<<moleculeBuildCommandVector[i][2]<<"< ."<<endl;
	    if (moleculeBuildCommandVector[i][1].length() == 0 ) {ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" Length of Compound name must be > 0"<<endl; ErrorManager::instance.treatError();}
        setCompoundName(moleculeBuildCommandVector[i][1]);}
	else if ((moleculeBuildCommandVector[i])[0].compare("defineAndSetChargedAtomType") == 0 ) {

            cout<<__FILE__<<":"<<__LINE__<<" Syntax : defineAndSetChargedAtomType <biotype name> <FF atom class index> <charge> "<<endl;
            cout<<__FILE__<<":"<<__LINE__<<"     or : defineAndSetChargedAtomType <residue e.g. Phosphate,?RNA> <generic atom name e.g. OP> <ordinality : Initialy, Any, or Final> <FF atom class index> <charge> "<<endl;

            //cout<<__FILE__<<":"<<__LINE__<<" Contents of moleculeBuildCommandVector[i] : >"<<moleculeBuildCommandVector[i][0]<<"<, >"<<moleculeBuildCommandVector[i][1]<<"<, >"<<moleculeBuildCommandVector[i][2]<<"< ."<<endl;

            if ((moleculeBuildCommandVector[i]).size() == 4) {
                cout<<__FILE__<<":"<<__LINE__<<endl;
		DuMM::ChargedAtomTypeIndex 	myChargedAtomTypeIndex = dumm.getNextUnusedChargedAtomTypeIndex (); 
		String myBiotypeName = moleculeBuildCommandVector[i][1];
		dumm.defineChargedAtomType(myChargedAtomTypeIndex, 
		    myBiotypeName,           // biotype name. This is actually not used.
		    DuMM::AtomClassIndex(atoi(moleculeBuildCommandVector[i][2].c_str())), // force field atom class index
		    atof(moleculeBuildCommandVector[i][3].c_str()));	              // charge
		//if (
		dumm.setBiotypeChargedAtomType(myChargedAtomTypeIndex, compoundObjectMapContainer.fetchBiotype(myBiotypeName).getIndex());
		// For example:
		//                      index just needs to be unused.    doesn't matter.Class in force field.     charge
		//defineChargedAtomType(DuMM::ChargedAtomTypeIndex(5000), "Methane C",   DuMM::AtomClassIndex(1),  0.04);
		//setBiotypeChargedAtomType(DuMM::ChargedAtomTypeIndex(5000), Biotype::MethaneC().getIndex());

            } else if ((moleculeBuildCommandVector[i]).size() == 6) {
                cout<<__FILE__<<":"<<__LINE__<<endl;
		DuMM::ChargedAtomTypeIndex 	myChargedAtomTypeIndex = dumm.getNextUnusedChargedAtomTypeIndex (); 
		//#String myBiotypeName = moleculeBuildCommandVector[i][1];
		dumm.defineChargedAtomType(myChargedAtomTypeIndex, 
		    String("myBiotypeName"),           // biotype name. This is actually not used.
		    DuMM::AtomClassIndex(atoi(moleculeBuildCommandVector[i][4].c_str())), // force field atom class index
		    atof(moleculeBuildCommandVector[i][5].c_str()));	              // charge
		String residueName      = moleculeBuildCommandVector[i][1];
                for (int k = 0; k < residueName.length(); k++) {
                    if (String(residueName[k]).compare("?") ==0) {
                        residueName[k] = (String(" "))[0];
                    }
                }
                cout<<__FILE__<<":"<<__LINE__<<" Residue name is now : >"<<residueName<<"< "<<endl;
		String genericAtomName  = moleculeBuildCommandVector[i][2];
                String ordinalityString = moleculeBuildCommandVector[i][3];
                enum  	ordinalityEnum { Any = 1, Initial = 2, Final = 3 };
                Ordinality::Residue myOrdinality;
                if (ordinalityString.compare("Initial") == 0) myOrdinality =  SimTK::Ordinality::Initial;
                else if (ordinalityString.compare("Final") == 0) myOrdinality =  SimTK::Ordinality::Final;
                else if (ordinalityString.compare("Any") == 0) myOrdinality =  SimTK::Ordinality::Any;
                else {
		    ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" Invalid ordinality : >"<<ordinalityString<<"< "<<endl; 
		    ErrorManager::instance.treatError();
                }
                dumm.setBiotypeChargedAtomType (myChargedAtomTypeIndex, 
		    Biotype::get(residueName, genericAtomName, 
                        myOrdinality
                        ).getIndex()                
                    );

            } else {
	        ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" Wrong number of parameters ("<<(moleculeBuildCommandVector[i]).size()<<")!"<<endl; 
                ErrorManager::instance.treatError();
            }

        }
	else if ((moleculeBuildCommandVector[i])[0].compare("setDefaultBondAngle") == 0 ) {
            cout<<__FILE__<<":"<<__LINE__<<" Syntax : setDefaultBondAngle <Angle (in degrees)> <atom name 1> <atom name 2 (central atom)> <atom name 3>"<<endl;
            cout<<__FILE__<<":"<<__LINE__<<" example: setDefaultBondAngle 104.52 HW1 OW HW2 "<<endl; 
            double myAngle = atof(moleculeBuildCommandVector[i][1].c_str());
            String atomName1 = moleculeBuildCommandVector[i][2]; 
            String atomName2 = moleculeBuildCommandVector[i][3]; 
            String atomName3 = moleculeBuildCommandVector[i][4]; 
            setDefaultBondAngle(myAngle*Deg2Rad,atomName1,atomName2,atomName3);
        }
	else if ((moleculeBuildCommandVector[i])[0].compare("setBiotypeChargedAtomType") == 0 ) {
                cout<<__FILE__<<":"<<__LINE__<<" Syntax: setBiotypeChargedAtomType <residue name> <generic atom name> <ordinality> <FF atom class index> <charge>"<<endl;
                cout<<__FILE__<<":"<<__LINE__<<" Or    : setBiotypeChargedAtomType <(DuMM) chargedAtomTypeIndex> <biotypeIndex>"<<endl;
                /*if ((moleculeBuildCommandVector[i]).size() != 6) {
	            ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" Incorrect number of parameters!"<<endl; ErrorManager::instance.treatError();
                }*/
                DuMM::ChargedAtomTypeIndex      myChargedAtomTypeIndex;
                cout<<__FILE__<<":"<<__LINE__<<" myChargedAtomTypeIndex = >"<<myChargedAtomTypeIndex<<"< "<<endl;
                BiotypeIndex myBiotypeIndex;
                if ((moleculeBuildCommandVector[i]).size() == 6) {
			String residueName      = moleculeBuildCommandVector[i][1];
			for (int k = 0; k < residueName.length(); k++) {
			    if (String(residueName[k]).compare("?") ==0) {
				residueName[k] = (String(" "))[0];
			    } // of if 
			} // of for k
			cout<<__FILE__<<":"<<__LINE__<<" Residue name is now : >"<<residueName<<"< "<<endl;
			String genericAtomName  = moleculeBuildCommandVector[i][2];
			String ordinalityString = moleculeBuildCommandVector[i][3];
			DuMM::AtomClassIndex myAtomClassIndex ( atoi(moleculeBuildCommandVector[i][4].c_str())); // force field atom class index
			double charge = atof(moleculeBuildCommandVector[i][5].c_str());// charge
			enum  	ordinalityEnum { Any = 1, Initial = 2, Final = 3 };
			Ordinality::Residue myOrdinality;
			if (ordinalityString.compare("Initial") == 0) myOrdinality =  SimTK::Ordinality::Initial;
			else if (ordinalityString.compare("Final") == 0) myOrdinality =  SimTK::Ordinality::Final;
			else if (ordinalityString.compare("Any") == 0) myOrdinality =  SimTK::Ordinality::Any;
			else {
			    ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" Invalid ordinality : >"<<ordinalityString<<"< "<<endl; 
			    ErrorManager::instance.treatError();
			}
			myBiotypeIndex = Biotype::get(residueName, genericAtomName, 
				myOrdinality
				).getIndex();                
                        cout<<__FILE__<<":"<<__LINE__<<" myChargedAtomTypeIndex = >"<<myChargedAtomTypeIndex<<"< "<<endl;
			myChargedAtomTypeIndex = dumm.getNextUnusedChargedAtomTypeIndex (); 
                        cout<<__FILE__<<":"<<__LINE__<<" myChargedAtomTypeIndex = >"<<myChargedAtomTypeIndex<<"< "<<endl;
			String myBiotypeName = "blah";
		   
			dumm.defineChargedAtomType(myChargedAtomTypeIndex, 
				myBiotypeName,         // biotype name. This is actually not used.
				myAtomClassIndex,      // force field atom class index
			charge);	       // charge
                        dumm.setBiotypeChargedAtomType(myChargedAtomTypeIndex, myBiotypeIndex);
                } else if ((moleculeBuildCommandVector[i]).size() == 3) { cout<<__FILE__<<":"<<__LINE__<<endl;//exit(1);
                        if (dumm.hasChargedAtomType ((moleculeBuildCommandVector[i])[1])) {
                            cout<<__FILE__<<":"<<__LINE__<<" Found chargedAtomType for >"<<(moleculeBuildCommandVector[i])[1]<<"< "<<endl;
                        }
                        //myChargedAtomTypeIndex= atoi((moleculeBuildCommandVector[i])[1].c_str());
                        //myBiotypeIndex= atoi((moleculeBuildCommandVector[i])[2].c_str());
                        dumm.setBiotypeChargedAtomType(myChargedAtomTypeIndex, myBiotypeIndex);
                } else {    ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" Incorrect number of parameters!"<<endl; ErrorManager::instance.treatError();}
        }
        else {
            /*cout<<__FILE__<<":"<<__LINE__<<" Contents of moleculeBuildCommandVector[i] : "<<endl;;
            for (int j = 0; j < moleculeBuildCommandVector[i].size(); j++) {
                cout<<">"<<moleculeBuildCommandVector[i][j]<<"< ";
            }
            cout<<endl;*/
            ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" Unknown command!"<<endl; ErrorManager::instance.treatError();
        }
    }
};

/*void MoleculeClass::setChainID(String chainID) {
    molecule.setPdbChainId(chainID);
}*/

void  MoleculeClass::setPdbResidueName() {
    if (residueName.length() > 3) {
        ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" Error! residueName >"<<residueName <<"< is too long!" <<endl; ErrorManager::instance.treatError();
    }
    cout<<__FILE__<<":"<<__LINE__<<" setting residueName >"<<residueName <<"<"<<endl;
    molecule.setPdbResidueName(residueName);
};

void MoleculeClass::includeAllAtoms( DuMMForceFieldSubsystem & dumm) {
    for (Compound::AtomIndex i  = Compound::AtomIndex(0) ; i <molecule.getNumAtoms(); i++) {
	dumm.includeNonbondAtom(molecule.getDuMMAtomIndex(i));
    } 
}
void MoleculeClass::addRingClosingBond(CovalentBondClass myBond){
        cout<<__FILE__<<":"<<__LINE__<<" about start MoleculeClass::addRingClosingBond. molecule.getNumAtoms() = "<<molecule.getNumAtoms()<<endl;
        if (molecule.getNumAtoms() == 0) {
            ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" Error! Molecule has  "<<molecule.getNumAtoms() <<" atoms! " <<endl; 
            ErrorManager::instance.treatError();
        }
        const Compound::BondCenterPathName & centerName1 = /*String ("1/")+*/String(myBond.getAtomName1())+String('/')+String(myBond.getBondCenterName1());
        const Compound::BondCenterPathName & centerName2 = String(myBond.getAtomName2())+String('/')+String(myBond.getBondCenterName2());
        if (!(molecule.hasBondCenter(centerName1))){
            ErrorManager::instance<<__FILE__<<":"<<__LINE__<<" Unable to find bond center "<<centerName1<<std::endl;
            ErrorManager::instance.treatError();
    }
        if (!(molecule.hasBondCenter(centerName2))){
            ErrorManager::instance<<__FILE__<<":"<<__LINE__<<" Unable to find bond center "<<centerName2<<std::endl;
            ErrorManager::instance.treatError();
    }
    double bondLength = 111.111 ; // This doesn't matter, so I set to an absurd value
    double dihedralAngle = 0.0; // This doesn't matter either.
    
    molecule.addRingClosingBond( centerName1,    centerName2 , bondLength, dihedralAngle, SimTK::BondMobility::Rigid);
}
void MoleculeClassContainer::add(String myChainID,  String myResidueName, MoleculeClass & myMoleculeClass) {
    //myMoleculeClass.setChainID(chainID);
    myMoleculeClass.setChainID(myChainID);
    myMoleculeClass.setResidueName(myResidueName);
    if (moleculeClassMap.count(myChainID) > 0 ) {
        ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" Error! "<<moleculeClassMap.count(myChainID) <<" MoleculeClass's found with chain ID "<<myChainID <<endl; ErrorManager::instance.treatError();
    }
    moleculeClassMap.insert(std::pair <const String , MoleculeClass>  (myChainID,myMoleculeClass));
    //cout<<__FILE__<<":"<<__LINE__<<" Just inserted MoleculeClass with chain ID >"<<myChainID<<"< "<<endl;
    //cout<<__FILE__<<":"<<__LINE__<<" The inserted MoleculeClass has chain ID >"<<myMoleculeClass.getChainID()<<"< "<<endl;
    //cout<<__FILE__<<":"<<__LINE__<<" The retrieved MoleculeClass has chain ID >"<<updMoleculeClass(myChainID).getChainID()<<"< "<<endl;
    
}

void MoleculeClassContainer::adoptCompounds(SimTK::CompoundSystem & mySystem){
    map<const String, MoleculeClass>::iterator it;
    map<const String, MoleculeClass>::iterator next;
    //int i = 0;
    next = moleculeClassMap.begin();
    while (next != moleculeClassMap.end())
    {
       it = next;
       
       cout<<__FILE__<<":"<<__LINE__<<" About to adopt CustomMolecule with chain ID >"<<(it->second).getChainID()<<"<"<<endl;
       //cout<<__FILE__<<":"<<__LINE__<<" Top level transform before adopting: "<<(it->second.molecule).getTopLevelTransform()<<endl;
       mySystem.adoptCompound(it->second.molecule);
       //cout<<__FILE__<<":"<<__LINE__<<" Top level transform after adopting: "<<(it->second.molecule).getTopLevelTransform()<<endl;
       next++;
    }
};
void MoleculeClassContainer::initializeCompounds(DuMMForceFieldSubsystem & dumm){
    map<const String, MoleculeClass>::iterator it;
    map<const String, MoleculeClass>::iterator next;
    //int i = 0;
    next = moleculeClassMap.begin();
    while (next != moleculeClassMap.end())
    {
       it = next;
       (it->second).molecule = CustomMolecule((it->second).moleculeBuildCommandVector, dumm);
       //cout<<__FILE__<<":"<<__LINE__<<" About to set chain ID to >"<<(it->second).getChainID()<<"<"<<endl;
       (it->second).setPdbChainID(); // transfer MoleculeClass.chainID to the molecule member's PDB chain ID
       (it->second).setPdbResidueName(); // ditto for residue Name.
       (it->second).molecule.setPdbResidueNumber(1); // ditto for residue Name.
       next++;
    }
};

MoleculeClass &   MoleculeClassContainer::updMoleculeClass(String myChainID) {
    validateChainID(myChainID);
    map<const String, MoleculeClass>::iterator it;
    //if (moleculeClassMap.count(myChainID) == 1) {
    it = moleculeClassMap.find(myChainID);
    return it->second;
    //} else {
    //    ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" Error! "<<moleculeClassMap.count(myChainID) <<" MoleculeClass's found with chain ID "<<myChainID <<endl; ErrorManager::instance.treatError();
    //}
};

void  MoleculeClassContainer::validateChainID(String myChainID){
    if (moleculeClassMap.count(myChainID) != 1) {
        ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" Error! "<<moleculeClassMap.count(myChainID) <<" MoleculeClass's found with chain ID "<<myChainID <<endl; ErrorManager::instance.treatError();
    }
}

void MoleculeClassContainer::matchDefaultConfiguration(bool readPreviousFrameFile, String pdbFileName,bool matchExact, bool matchIdealized)
{
    cout<<__FILE__<<":"<<__LINE__<<" readPreviousFrameFile = "<<readPreviousFrameFile<<", pdbFileName = >"<<pdbFileName<<"< "<<endl;
    if (readPreviousFrameFile) {
        std::ifstream inputFile(pdbFileName.c_str(), ifstream::in);
        PdbStructure pdbStructure(inputFile);
	map<const String, MoleculeClass>::iterator it;
	map<const String, MoleculeClass>::iterator next;
	//int i = 0;
	next = moleculeClassMap.begin();
	while (next != moleculeClassMap.end())
	{
	   it = next;
	   Compound & myCompound = (it->second).molecule;
           cout <<__FILE__<<":"<<__LINE__<<" About to create atom targets from file "<<pdbFileName<<endl;
           cout <<__FILE__<<":"<<__LINE__<<myCompound.getPdbChainId()<<endl;
	   Compound::AtomTargetLocations atomTargets = myCompound.createAtomTargets(pdbStructure);
	   map<Compound::AtomIndex, Vec3>::iterator targetIt;
	   map<Compound::AtomIndex, Vec3>::iterator targetNext;
	   targetNext = atomTargets.begin();
	   while (targetNext != atomTargets.end())
	   { 
	      targetIt = targetNext;
              cout <<__FILE__<<":"<<__LINE__<<endl;//
              cout<<" "<<targetIt->second<<endl;
              targetNext++; 
           }
	   if (matchExact)
		    {
                    //cout<<__FILE__<<":"<<__LINE__<<" Top level transform before fitting: "<<myCompound.getTopLevelTransform()<<endl;
		    myCompound.matchDefaultConfiguration(atomTargets,   Compound::Match_Exact );
                    //cout<<__FILE__<<":"<<__LINE__<<" Top level transform after fitting: "<<myCompound.getTopLevelTransform()<<endl;
		    }
	   if (matchIdealized)
		    {myCompound.matchDefaultConfiguration(atomTargets,   Compound::Match_Idealized );} //planarity tolerance is in Radians, if Sherm's email is to be believed
	   next++;
        }
    }
}

void MoleculeClassContainer::includeAllAtoms( DuMMForceFieldSubsystem & dumm) {
	map<const String, MoleculeClass>::iterator it;
	map<const String, MoleculeClass>::iterator next;
	//int i = 0;
	next = moleculeClassMap.begin();
	while (next != moleculeClassMap.end())
	{
	   it = next;
           it->second.includeAllAtoms(dumm);
	   next++;
        }

}

bool MoleculeClassContainer::hasChainID( String chain) {
    cout<<__FILE__<<":"<<__LINE__<<" Inside MoleculeClassContainer::hasChainID("<<chain<<"). Found exactly "<<moleculeClassMap.count(chain)<<" molecules with chain ID >"<<chain<<"<."<<endl;
    if (moleculeClassMap.count(chain)== 1 ) {
        return true;}
    else if (moleculeClassMap.count(chain)== 0) {return false;} 
    else {
        ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" Error! "<<moleculeClassMap.count(chain) <<" MoleculeClass's found with chain ID "<<chain <<endl; ErrorManager::instance.treatError();

    }
}

void MoleculeClassContainer::addConstraintToGround(map<const String,double> myUserVariables,  const String chain, const String atomName, ConstraintToGroundContainer & constraintToGroundContainer){
    constraintToGroundContainer.addConstraintClassToVector(
        chain,
        ResidueID(1,(' ')), 
        atomName
        );
}



