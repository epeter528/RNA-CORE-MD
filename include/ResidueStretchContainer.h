#ifndef ResidueStretchContainer_H_
#define ResidueStretchContainer_H_
#include "BiopolymerClass.h"
#include "Utils.h"
#include "SimTKsimbody.h"
#include "SimTKmolmodel.h"
#include "RealVec.h"
#include "ReferenceNeighborList.h"

using namespace SimTK;

template <class ResidueStretchType>

class ResidueStretchContainer{
    //private:
    public:
    vector<ResidueStretchType> residueStretchVector;
    InterfaceContainer interfaceContainer;

    ResidueStretchContainer() {};

    void clear(){residueStretchVector.clear();interfaceContainer.clear(); };
    void validateResidueStretch(ResidueStretchType myResidueStretch, BiopolymerClassContainer & myBiopolymerClassContainer) {
	if (myBiopolymerClassContainer.updBiopolymerClass(myResidueStretch.getChain()).difference(myResidueStretch.getEndResidue() , myResidueStretch.getStartResidue()) < 0) {
	    ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" The end residue (currently "<<myResidueStretch.getEndResidue().outString()
		<<") must be greater than or equal to the start residue (currently "<<myResidueStretch.getStartResidue().outString()<<". "<<endl; ErrorManager::instance.treatError();
	}
	if ((myResidueStretch.getEndResidue() > myBiopolymerClassContainer.updBiopolymerClass(myResidueStretch.getChain()).getLastResidueID()) ) {
	    ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" The end residue (currently "<<myResidueStretch.getEndResidue().outString()
		<<") is greater than the last residue number of the chain."<<endl; ErrorManager::instance.treatError();
	}
	if ((myResidueStretch.getStartResidue() < myBiopolymerClassContainer.updBiopolymerClass(myResidueStretch.getChain()).getFirstResidueID()) ) {
	    ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" The start residue (currently "<<myResidueStretch.getStartResidue().outString()
		<<") is lesser than the first residue number of the chain."<<endl; ErrorManager::instance.treatError();
	}
	if (!(myBiopolymerClassContainer.hasChainID(myResidueStretch.getChain()))){
	    ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" Couldn't find chain "<<myResidueStretch.getChain()<<endl; ErrorManager::instance.treatError();
	}
    }

    void addStretch(ResidueStretchType newStretch) {
        residueStretchVector.push_back(newStretch);
    }    

    vector<ResidueStretchType> getResidueStretchVector() {return residueStretchVector;};

    void printResidueStretchVector() {
        for (int i = 0 ; i <residueStretchVector.size(); i++) {
            (residueStretchVector[i]).printStretch();
        }
    }
    /*void print(ResidueStretchType myResidueStretch)   {
        cout<<__FILE__<<":"<<__LINE__<<" Contents of residue stretch:  From residue : "<<myResidueStretch.getStartResidue().outString() <<" to residue: "<<myResidueStretch.getEndResidue().outString()<<", Chain : "<<myResidueStretch.getChain()<<endl;
    }*/
    void printResidueStretchVectorBreederFormat(std::ofstream & myofstream) {
        for (int i = 0 ; i <residueStretchVector.size(); i++) {
            printBreederFormat (residueStretchVector[i], myofstream);
        }
    }

    void printBreederFormat(ResidueStretchType myResidueStretch,std::ofstream & myofstream  )   { // print myResidueStretch in the format : C-NNNI , where C is the chain ID, NNN the residue number, and I the insertion code.
        myofstream<<myResidueStretch.getChain()<<"-"<<myResidueStretch.getStartResidue().outString()<<std::endl;
    }

    ResidueStretchType & getResidueStretch(int residueStretchIndex) {
        //print (residueStretchVector[residueStretchIndex]);
        return residueStretchVector[residueStretchIndex];
    };

    const int getNumResidueStretches() { 
        //cout<<__FILE__<<":"<<__LINE__<<" residueStretchVector.size() = "<<residueStretchVector.size()<<endl;
        return residueStretchVector.size();
    }

    vector <MobilizerStretch> getMobilizerStretchVector(String bondMobilityString)  {
		vector <MobilizerStretch> myMobilizerStretchVector;
		for (int i = 0; i < (int)residueStretchVector.size(); i++){
                        MobilizerStretch myMobilizerStretch ( residueStretchVector[i], bondMobilityString); 
			myMobilizerStretchVector.push_back(myMobilizerStretch);
              
                        cout<<__FILE__<<":"<<__LINE__<<" For chain "<<myMobilizerStretch.getChain()<< " start res "<<myMobilizerStretch.getStartResidue().outString()<<" end res "<< myMobilizerStretch.getEndResidue().outString() <<" mobility " << myMobilizerStretch.getBondMobility()<<endl;//" and chain "<< targetChain<<" residue "<<targetResidue.outString()<<" distance is "<<myDistance<<" nm"<<endl;
                        cout<< __FILE__<<":"<<__LINE__<<" Note that in prior releases of MMB we took dead lengths in Å.  For consistency with molmodel we are going back to nm, kJ/mol, ps, with apologies for the confusion."<<endl;

		}
		return myMobilizerStretchVector;
	};

    bool vectorHasResidueStretch(ResidueStretchType & residueStretch){
                //auto it = residueStretchVector.find(residueStretch); 
                // Can't just use find, perhaps because "<" and ">" operators are not defined.
                for (int i = 0 ; i <residueStretchVector.size(); i++) {
                    if (residueStretchVector[i] == residueStretch) return true;
                }
                return false;
                /*typename vector<ResidueStretchType>::iterator it;//residueStretchVectorIterator; 
                
		it = find(residueStretchVector.begin(), residueStretchVector.end(),residueStretch) ;
                if (it == residueStretchVector.end()) {return false;} else {return true;}*/
	};

    typename vector<ResidueStretchType>::iterator findResidueStretch(ResidueStretchType & residueStretch)
    {
        typename vector<ResidueStretchType>::iterator it;
        for(it = residueStretchVector.begin(); it != residueStretchVector.end(); it++)
        {
            if( *it == residueStretch)
                return it;
        }
        return residueStretchVector.end();
    }

	void addResidueStretchToVector(ResidueStretchType & residueStretch)
    {
        residueStretchVector.push_back(residueStretch);
    }

    void removeResidueStretchFromVector(ResidueStretchType & residueStretch)
    {
        typename vector<ResidueStretchType>::iterator it = findResidueStretch(residueStretch);
        if(it != residueStretchVector.end())
            residueStretchVector.erase(it);
        else
        {
            ErrorManager::instance <<__FILE__<<":"<<__LINE__<<": you tried to delete a non existing ResidueStretch: "<<"Stretch chain ="<< residueStretch.getChain()<<", first residue ="<<residueStretch.getStartResidue().outString()<<", last residue ="<<residueStretch.getEndResidue().outString()<< endl;
            ErrorManager::instance.treatError();
        }
    }

    void deleteResidueStretch(int id)
    {
        if(id >= 0 && id < residueStretchVector.size())
            residueStretchVector.erase(residueStretchVector.begin()+id);
        else
        {
            ErrorManager::instance <<__FILE__<<":"<<__LINE__<<": you tried to delete a non existing ResidueStretch." << endl;
            ErrorManager::instance.treatError();
        }
    }

    bool residuesAreWithinCutoff(double radius, BiopolymerClass & referenceBiopolymerClass, ResidueID referenceResidue, BiopolymerClass & targetBiopolymerClass, ResidueID targetResidue, BiopolymerClassContainer myBiopolymerClassContainer) {
                String referenceChain = referenceBiopolymerClass.getChainID();
                String targetChain = targetBiopolymerClass.getChainID();
		double myDistance = residueDistance(referenceBiopolymerClass, referenceResidue, targetBiopolymerClass, targetResidue, myBiopolymerClassContainer);
		if (myDistance <= radius) {return true;} else {return false;}
		
	};
				

    double  residueDistance( BiopolymerClass & referenceBiopolymerClass, ResidueID referenceResidue, BiopolymerClass & targetBiopolymerClass, ResidueID targetResidue, BiopolymerClassContainer myBiopolymerClassContainer) {
                //This is deprecated, not needed with new neighborlisting algorithm. Should extend for the case of mobilizerWithin, etc.
		//cout<<__FILE__<<":"<<__LINE__<<" Error!"<<endl; exit(1);
                String referenceChain = referenceBiopolymerClass.getChainID();
                String targetChain = targetBiopolymerClass.getChainID();
		double myDistance = (double)(
	                    	targetBiopolymerClass.calcDefaultAtomLocationInGroundFrame  ( targetResidue, targetBiopolymerClass.getRepresentativeAtomName())
                    		- referenceBiopolymerClass.calcDefaultAtomLocationInGroundFrame(referenceResidue, referenceBiopolymerClass.getRepresentativeAtomName())
                    	).norm(); // *10; // convert to Å -- no longer done, now using nm directly
                cout<<__FILE__<<":"<<__LINE__<<" For chain "<<referenceChain<<" residue "<<referenceResidue.outString()<<" and chain "<< targetChain<<" residue "<<targetResidue.outString()<<" distance is "<<myDistance<<" nm "<<endl;
                cout<< __FILE__<<":"<<__LINE__<<" Note that in prior releases of MMB we took distances in Å.  For consistency with molmodel we are going back to nm, kJ/mol, ps, with apologies for the confusion."<<endl;
                return myDistance;
	};
				




    void addTargetChainResidue(double radius, BiopolymerClass & referenceBiopolymerClass, ResidueID referenceResidue, BiopolymerClass & targetBiopolymerClass, ResidueID targetResidue,BiopolymerClassContainer myBiopolymerClassContainer) {
                String referenceChain = referenceBiopolymerClass.getChainID(); 
                String targetChain = targetBiopolymerClass.getChainID(); 
		if (referenceBiopolymerClass.getChainID().compare(targetBiopolymerClass.getChainID()) == 0) {
			ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" Chain "<<referenceChain<<" can't have an interface with itself!"<<endl;
			ErrorManager::instance.treatError();
		};
			residueStretchVector.push_back(ResidueStretchType(targetChain, targetResidue, targetResidue));
                
		cout<<__FILE__<<":"<<__LINE__<<" Chain "<<targetChain<<" residue "<<targetResidue.outString()<<" added to interface mobility zone. "<<endl;
	};


    void addMutualChainResidues(double radius, BiopolymerClass & referenceBiopolymerClass, BiopolymerClass & targetBiopolymerClass,BiopolymerClassContainer myBiopolymerClassContainer) {

		for (ResidueID referenceResidueID = referenceBiopolymerClass.getFirstResidueID(); referenceResidueID <= referenceBiopolymerClass.getLastResidueID(); referenceBiopolymerClass.incrementResidueID(referenceResidueID)) {
			for (ResidueID targetResidueID = targetBiopolymerClass.getFirstResidueID(); targetResidueID <= targetBiopolymerClass.getLastResidueID(); targetBiopolymerClass.incrementResidueID(targetResidueID))
			{
				addTargetChainResidue(radius, referenceBiopolymerClass, referenceResidueID, targetBiopolymerClass,        targetResidueID,       myBiopolymerClassContainer );	
                                addTargetChainResidue(radius, targetBiopolymerClass,        targetResidueID,       referenceBiopolymerClass, referenceResidueID, myBiopolymerClassContainer );
                                if (targetResidueID == targetBiopolymerClass.getLastResidueID()) {break;  }
			}
                        if (referenceResidueID == referenceBiopolymerClass.getLastResidueID()) {break;  }
                }
	};


        ////////////////////////////////////////////////
    #ifdef USE_OPENMM 
	// loop through all specified Interfaces, add interface residues to residueStretchVector
    void addStretchesToVectorFromInterfaceContainer(BiopolymerClassContainer  & myBiopolymerClassContainer) {
		for (int i = 0; i < interfaceContainer.numInterfaces(); i++) {
			Interface tempInterface = interfaceContainer.getInterface(i);
			addAllMutualChainResidues(tempInterface.Depth, tempInterface.Chains, tempInterface.PartnerChains, myBiopolymerClassContainer);
			cout<<__FILE__<<":"<<__LINE__<<endl;
		}
	};

        ////////////////////////////////////////////////
    #endif

    /*bool vectorCompare(String myString, vector<String> & comparisonStringVector) {
            if (comparisonStringVector.size() == 0) {return true;} // If we are comparing to an empty vector, return true.  This is in case no partner chains have been specified, in which case any chain will pass.
	    for (int i = 0; i < comparisonStringVector.size(); i++) {
                //cout<<__FILE__<<":"<<__LINE__<<" comparing "<<myString<< " to comparisonStringVector["<<i<<"] : "<<comparisonStringVector[i];
		if (comparisonStringVector[i].compare(myString ) == 0) { //cout<<", returning TRUE "<<endl ;  
                    return true; } else {//cout<<", returning FALSE";
                    }
                //cout<<endl;
	    }
	    return false; // If no String in comparisonStringVector is the same as myString
	}*/

        #ifdef USE_OPENMM
	void addAllMutualChainResidues(double radius, vector<String> referenceChains, vector<String> partnerChains,BiopolymerClassContainer & myBiopolymerClassContainer) { // This polymorphism requires that the user specify two sets of chains.  Only residues at the interface between the two sets will be included.  This lets the user leave out other chains (e.g. threading templates) which are in the system but which shouldn't be flexibilized.
            vector<MMBAtomInfo> concatenatedAtomInfoVector = myBiopolymerClassContainer.getConcatenatedAtomInfoVector();
	    vector<openmmVecType> particleList(concatenatedAtomInfoVector.size());
            for (int i = 0; i < concatenatedAtomInfoVector.size() ; i++) {
 		particleList[i] = concatenatedAtomInfoVector[i].position;
	    }
            cout<<__FILE__<<":"<<__LINE__<<endl;
	    vector<set<int> > exclusions( particleList.size() );
	    OpenMM::NeighborList neighborList;
	    //openmmVecType * boxSize ; 
            openmmVecType boxSize = openmmVecType(10000,10000,10000);
            cout<<__FILE__<<":"<<__LINE__<<" neighborList size is : "<<neighborList.size()<<endl;
	    computeNeighborListVoxelHash(neighborList, particleList.size() , particleList, exclusions, &boxSize, false, radius  , 0.0);
            //cout<<__FILE__<<":"<<__LINE__<<" neighborList size is : "<<neighborList.size()<<endl;
            //cout<<__FILE__<<":"<<__LINE__<<" neighborList size is : "<<neighborList.size()<<endl;
            for ( int j = 0 ; j < neighborList.size(); j++) {
                //cout<<__FILE__<<":"<<__LINE__<<" concatenatedAtomInfoVector[neighborList[j].first].chain ="<<concatenatedAtomInfoVector[neighborList[j].first].chain<< " ";
                //cout<<__FILE__<<":"<<__LINE__<<" concatenatedAtomInfoVector[neighborList[j].second].chain ="<<concatenatedAtomInfoVector[neighborList[j].second].chain<<endl;
                if (((( vectorCompare(concatenatedAtomInfoVector[neighborList[j].first].chain , (referenceChains))) == 1) &&
                     (( vectorCompare(concatenatedAtomInfoVector[neighborList[j].second].chain ,(  partnerChains))) == 1))  != //Use an XOR here. This means if the 'partnerChains' evaluation is later set to return 1 when partnerChains is empty, this will still work.
                    ((( vectorCompare(concatenatedAtomInfoVector[neighborList[j].second].chain ,(referenceChains))) == 1) &&
                     (( vectorCompare(concatenatedAtomInfoVector[neighborList[j].first].chain  ,(  partnerChains))) == 1))     //Make sure that exactly one residue is in the 'referenceChains', and the other residue is in the 'partnerChains' .. thus only the desired interface is included
															  )  
                 
		{
		    //cout<<__FILE__<<":"<<__LINE__<<endl;
                    ResidueStretchType myResidueStretch1;
                    myResidueStretch1.setStartResidue(concatenatedAtomInfoVector[neighborList[j].first].residueID);
                    myResidueStretch1.setEndResidue  (concatenatedAtomInfoVector[neighborList[j].first].residueID);
                    myResidueStretch1.setChain(concatenatedAtomInfoVector[neighborList[j].first].chain);
                    // It's not possible to set bondMobility here, also it's not necessary -- this is done in the calling function, MobilizerContainer::addMobilizerStretchesToVector
                    ResidueStretchType myResidueStretch2;
                    myResidueStretch2.setStartResidue(concatenatedAtomInfoVector[neighborList[j].second].residueID);
                    myResidueStretch2.setEndResidue(concatenatedAtomInfoVector[neighborList[j].second].residueID);
                    myResidueStretch2.setChain(concatenatedAtomInfoVector[neighborList[j].second].chain);
                    if (!(vectorHasResidueStretch(myResidueStretch1))) {addResidueStretchToVector(myResidueStretch1);
			cout<<__FILE__<<":"<<__LINE__<<" Added first residue stretch"<<endl; 
                        myResidueStretch1.printStretch();
                        //print (myResidueStretch1);
		    }
                    if (!(vectorHasResidueStretch(myResidueStretch2))) {addResidueStretchToVector(myResidueStretch2);
			cout<<__FILE__<<":"<<__LINE__<<" Added second residue stretch"<<endl; 
                        myResidueStretch2.printStretch() ;
                    }
	  	}
                else {
                    // Actually it seems unnecessary to actually delete the neighborList elements.  It should be enough not to call addResidueStretchToVector.
                }
            }
	};
        #endif
	
/*
	void addSingleWeldConstraintPerInterfaceChainPair( double radius, //vector<String> referenceChains, vector<String> partnerChains , ConstraintContainer & myConstraintContainer,  
            BiopolymerClassContainer & myBiopolymerClassContainer) { // This polymorphism requires that the user specify two sets of chains.  Only residues at the interface between the two sets will be included.  This lets the user leave out other chains (e.g. threading templates) which are in the system but which shouldn't be flexibilized.
            vector<MMBAtomInfo> concatenatedAtomInfoVector = myBiopolymerClassContainer.getConcatenatedAtomInfoVector();
	    vector<openmmVecType> particleList(concatenatedAtomInfoVector.size());
            for (int i = 0; i < concatenatedAtomInfoVector.size() ; i++) {
 		particleList[i] = concatenatedAtomInfoVector[i].position;
	    }
            cout<<__FILE__<<":"<<__LINE__<<endl;
	    vector<set<int> > exclusions( particleList.size() );
	    OpenMM::NeighborList neighborList;
            openmmVecType boxSize = openmmVecType(10000,10000,10000);
            cout<<__FILE__<<":"<<__LINE__<<" neighborList size is : "<<neighborList.size()<<endl;
	    computeNeighborListVoxelHash(neighborList, particleList.size() , particleList, exclusions, &boxSize, false, radius  , 0.0);
            //cout<<__FILE__<<":"<<__LINE__<<" neighborList size is : "<<neighborList.size()<<endl;
            //cout<<__FILE__<<":"<<__LINE__<<" neighborList size is : "<<neighborList.size()<<endl;
            for ( int j = 0 ; j < neighborList.size(); j++) {
                //cout<<__FILE__<<":"<<__LINE__<<" concatenatedAtomInfoVector[neighborList[j].first].chain ="<<concatenatedAtomInfoVector[neighborList[j].first].chain<< " ";
                //cout<<__FILE__<<":"<<__LINE__<<" concatenatedAtomInfoVector[neighborList[j].second].chain ="<<concatenatedAtomInfoVector[neighborList[j].second].chain<<endl;
                if (((( vectorCompare(concatenatedAtomInfoVector[neighborList[j].first].chain , (referenceChains))) == 1) &&
                     (( vectorCompare(concatenatedAtomInfoVector[neighborList[j].second].chain ,(  partnerChains))) == 1))  != //Use an XOR here. This means if the 'partnerChains' evaluation is later set to return 1 when partnerChains is empty, this will still work.
                    ((( vectorCompare(concatenatedAtomInfoVector[neighborList[j].second].chain ,(referenceChains))) == 1) &&
                     (( vectorCompare(concatenatedAtomInfoVector[neighborList[j].first].chain  ,(  partnerChains))) == 1))     //Make sure that exactly one residue is in the 'referenceChains', and the other residue is in the 'partnerChains' .. thus only the desired interface is included
															  )  
                 
		{
		    //cout<<__FILE__<<":"<<__LINE__<<endl;
                    ResidueID residueID1(concatenatedAtomInfoVector[neighborList[j].first].residueID);
                    String chain1(concatenatedAtomInfoVector[neighborList[j].first].chain);
                    String atom1(concatenatedAtomInfoVector[neighborList[j].first].atomName);
                    ResidueID residueID2(concatenatedAtomInfoVector[neighborList[j].second].residueID);
                    String chain2(concatenatedAtomInfoVector[neighborList[j].second].chain);
                    String atom2(concatenatedAtomInfoVector[neighborList[j].second].atomName);
	   	    if (!(myConstraintContainer.hasConstraintClass(myResidueStretch1.getChain(),myResidueStretch2.getChain()))) {
                        ConstraintClass myConstraintClass(chain1 ,residueID1,atom1,chain2, residueID2,atom2, ConstraintType::WeldToAtom);
		        myConstraintContainer.addConstraintClassToVector(myConstraintClass);
		        cout<<__FILE__<<":"<<__LINE__<<" Added ConstraintClass :"<<endl; 
		        myConstraintClass.print();
                        // Right here, should probably delete particleList  element.
		    }
	  	}
                else {
                    // Right here, should probably delete particleList element.
                }
            }
	};
*/
    #ifdef USE_OPENMM
    void addIntraChainInterfaceResidues(double radius, String referenceChain, BiopolymerClassContainer & myBiopolymerClassContainer) { // This polymorphism requires that the user specify ONE  chain.  Only residues at the interfaces between BODIES on that chain will be included.  This does NOT include interfaces with OTHER chains. 
            vector<MMBAtomInfo> chainAtomInfoVector = myBiopolymerClassContainer.updBiopolymerClass(referenceChain).getAtomInfoVector();
	    vector<openmmVecType> particleList(chainAtomInfoVector.size());
            for (int i = 0; i < chainAtomInfoVector.size() ; i++) {
                
 		particleList[i] = chainAtomInfoVector[i].position;
                //chainAtomInfoVector[i].print();
	    }
	    vector<set<int> > exclusions( particleList.size() );
            cout<<__FILE__<<":"<<__LINE__<<endl;
	    OpenMM::NeighborList neighborList;
	    ////openmmVecType * boxSize ;
            //*boxSize = openmmVecType(10000,10000,10000);
            openmmVecType boxSize = openmmVecType(10000,10000,10000);
            cout<<__FILE__<<":"<<__LINE__<<" neighborList size is : "<<neighborList.size()<<endl;
	    computeNeighborListVoxelHash(neighborList, particleList.size() , particleList, exclusions, &boxSize, false, radius  , 0.0);
            cout<<__FILE__<<":"<<__LINE__<<" in addIntraChainInterfaceResidues. About to add residues to physics zone." <<endl;
            cout<<__FILE__<<":"<<__LINE__<<" neighborList size is : "<<neighborList.size()<<endl;
            cout<<__FILE__<<":"<<__LINE__<<" depth = "<< radius <<endl;
            cout<<__FILE__<<":"<<__LINE__<<" chain = "<< referenceChain <<endl;
            for ( int j = 0 ; j < neighborList.size(); j++) {
                if (chainAtomInfoVector[neighborList[j].first].mobilizedBodyIndex < 0) {
		    ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" Bad chainAtomInfoVector[neighborList[j].first].mobilizedBodyIndex = " <<chainAtomInfoVector[neighborList[j].first].mobilizedBodyIndex <<endl; 
		    ErrorManager::instance.treatError(); }
                if (chainAtomInfoVector[neighborList[j].second].mobilizedBodyIndex < 0) {
		    ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" Bad chainAtomInfoVector[neighborList[j].second].mobilizedBodyIndex = " <<chainAtomInfoVector[neighborList[j].second].mobilizedBodyIndex <<endl; 
		    ErrorManager::instance.treatError(); }
                if (( chainAtomInfoVector[neighborList[j].first].mobilizedBodyIndex != chainAtomInfoVector[neighborList[j].second].mobilizedBodyIndex )) // &&  // if the two atoms belong to separate bodies, include them in the physics zone
                    //( chainAtomInfoVector[neighborList[j].first].chain.compare( referenceChain ) == 0) &&                                                   // the two atoms should also be in referenceChain
                    //( chainAtomInfoVector[neighborList[j].second].chain.compare( referenceChain ) == 0) ) 
		{ 
                    // To-Do: We need to start using a vector of atoms rather than residues.  
                    ResidueStretchType myResidueStretch1;
                    myResidueStretch1.setStartResidue(chainAtomInfoVector[neighborList[j].first].residueID);
                    //myResidueStretch1.setEndResidue  (chainAtomInfoVector[neighborList[j].first].residueID);
                    myResidueStretch1.setChain(chainAtomInfoVector[neighborList[j].first].chain);
                    ResidueStretchType myResidueStretch2;
                    myResidueStretch2.setStartResidue(chainAtomInfoVector[neighborList[j].second].residueID);
                    //myResidueStretch2.setEndResidue(chainAtomInfoVector[neighborList[j].second].residueID);
                    myResidueStretch2.setChain(chainAtomInfoVector[neighborList[j].second].chain);
                    if (!(vectorHasResidueStretch(myResidueStretch1))) {addResidueStretchToVector(myResidueStretch1);
			cout<<__FILE__<<":"<<__LINE__<<" Added first residue stretch"<<endl; 
                        myResidueStretch1.printStretch();
                        cout<<__FILE__<<":"<<__LINE__<<" based on chainAtomInfoVector["<<neighborList[j].first<<"].print() "<<endl; chainAtomInfoVector[neighborList[j].first].print();
                        cout<<__FILE__<<":"<<__LINE__<<" based on chainAtomInfoVector["<<neighborList[j].second<<"].print() "<<endl; chainAtomInfoVector[neighborList[j].second].print();
                        if (myResidueStretch1.getChain().compare(referenceChain) != 0) {
			    ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" Added undesired chain = "<<myResidueStretch1.getChain()<<endl;
			    ErrorManager::instance.treatError(); }
		    } // of if
                    if (!(vectorHasResidueStretch(myResidueStretch2))) {addResidueStretchToVector(myResidueStretch2);
			cout<<__FILE__<<":"<<__LINE__<<" Added second residue stretch"<<endl; 
                        cout<<__FILE__<<":"<<__LINE__<<" based on chainAtomInfoVector["<<neighborList[j].first<<"].print() "<<endl; chainAtomInfoVector[neighborList[j].first].print();
                        cout<<__FILE__<<":"<<__LINE__<<" based on chainAtomInfoVector["<<neighborList[j].second<<"].print() "<<endl; chainAtomInfoVector[neighborList[j].second].print();
                        myResidueStretch2.printStretch();
                        if (myResidueStretch2.getChain().compare(referenceChain) != 0) {
			    ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" Added undesired chain = "<<myResidueStretch2.getChain()<<endl;
			    ErrorManager::instance.treatError(); }
                    }
	  	}
                else {
                    //  It seems unnecessary to actually delete the neighborList elements.  It should be enough not to call addResidueStretchToVector.
                }
            }
            cout<<__FILE__<<":"<<__LINE__<<" done with addIntraChainInterfaceResidues. " <<endl;
	};
    #endif

};

typedef ResidueStretchContainer  <IncludeResidue> PhysicsContainer;

#endif

