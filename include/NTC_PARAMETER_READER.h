#ifndef NTC_PARAMETER_READER_H_
#define NTC_PARAMETER_READER_H_
/* -------------------------------------------------------------------------- *
 *                           MMB (MacroMoleculeBuilder)                       *
 * -------------------------------------------------------------------------- *
 *                                                                            *
 * Copyright (c) 2011-12 by the Author.                                       *
 * Author: Samuel Flores                                                      *
 *                                                                            *
 * See RNABuilder.cpp for the copyright and usage agreement.                  *
 * -------------------------------------------------------------------------- */


#include "SimTKmolmodel.h"
//#include "SimTKsimbody_aux.h"
#include "SimTKsimbody.h"
#include <ios>
#include <fstream>
#include <vector>
#include <iostream>
#include <string>
#include "Utils.h"
#include "NtC_Class_Container.h"

using namespace std;
using namespace SimTK;

/**
 * 
 * 
 * /param 
 * myPdbResidueName1,2 must be one of "A","C","G","U".
 * bondingEdge1,2 must be one of "WatsonCrick","Hoogsteen","Sugar","Bifurcated".
 * dihedraltype must be either "Cis" or "Trans".
 *
 */
struct NTC_PAR_BondRow   {
	String pdbResidueName1;
        String bondingEdge1; 
	String pdbResidueName2;
	String bondingEdge2;
	String dihedraltype;
        String residue1Atom[4];
        String residue2Atom[4];
        String atom_shift[4];
	double  bondLength[4];
	double  springConstant[4];
    double  CONFALVALUE;
	double torqueConstant;
	Vec3   attachmentPoint;
	double  rotationAngle;
	Vec3   rotationAxis;
        String isTwoTransformForce;
        double distanceC1pC1p;
        int    br;
};

class  NTC_PAR_BondKey {
    public:
	String pdbResidueName1;
	String pdbResidueName2;
        String bondingEdge1; 
	String bondingEdge2;
	String dihedraltype;
        String isTwoTransformForce;
        NTC_PAR_BondKey(String myPdbResidueName1, String myPdbResidueName2,String myBondingEdge1, String myBondingEdge2, String mydihedraltype, String myIsTwoTransformForce); 
        NTC_PAR_BondKey(NTC_PAR_BondRow myNTC_PAR_BondRow) ; 
};
struct NTC_PAR_BondKeyCmp {
    bool operator()( const NTC_PAR_BondKey ti1, const NTC_PAR_BondKey ti2 ) const {
        if (ti1.pdbResidueName1 < ti2.pdbResidueName1) return 1;
        else if (ti1.pdbResidueName1 > ti2.pdbResidueName1) return 0;
        else if ((ti1.pdbResidueName2 < ti2.pdbResidueName2)) return  1;
        else if ((ti1.pdbResidueName2 > ti2.pdbResidueName2)) return  0;
        else if ((ti1.bondingEdge1 < ti2.bondingEdge1)) return  1;
        else if ((ti1.bondingEdge1 > ti2.bondingEdge1)) return  0;
        else if ((ti1.bondingEdge2 < ti2.bondingEdge2)) return  1;
        else if ((ti1.bondingEdge2 > ti2.bondingEdge2)) return  0;
        else if ((ti1.dihedraltype < ti2.dihedraltype)) return  1;
        else if ((ti1.dihedraltype > ti2.dihedraltype)) return  0;
        else if ((ti1.isTwoTransformForce < ti2.isTwoTransformForce)) return  1;
        else if ((ti1.isTwoTransformForce > ti2.isTwoTransformForce)) return  0;
        else return 0;
    }
};

    static map <const NTC_PAR_BondKey, NTC_PAR_BondRow, NTC_PAR_BondKeyCmp> NTC_PAR_Map;


struct NTC_PAR_BondMatrix {
	vector<NTC_PAR_BondRow> myNTC_PAR_BondRow;
};
class NTC_PAR_Class  { 
public:
    NTC_PAR_BondMatrix myNTC_PAR_BondMatrix;
    int  initialize ( String inFileName) ;
    Transform getNTC_PAR_Transform(NTC_PAR_BondRow myNTC_PAR_BondRow) const;

    NTC_PAR_BondRow getNearestNTC_PAR_BondRow(String myPdbResidueName1,  String myPdbResidueName2, Transform residue1Transform, Transform residue2Transform)  const;
    void printNTC_PAR_BondRows ();

    int getNTC_PAR_BondRowIndex ( String myPdbResidueName1, String myPdbResidueName2, String Classtype, String dihedraltype,String myBasePairIsTwoTransformForce,NTC_Classes NTC) const;
    
    NTC_PAR_BondRow getNTC_PAR_BondRow(ResidueID myResidueNumber1,ResidueID myResidueNumber2, String myPdbResidueName1, String myBondingEdge1, String myPdbResidueName2,String myBondingEdge2, String mydihedraltype,String myBasePairIsTwoTransformForce) const  ;

};

#endif //      BaseInteractionParameterReader_H_
