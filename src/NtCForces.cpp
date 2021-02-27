/* -------------------------------------------------------------------------- *
 *                           MMB (MacroMoleculeBuilder)                       *
 * -------------------------------------------------------------------------- *
 *                                                                            *
 * Copyright (c) 2011-12 by the Author.                                       *
 * Author: Samuel Flores                                                      *
 *                                                                            *
 * See RNABuilder.cpp for the copyright and usage agreement.                  *
 * -------------------------------------------------------------------------- */
#include "NtCForces.h"  
#include <string.h>
#include <sstream>
#include <Utils.h>
#include "ParameterReader.h"

    NTC_Torque::NTC_Torque (SimbodyMatterSubsystem& matter,ParameterReader& myParameterReader,  NTC_PAR_Class& myNTC_PAR_Class, BiopolymerClassContainer & myBiopolymerClassContainer, std::ostream& outputStream ) : matter(matter),myParameterReader(myParameterReader), myNTC_PAR_Class (myNTC_PAR_Class), myBiopolymerClassContainer(myBiopolymerClassContainer), outputStream(outputStream)
        { 
    };    
    
    void NTC_Torque::calcForce(const State& state, Vector_<SpatialVec>& bodyForces,  
            Vector_<Vec3>& particleForces, Vector& mobilityForces) const 
        {  
        double energy = 0.0;        
        MobilizedBody body1;
        MobilizedBody body2;
        MobilizedBody body3;
        MobilizedBody body4;
        
        MobilizedBody body5;
        
        Transform transform1;
        Transform transform2;
        double forceConstant;
        double torqueConstant;
        double dutyCycle; //must be between 0 and 1.  at 1, force is applied all the time.  at 0, basically never applied.    
        double scrubberPeriod;
        double cutoffRadius;
        double pot_angle;
        double angle;
        double x_d1,x_d2,x_d3,x_d4;
        double y_d1,y_d2,y_d3,y_d4;
        double z_d1,z_d2,z_d3,z_d4;                   
        double d_d1_x,d_d1_y,d_d1_z;
        double d_d2_x,d_d2_y,d_d2_z;
        double d_d3_x,d_d3_y,d_d3_z;
        double cross_1_x,cross_1_y,cross_1_z;
        double cross_2_x,cross_2_y,cross_2_z;
        double d_t,d_t2;
        double cross_3_x,cross_3_y,cross_3_z;
        double direction_x,direction_y,direction_z;
        double rgsq,fg,hg,gaa,gbb,fga,hgb;
        double dfgx,dfgy,dfgz,dthx,dthy,dthz;
        double dtfx,dtfy,dtfz;
        double s_x2,s_y2,s_z2;
        double force_a1_x,force_a1_y,force_a1_z;
        double force_a2_x,force_a2_y,force_a2_z;        
        double force_a3_x,force_a3_y,force_a3_z;
        double force_a4_x,force_a4_y,force_a4_z;
        double PI = 3.14159265359;     
        double dih,bias;
        int    value,i;
        double prob[361];
        double diff_vec12 = 0.0;
        double diff_vec13 = 0.0;
        double diff_vec23 = 0.0;
        double diff_angle_123 = 0.0;       
        double fraction = 2.0;
        double fraction2 = 1.0;
        
        double  p_2body_x;
        double  p_2body_y;
        double  p_2body_z;    

        double  p_3body_x;
        double  p_3body_y;
        double  p_3body_z;         
        
        double alpha = 1E-6;
        double beta  = 1E-2;
        
        double alpha2,beta2,epsilon1,epsilon2;
        
        epsilon1 = 100.0;
        epsilon2 = 10.0;
        
        double ranA, ran2, ran, ran44;
        
        double b_norm,t_norm,dot_AB,cosgamma=0.0;
        
        double sam_x,sam_y,sam_z;
        
        AtomSpringContainer atomSpringContainer;
        
        ran44 = rand()%200/200.0;
        
   for (int r=0;r<myParameterReader.ntc_class_container.numNTC_Torsions();r++) 
        { 
        
            String chainId1=(myParameterReader.ntc_class_container.getNTC_Class(r)).NtC_FirstBPChain;            
            NTC_PAR_BondRow myNTC_PAR_BondRow = myNTC_PAR_Class.myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[(myParameterReader.ntc_class_container.getNTC_Class(r)).NTC_PAR_BondRowIndex];
            String basePairIsTwoTransformForce="ntcstep";
            ResidueID residueNumber1=(myParameterReader.ntc_class_container.getNTC_Class(r)).FirstBPResidue;
            ResidueID residueNumber2=(myParameterReader.ntc_class_container.getNTC_Class(r)).SecondBPResidue;
            ResidueID myResidueNumber,myResidueNumber2;

	    if ( myNTC_PAR_BondRow.bondLength[0] == 0.0 ) {
//		cout << "torsion " << r << " is real torsion" << endl;
            
            myBiopolymerClassContainer.cos_body_1[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]   = 0.0;
            myBiopolymerClassContainer.count_body_1[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 0.0;
            myBiopolymerClassContainer.cos_body_2[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]   = 0.0;
            myBiopolymerClassContainer.count_body_2[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 0.0;            
            myBiopolymerClassContainer.cos_body_3[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]   = 0.0;
            myBiopolymerClassContainer.count_body_3[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 0.0;         
         //   myResidueNumber = residueNumber1;            
         //   myResidueNumber.ResidueNumber += stoi(myNTC_PAR_BondRow.atom_shift[0]);             
             Vec3  state_1,state_2,state_3,state_4;
        
            if(stoi(myNTC_PAR_BondRow.atom_shift[0]) == 0) myResidueNumber = residueNumber1;
            if(stoi(myNTC_PAR_BondRow.atom_shift[0]) == 1) myResidueNumber = residueNumber2;
            
            body1 = myBiopolymerClassContainer.updAtomMobilizedBody(matter,chainId1,myResidueNumber,myNTC_PAR_BondRow.residue1Atom[0]);
            Vec3 stationA = myBiopolymerClassContainer.getAtomLocationInMobilizedBodyFrame(chainId1,myResidueNumber,myNTC_PAR_BondRow.residue1Atom[0]);
            state_1 = myBiopolymerClassContainer.calcAtomLocationInGroundFrame(state,chainId1,myResidueNumber,myNTC_PAR_BondRow.residue1Atom[0]); 
            
            if(stoi(myNTC_PAR_BondRow.atom_shift[1]) == 0) myResidueNumber = residueNumber1;
            if(stoi(myNTC_PAR_BondRow.atom_shift[1]) == 1) myResidueNumber = residueNumber2;
            
            body2 = myBiopolymerClassContainer.updAtomMobilizedBody(matter,chainId1,myResidueNumber,myNTC_PAR_BondRow.residue1Atom[1]);
            Vec3 stationB = myBiopolymerClassContainer.getAtomLocationInMobilizedBodyFrame(chainId1,myResidueNumber,myNTC_PAR_BondRow.residue1Atom[1]);
            state_2 = myBiopolymerClassContainer.calcAtomLocationInGroundFrame(state,chainId1,myResidueNumber,myNTC_PAR_BondRow.residue1Atom[1]); 
            
            if(stoi(myNTC_PAR_BondRow.atom_shift[2]) == 0) myResidueNumber = residueNumber1;
            if(stoi(myNTC_PAR_BondRow.atom_shift[2]) == 1) myResidueNumber = residueNumber2;
            
            body3 = myBiopolymerClassContainer.updAtomMobilizedBody(matter,chainId1,myResidueNumber,myNTC_PAR_BondRow.residue1Atom[2]);
            Vec3 stationC = myBiopolymerClassContainer.getAtomLocationInMobilizedBodyFrame(chainId1,myResidueNumber,myNTC_PAR_BondRow.residue1Atom[2]);
            state_3 = myBiopolymerClassContainer.calcAtomLocationInGroundFrame(state,chainId1,myResidueNumber,myNTC_PAR_BondRow.residue1Atom[2]); 

            if(stoi(myNTC_PAR_BondRow.atom_shift[3]) == 0) myResidueNumber = residueNumber1;
            if(stoi(myNTC_PAR_BondRow.atom_shift[3]) == 1) myResidueNumber = residueNumber2;

            body4 = myBiopolymerClassContainer.updAtomMobilizedBody(matter,chainId1,myResidueNumber,myNTC_PAR_BondRow.residue1Atom[3]);
            Vec3 stationD = myBiopolymerClassContainer.getAtomLocationInMobilizedBodyFrame(chainId1,myResidueNumber,myNTC_PAR_BondRow.residue1Atom[3]);
            state_4 = myBiopolymerClassContainer.calcAtomLocationInGroundFrame(state,chainId1,myResidueNumber,myNTC_PAR_BondRow.residue1Atom[3]);               
            
            torqueConstant = myNTC_PAR_BondRow.torqueConstant;            
            
            Vec3 d_d1,d_d2,d_d3;
            
            d_d1 = state_2 - state_1;
            d_d2 = state_3 - state_2;
            d_d3 = state_4 - state_3;
            
            Vec3 cross_1, cross_2;
            
            cross_1 = d_d1 % d_d2;
            cross_2 = d_d2 % d_d3;
            
            cross_1 = cross_1 / cross_1.norm();
            cross_2 = cross_2 / cross_2.norm();
        
            Vec3 cross_3;
            
            cross_3 = cross_1 % cross_2;
            
            
     angle = return_angle(cross_1,cross_2,cross_3,d_d2);
     
     double dist_ang = return_dist_ang(angle,myNTC_PAR_BondRow.rotationAngle);

    // pot_angle = torqueConstant*(dist_ang/57.295779513);//(exp(-(pow(dist_ang,2)/(2.0*pow(l_param,2)))))*(dist_ang/57.295779513)/(2.0*pow(l_param/57.295779513,2))*(sin(dist_ang/57.295779513));       
     
       if((myParameterReader.ntc_class_container.getNTC_Class(r)).meta == 0) pot_angle = torqueConstant*(myParameterReader.ntc_class_container.getNTC_Class(r)).weight*(-sin((dist_ang + 180.0)/57.295779513))*(360.0/57.295779513 + 1.0)/(1.0 + myNTC_PAR_BondRow.CONFALVALUE)/(360.0/57.295779513);
     
     Vec3 torque;    
     
       if((myParameterReader.ntc_class_container.getNTC_Class(r)).meta == 0) 
          
       { 
           
      //   cout << pot_angle << "  -- pot_angle " << endl;  
           
         torque = d_d2/d_d2.norm()*pot_angle;
   
         double virgin1,virgin2,virgin3,virgin4;
         
         virgin1 = sqrt(pow(bodyForces[body1.getMobilizedBodyIndex()][0][0],2) + pow(bodyForces[body1.getMobilizedBodyIndex()][0][1],2) + pow(bodyForces[body1.getMobilizedBodyIndex()][0][2],2));
         
         virgin2 = sqrt(pow(bodyForces[body2.getMobilizedBodyIndex()][0][0],2) + pow(bodyForces[body2.getMobilizedBodyIndex()][0][1],2) + pow(bodyForces[body2.getMobilizedBodyIndex()][0][2],2));

         virgin3 = sqrt(pow(bodyForces[body3.getMobilizedBodyIndex()][0][0],2) + pow(bodyForces[body3.getMobilizedBodyIndex()][0][1],2) + pow(bodyForces[body3.getMobilizedBodyIndex()][0][2],2));        
         
         virgin4 = sqrt(pow(bodyForces[body4.getMobilizedBodyIndex()][0][0],2) + pow(bodyForces[body4.getMobilizedBodyIndex()][0][1],2) + pow(bodyForces[body4.getMobilizedBodyIndex()][0][2],2));         
         
         double estranged;
         
         estranged = sqrt(pow(torque[0],2) + pow(torque[1],2) + pow(torque[2],2));
         
     if( (estranged) > 0.0 && virgin1 > 0.0) {     
         
         cosgamma = (bodyForces[body1.getMobilizedBodyIndex()][0][0]*(torque[0]) + 
                     bodyForces[body1.getMobilizedBodyIndex()][0][1]*(torque[1]) +
                     bodyForces[body1.getMobilizedBodyIndex()][0][2]*(torque[2]))/(virgin1*estranged);
         
     myBiopolymerClassContainer.l_max_tt[0] += cosgamma;
      myBiopolymerClassContainer.counter_l_max_tt[0] += 1.0;                     
                     
         torque[0] *= cosgamma;
         torque[1] *= cosgamma;         
         torque[2] *= cosgamma;
         
         bodyForces[body1.getMobilizedBodyIndex()] += SpatialVec(torque,Vec3(0));
         
       }
       
         estranged = sqrt(pow(torque[0],2) + pow(torque[1],2) + pow(torque[2],2));      
        
     if( (estranged) > 0.0 && virgin2 > 0.0) {         
         
         cosgamma = (bodyForces[body2.getMobilizedBodyIndex()][0][0]*(torque[0]) + 
                     bodyForces[body2.getMobilizedBodyIndex()][0][1]*(torque[1]) +
                     bodyForces[body2.getMobilizedBodyIndex()][0][2]*(torque[2]))/(virgin2*estranged);         
         
     myBiopolymerClassContainer.l_max_tt[0] += cosgamma;
      myBiopolymerClassContainer.counter_l_max_tt[0] += 1.0;                     
                     
         torque[0] *= cosgamma;
         torque[1] *= cosgamma;         
         torque[2] *= cosgamma;                     
                     
         bodyForces[body2.getMobilizedBodyIndex()] += SpatialVec(torque,Vec3(0));         

       }
         
         estranged = sqrt(pow(torque[0],2) + pow(torque[1],2) + pow(torque[2],2));    
        
     if( (estranged) > 0.0 && virgin3 > 0.0) {          
         
         cosgamma = (bodyForces[body3.getMobilizedBodyIndex()][0][0]*(torque[0]) + 
                     bodyForces[body3.getMobilizedBodyIndex()][0][1]*(torque[1]) +
                     bodyForces[body3.getMobilizedBodyIndex()][0][2]*(torque[2]))/(virgin3*estranged);         
        
     myBiopolymerClassContainer.l_max_tt[0] += cosgamma;
      myBiopolymerClassContainer.counter_l_max_tt[0] += 1.0;                     
                     
         torque[0] *= cosgamma;
         torque[1] *= cosgamma;         
         torque[2] *= cosgamma;                       
                     
         bodyForces[body3.getMobilizedBodyIndex()] -= SpatialVec(torque,Vec3(0));             
         
       }
         
         estranged = sqrt(pow(torque[0],2) + pow(torque[1],2) + pow(torque[2],2));          
        
     if( (estranged) > 0.0 && virgin4 > 0.0) {          
         
         cosgamma = (bodyForces[body4.getMobilizedBodyIndex()][0][0]*(torque[0]) + 
                     bodyForces[body4.getMobilizedBodyIndex()][0][1]*(torque[1]) +
                     bodyForces[body4.getMobilizedBodyIndex()][0][2]*(torque[2]))/(virgin4*estranged);          
         
     myBiopolymerClassContainer.l_max_tt[0] += cosgamma;
      myBiopolymerClassContainer.counter_l_max_tt[0] += 1.0;                     
                     
         torque[0] *= cosgamma;
         torque[1] *= cosgamma;         
         torque[2] *= cosgamma;                      
                     
         bodyForces[body4.getMobilizedBodyIndex()] -= SpatialVec(torque,Vec3(0));        
         
       }; 
      };
       
       if((myParameterReader.ntc_class_container.getNTC_Class(r)).meta == 1) {
           
        dih = 0;
        value = -1;
        bias = 0;
           
        angle *= 57.295779513;
         
       if(isfinite(angle) == 1) { 
        
        i = (int) round(angle);
        
        myBiopolymerClassContainer.hist[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][i] += 1.0;
        value = i;
                  
        myBiopolymerClassContainer.counter[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] += 1.0;
             
        if(myBiopolymerClassContainer.counter[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] > 0.0) myBiopolymerClassContainer.prob[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][i] = myBiopolymerClassContainer.hist[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][i]/myBiopolymerClassContainer.counter[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]; 
               
        if(value > 0.0 && myBiopolymerClassContainer.prob[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][value] > 1E-0 && myBiopolymerClassContainer.prob[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][value+1] > 1E-0 && value < 360) bias = 2.479*(log(myBiopolymerClassContainer.prob[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][value]/(1E-0))-log(myBiopolymerClassContainer.prob[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][value+1]/(1E-0)));
        if(prob[value] > 1E-0 && myBiopolymerClassContainer.prob[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][value+1] > 1E-0 && value == 360) bias = 2.479*(log(myBiopolymerClassContainer.prob[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][value]/(1E-0))-log(myBiopolymerClassContainer.prob[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][1]/(1E-0))); 
           
       /* if( isfinite(log(prob[value]/(1E-0))) == 1 && isfinite(log(prob[value+1]/(1E-0))) == 1) {   
           
         cout << bias << " value " << value << " bias "<< endl;
         cout << pot_angle << " pot angle " << endl;  
         cout << log(prob[value]/(1E-0)) << log(prob[value+1]/(1E-0)) << " log prob " << endl; 
         
        };*/
         
         pot_angle = torqueConstant*(-sin((dist_ang + 180.0)/57.295779513))*(360.0/57.295779513 + 1.0)/(1.0 + myNTC_PAR_BondRow.CONFALVALUE)/(360.0/57.295779513);
       
         torque = d_d2/d_d2.norm()*(pot_angle)/(1.0+(myParameterReader.ntc_class_container.getNTC_Class(r)).weight2)*(myParameterReader.ntc_class_container.getNTC_Class(r)).weight;
   
         bodyForces[body1.getMobilizedBodyIndex()] += SpatialVec(torque,Vec3(0));
         bodyForces[body4.getMobilizedBodyIndex()] -= SpatialVec(torque,Vec3(0));
         bodyForces[body2.getMobilizedBodyIndex()] += SpatialVec(torque,Vec3(0));
         bodyForces[body3.getMobilizedBodyIndex()] -= SpatialVec(torque,Vec3(0));          
         
        // cout << (pot_angle)/(1.0+(myParameterReader.ntc_class_container.getNTC_Class(r)).weight) << " pot angle " << endl;
        // cout << (myParameterReader.ntc_class_container.getNTC_Class(r)).weight << " weight " << endl;
        // cout << pot_angle << " pot angle " << endl;
         
       if( isfinite(log(myBiopolymerClassContainer.prob[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][value]/(1E-0))) == 1 && isfinite(log(myBiopolymerClassContainer.prob[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][value+1]/(1E-0))) == 1) {       
         
       if( isfinite(bias) == 1 && sqrt(pow(bias,2)) > 0.0) { 
           
        torque = -d_d2/d_d2.norm()*(bias)*sqrt(pow(pot_angle,2))/sqrt(pow(bias,2))*(myParameterReader.ntc_class_container.getNTC_Class(r)).weight2*(myParameterReader.ntc_class_container.getNTC_Class(r)).weight;
         
         bodyForces[body1.getMobilizedBodyIndex()] += SpatialVec(torque,Vec3(0));
         bodyForces[body4.getMobilizedBodyIndex()] -= SpatialVec(torque,Vec3(0));
         bodyForces[body2.getMobilizedBodyIndex()] += SpatialVec(torque,Vec3(0));
         bodyForces[body3.getMobilizedBodyIndex()] -= SpatialVec(torque,Vec3(0));  
  
        };
       };
     };
         
   }; 
       
// end real torsions
	}
// bonds
	else {

            Vec3  state_1;
            Vec3  state_2;
        
            body1 = myBiopolymerClassContainer.updAtomMobilizedBody(matter,chainId1,residueNumber1,myNTC_PAR_BondRow.residue1Atom[0]);
            Vec3 stationA = myBiopolymerClassContainer.getAtomLocationInMobilizedBodyFrame(chainId1,residueNumber1,myNTC_PAR_BondRow.residue1Atom[0]);
            state_1 = myBiopolymerClassContainer.calcAtomLocationInGroundFrame(state,chainId1,residueNumber1,myNTC_PAR_BondRow.residue1Atom[0]); 
            
            body2 = myBiopolymerClassContainer.updAtomMobilizedBody(matter,chainId1,residueNumber2,myNTC_PAR_BondRow.residue1Atom[1]);
            Vec3 stationB = myBiopolymerClassContainer.getAtomLocationInMobilizedBodyFrame(chainId1,residueNumber2,myNTC_PAR_BondRow.residue1Atom[1]);
            state_2 = myBiopolymerClassContainer.calcAtomLocationInGroundFrame(state,chainId1,residueNumber2,myNTC_PAR_BondRow.residue1Atom[1]); 
            
            Vec3 ptp = state_2 - state_1;
            double d = ptp.norm(); //sqrt(pow(x_d2-x_d1,2) + pow(y_d2-y_d1,2) + pow(z_d2-z_d1,2));
            double frc;
            Vec3   frcVec;
            
    if((myParameterReader.ntc_class_container.getNTC_Class(r)).meta == 0) {
      
    //  frc = (1.0-exp(-(2.0*myNTC_PAR_BondRow.CONFALVALUE)*(d-myNTC_PAR_BondRow.bondLength[0])))*(-exp(-(2.0*myNTC_PAR_BondRow.CONFALVALUE)*(d-myNTC_PAR_BondRow.bondLength[0])))*myNTC_PAR_BondRow.springConstant[0]*(myParameterReader.ntc_class_container.getNTC_Class(r)).weight;
        
        frc =  (1.0-1.0/(1.0 + exp(-(2.0*myNTC_PAR_BondRow.CONFALVALUE)*(d-myNTC_PAR_BondRow.bondLength[0]))))*(1.0/(1.0 + exp(-(2.0*myNTC_PAR_BondRow.CONFALVALUE)*(d-myNTC_PAR_BondRow.bondLength[0]))))*myNTC_PAR_BondRow.springConstant[0]*(myParameterReader.ntc_class_container.getNTC_Class(r)).weight;         
        
      frcVec = (frc)*ptp/d;
      
       bodyForces[body1.getMobilizedBodyIndex()] -=  SpatialVec(frcVec,Vec3(1));     
	   bodyForces[body2.getMobilizedBodyIndex()] +=  SpatialVec(frcVec,Vec3(1));       
      
    };
  
    
  if((myParameterReader.ntc_class_container.getNTC_Class(r)).meta == 1) {  
  
      double prob_d[31];
      
      double dist = 0.0;
      
      bias = 0.0;
      
       if(d < 3.0) {
           
       i = (int)round((d)*10.0);
       
       myBiopolymerClassContainer.hist_d[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][i] += 1.0;
       value = i;
                   
       myBiopolymerClassContainer.counter_d[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] += 1.0;
           
       if(myBiopolymerClassContainer.counter_d[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] > 0.0) myBiopolymerClassContainer.prob_d[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][i] = myBiopolymerClassContainer.hist_d[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][i]/myBiopolymerClassContainer.counter_d[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]; 
           
        if(value > 0.0 && myBiopolymerClassContainer.prob_d[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][value] > 1E-0 && myBiopolymerClassContainer.prob_d[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][value+1] > 1E-0 && value < 31) bias = 2.479*(log(myBiopolymerClassContainer.prob_d[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][value]/(1E-0))-log(myBiopolymerClassContainer.prob_d[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][value+1]/(1E-0)));

        frc =  (1.0-1.0/(1.0 + exp(-(2.0*myNTC_PAR_BondRow.CONFALVALUE)*(d-myNTC_PAR_BondRow.bondLength[0]))))*(1.0/(1.0 + exp(-(2.0*myNTC_PAR_BondRow.CONFALVALUE)*(d-myNTC_PAR_BondRow.bondLength[0]))))*myNTC_PAR_BondRow.springConstant[0]; 
        frc =  frc/(1.0 + (myParameterReader.ntc_class_container.getNTC_Class(r)).weight2);
        
        frcVec = (frc)*ptp/d*(myParameterReader.ntc_class_container.getNTC_Class(r)).weight;
        
	    bodyForces[body1.getMobilizedBodyIndex()] -=  SpatialVec(frcVec,Vec3(1));     
	    bodyForces[body2.getMobilizedBodyIndex()] +=  SpatialVec(frcVec,Vec3(1));        
       
     if( isfinite(bias) == 1 && sqrt(pow(bias,2)) > 0.0) {   
        
        bias = bias*(sqrt(pow(frc,2))/(sqrt(pow(bias,2))))*(myParameterReader.ntc_class_container.getNTC_Class(r)).weight2; 
        frcVec = (bias)*ptp/d*(myParameterReader.ntc_class_container.getNTC_Class(r)).weight;
        
	    bodyForces[body1.getMobilizedBodyIndex()] +=  SpatialVec(frcVec,Vec3(1));     
	    bodyForces[body2.getMobilizedBodyIndex()] -=  SpatialVec(frcVec,Vec3(1));          
        
      };
     };
     
   }; 
  };  
};


double prefactor_12=-1.0,prefactor_13=-1.0,prefactor_23=-1.0;
      
Vec3 d12;
Vec3 d13;
Vec3 d23;
double d_t_12;
double d_t_13;
double d_t_23;
int    pos_12_x, pos_12_y, pos_12_z;
int    pos_12b_x, pos_12b_y, pos_12b_z;
int    pos_13_x, pos_13_y, pos_13_z;
int    pos_23_x, pos_23_y, pos_23_z;
int    pos_angle;
Vec3 d21;
Vec3 frcVec;
double d_0;
double d_0_x,d_0_y,d_0_z;
double r_1;
double r_2;
double d_3;
int    k;

double dL_norm_2;
      
double hist_d_max;
double Z_state_X,Z_state_Y,Z_state_Z;
double Z_state_angle;

double sum_grad_angle=0.0;

double sum_grad_X=0.0, sum_grad_Y=0.0, sum_grad_Z=0.0;

Vec3   frcVec_B;

double delta[101];
double delta_x[101],delta_y[101],delta_z[101];

Vec3 state_1,state_2;

double d_x,d_x2,d_x3,d_y,d_y2,d_y3,d_z,d_z2,d_z3;
double dL_norm, dL_norm2, dL_norm3;
double delta_Lx_Z, delta_Ly_Z, delta_Lz_Z;

double d_tt,d_tt2;
double h3,h4;

double range_x,range_y,range_z;
double sigma_x, sigma_y, sigma_z;
int    pos_x,pos_y,pos_z;
double f_max_threebody_x,f_max_threebody_y,f_max_threebody_z;
double f_max_twobody_x,f_max_twobody_y,f_max_twobody_z;

int    all_refs_true;

double factorize = 1.0;//10.0*(double)myParameterReader.ntc_class_container.numNTC_Torsions();

int    r,n22;

for (int r=0;r<myParameterReader.ntc_class_container.numNTC_Torsions();r++) {  
    
  if((myParameterReader.ntc_class_container.getNTC_Class(r)).count != (myParameterReader.ntc_class_container.getNTC_Class(r+1)).count) {
     
     myBiopolymerClassContainer.counter_gen[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] += 1.0;  
      
    }  
    
   if( myBiopolymerClassContainer.counter_gen[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] == 1.0  || 
       myBiopolymerClassContainer.reset_3body[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] == 1) {
    
    myBiopolymerClassContainer.reset_3body[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 0;    
       
    myBiopolymerClassContainer.tau_2_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 1000.0;
    myBiopolymerClassContainer.tau_2_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 1000.0;
    myBiopolymerClassContainer.tau_2_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 1000.0; 
     
    myBiopolymerClassContainer.tau_1_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 1000.0;
    myBiopolymerClassContainer.tau_1_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 1000.0;
    myBiopolymerClassContainer.tau_1_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 1000.0; 
    
     myBiopolymerClassContainer.bias_angle[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = myBiopolymerClassContainer.dL_angle1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     /myBiopolymerClassContainer.counter_threebody_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];     
   
    myBiopolymerClassContainer.bias_3body_x_last[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] =
    (myBiopolymerClassContainer.dL12_x1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     /myBiopolymerClassContainer.counter_threebody_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     *myBiopolymerClassContainer.dL13_x1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     /myBiopolymerClassContainer.counter_threebody_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     *myBiopolymerClassContainer.dL23_x1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     /myBiopolymerClassContainer.counter_threebody_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count])
     *myBiopolymerClassContainer.bias_angle[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];
                
    myBiopolymerClassContainer.bias_3body_y_last[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] =
         (myBiopolymerClassContainer.dL12_y1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     /myBiopolymerClassContainer.counter_threebody_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     *myBiopolymerClassContainer.dL13_y1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     /myBiopolymerClassContainer.counter_threebody_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     *myBiopolymerClassContainer.dL23_y1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     /myBiopolymerClassContainer.counter_threebody_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count])
     *myBiopolymerClassContainer.bias_angle[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];

    myBiopolymerClassContainer.bias_3body_z_last[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] =
         (myBiopolymerClassContainer.dL12_z1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     /myBiopolymerClassContainer.counter_threebody_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     *myBiopolymerClassContainer.dL13_z1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     /myBiopolymerClassContainer.counter_threebody_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     *myBiopolymerClassContainer.dL23_z1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     /myBiopolymerClassContainer.counter_threebody_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count])
     *myBiopolymerClassContainer.bias_angle[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];               
                
    if(myBiopolymerClassContainer.counter_gen[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] == 1.0) {
               
       myBiopolymerClassContainer.bias_3body_x_last[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 1.0;
       myBiopolymerClassContainer.bias_3body_y_last[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 1.0;    
       myBiopolymerClassContainer.bias_3body_z_last[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 1.0;                  
                 
    };     
   
    myBiopolymerClassContainer.bias_3body_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 0.0;
    myBiopolymerClassContainer.bias_3body_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 0.0;    
    myBiopolymerClassContainer.bias_3body_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 0.0;    
    
    myBiopolymerClassContainer.bias_2body_norm[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 0.0;
    
    myBiopolymerClassContainer.bias_3body_norm[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 0.0;  
    
    myBiopolymerClassContainer.counter_threebody_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 0;
    myBiopolymerClassContainer.counter_threebody_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 0;    
    myBiopolymerClassContainer.counter_threebody_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 0;
    
    myBiopolymerClassContainer.av_angle[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 0.0;
    
    myBiopolymerClassContainer.dL_angle1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 0.0;
    
    myBiopolymerClassContainer.av_d12_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 0.0;
    myBiopolymerClassContainer.av_d12_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 0.0;    
    myBiopolymerClassContainer.av_d12_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 0.0;

    myBiopolymerClassContainer.dL12_x1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 0.0;
    myBiopolymerClassContainer.dL12_y1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 0.0;    
    myBiopolymerClassContainer.dL12_z1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 0.0;    
    
    myBiopolymerClassContainer.av_d13_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 0.0;
    myBiopolymerClassContainer.av_d13_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 0.0;    
    myBiopolymerClassContainer.av_d13_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 0.0;

    myBiopolymerClassContainer.dL13_x1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 0.0;
    myBiopolymerClassContainer.dL13_y1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 0.0;    
    myBiopolymerClassContainer.dL13_z1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 0.0;      

    myBiopolymerClassContainer.av_d23_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 0.0;
    myBiopolymerClassContainer.av_d23_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 0.0;    
    myBiopolymerClassContainer.av_d23_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 0.0;

    myBiopolymerClassContainer.dL23_x1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 0.0;
    myBiopolymerClassContainer.dL23_y1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 0.0;    
    myBiopolymerClassContainer.dL23_z1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 0.0;      
    
   for(k=1;k<=100;k++) {        
     
       myBiopolymerClassContainer.hist_global_12_3body_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][k] = 0.0;
       myBiopolymerClassContainer.hist_global_12_3body_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][k] = 0.0;       
       myBiopolymerClassContainer.hist_global_12_3body_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][k] = 0.0;
       
       myBiopolymerClassContainer.probability_x3_threebody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][k] = 0.0;
       myBiopolymerClassContainer.probability_y3_threebody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][k] = 0.0;       
       myBiopolymerClassContainer.probability_z3_threebody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][k] = 0.0;       
       
       myBiopolymerClassContainer.weight_12_3body_x[k] = 0.0;
       myBiopolymerClassContainer.weight_12_3body_y[k] = 0.0;       
       myBiopolymerClassContainer.weight_12_3body_z[k] = 0.0;
       
       myBiopolymerClassContainer.alpha_2_x_threebody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][k] = 0.0;
       myBiopolymerClassContainer.alpha_2_y_threebody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][k] = 0.0;                         
       myBiopolymerClassContainer.alpha_2_z_threebody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][k] = 0.0;
                         
       
      };     
    
    p_3body_x = 1.0;
    p_3body_y = 1.0;
    p_3body_z = 1.0;      
  
    p_2body_x = 1.0;
    p_2body_y = 1.0;
    p_2body_z = 1.0;     
    
   };      
      
  Vec3 state_3,stationA,stationB,stationC; 
   
    
  if((myParameterReader.ntc_class_container.getNTC_Class(r)).threebody == 1) {
      
 //  if((myParameterReader.ntc_class_container.getNTC_Class(r)).count != (myParameterReader.ntc_class_container.getNTC_Class(r+1)).count) {
       
    for(n22 = 0; n22 <= 0; n22 ++) {   
      
      String chainId1=(myParameterReader.ntc_class_container.getNTC_Class(r)).NtC_FirstBPChain; 
      ResidueID residueNumber1=(myParameterReader.ntc_class_container.getNTC_Class(r)).FirstBPResidue;      
      ResidueID residueNumber2=(myParameterReader.ntc_class_container.getNTC_Class(r)).sec_res_threebody;
      
      ResidueID residueNumber3=(myParameterReader.ntc_class_container.getNTC_Class(r)).third_res_threebody;
      
      NTC_PAR_BondRow myNTC_PAR_BondRow = myNTC_PAR_Class.myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[(myParameterReader.ntc_class_container.getNTC_Class(r)).NTC_PAR_BondRowIndex];

    all_refs_true = 0;  
      
    if( myBiopolymerClassContainer.updBiopolymerClass(chainId1).hasAtom(residueNumber1,"C2") == true) {      
      
      body1 = myBiopolymerClassContainer.updAtomMobilizedBody(matter,chainId1,residueNumber1,"C2");
      stationA = myBiopolymerClassContainer.getAtomLocationInMobilizedBodyFrame(chainId1,residueNumber1,"C2");
      state_1 = myBiopolymerClassContainer.calcAtomLocationInGroundFrame(state,chainId1,residueNumber1,"C2");   
      all_refs_true +=1;
      
    }
      
      String chainId2=(myParameterReader.ntc_class_container.getNTC_Class(r)).sec_res_chain_threebody;
     
    if( myBiopolymerClassContainer.updBiopolymerClass(chainId2).hasAtom(residueNumber2,"C2") == true) {      
      
      body2 = myBiopolymerClassContainer.updAtomMobilizedBody(matter,chainId2,residueNumber2,"C2");
      stationB = myBiopolymerClassContainer.getAtomLocationInMobilizedBodyFrame(chainId2,residueNumber2,"C2");
      state_2 = myBiopolymerClassContainer.calcAtomLocationInGroundFrame(state,chainId2,residueNumber2,"C2");     
      all_refs_true +=1;

    }
      
      String chainId3=(myParameterReader.ntc_class_container.getNTC_Class(r)).third_res_chain_threebody;
     
    if( myBiopolymerClassContainer.updBiopolymerClass(chainId3).hasAtom(residueNumber3,"C2") == true) {      
      
      body3 = myBiopolymerClassContainer.updAtomMobilizedBody(matter,chainId3,residueNumber3,"C2");
      stationC = myBiopolymerClassContainer.getAtomLocationInMobilizedBodyFrame(chainId3,residueNumber3,"C2");
      state_3 = myBiopolymerClassContainer.calcAtomLocationInGroundFrame(state,chainId3,residueNumber3,"C2"); 
      all_refs_true +=1;
      
    }
    
  if(all_refs_true == 3) {  
      
      Vec3 d12 = state_2 - state_1;
      Vec3 d13 = state_3 - state_1;
      Vec3 d23 = state_3 - state_2;
      
      d_t_12 = d12.norm();
      d_t_13 = d13.norm();
      d_t_23 = d23.norm();
      
  if((myParameterReader.ntc_class_container.getNTC_Class(r)).count != (myParameterReader.ntc_class_container.getNTC_Class(r+1)).count && n22 == 0){    
      
      diff_vec13 += d_t_13;
      diff_vec12 += d_t_12;
      diff_vec23 += d_t_23;     
      
  }
      
    d_3 = 0.0;      
      
    if(d12.norm() > 0.0 && d23.norm() > 0.0) 
    
    {    
        
    d_3 = dot(d12,d23)/(d12.norm()*d23.norm()); 
      
    angle = acos(d_3);

    angle *= 180.0/PI;
     
    };
    
   if((myParameterReader.ntc_class_container.getNTC_Class(r)).count != (myParameterReader.ntc_class_container.getNTC_Class(r+1)).count && n22 == 0) diff_angle_123 += angle;     
  
    if((myParameterReader.ntc_class_container.getNTC_Class(r)).count != (myParameterReader.ntc_class_container.getNTC_Class(r+1)).count && n22 == 0) { myBiopolymerClassContainer.counter_threebody_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] += 1.0;
    myBiopolymerClassContainer.counter_threebody_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] += 1.0;
    myBiopolymerClassContainer.counter_threebody_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] += 1.0;
    }
    
  //   std::cout << myBiopolymerClassContainer.counter_threebody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] << " counter 3 body " <<
  //   (myParameterReader.ntc_class_container.getNTC_Class(r)).count << "\n";
     
     myBiopolymerClassContainer.diff_angle[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]  = (angle);
     
     myBiopolymerClassContainer.diff_x12[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = (d12[0]);
     myBiopolymerClassContainer.diff_y12[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = (d12[1]);
     myBiopolymerClassContainer.diff_z12[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = (d12[2]);
     
     myBiopolymerClassContainer.diff_x13[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = (d13[0]); 
     myBiopolymerClassContainer.diff_y13[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = (d13[1]);
     myBiopolymerClassContainer.diff_z13[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = (d13[2]);      
     
     myBiopolymerClassContainer.diff_x23[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = (d23[0]);
     myBiopolymerClassContainer.diff_y23[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = (d23[1]);
     myBiopolymerClassContainer.diff_z23[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = (d23[2]);   

     myBiopolymerClassContainer.delta_x12[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = (d12[0] + d13[0] - d23[0])/3.0 - myBiopolymerClassContainer.d_12_last_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];
     myBiopolymerClassContainer.delta_y12[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = (d12[1] + d13[1] - d23[1])/3.0 - myBiopolymerClassContainer.d_12_last_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];     
     myBiopolymerClassContainer.delta_z12[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = (d12[2] + d13[2] - d23[2])/3.0 - myBiopolymerClassContainer.d_12_last_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];     

     //   myBiopolymerClassContainer.diff_x13[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = (d13[0]) - myBiopolymerClassContainer.d_13_last_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];
     //   myBiopolymerClassContainer.diff_y13[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = (d13[1]) - myBiopolymerClassContainer.d_13_last_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];     
     //   myBiopolymerClassContainer.diff_z13[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = (d13[2]) - myBiopolymerClassContainer.d_13_last_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];       

     //  myBiopolymerClassContainer.diff_x23[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = (d23[0]) - myBiopolymerClassContainer.d_23_last_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];
     //  myBiopolymerClassContainer.diff_y23[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = (d23[1]) - myBiopolymerClassContainer.d_23_last_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];     
     //  myBiopolymerClassContainer.diff_z23[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = (d23[2]) - myBiopolymerClassContainer.d_23_last_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];      

     myBiopolymerClassContainer.d_12_last_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = (d12[0] + d13[0] - d23[0])/3.0;
     myBiopolymerClassContainer.d_12_last_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = (d12[1] + d13[1] - d23[1])/3.0;     
     myBiopolymerClassContainer.d_12_last_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = (d12[2] + d13[2] - d23[2])/3.0;     

     if(isnan(myBiopolymerClassContainer.diff_angle[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] == 0))
	     myBiopolymerClassContainer.av_angle[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] += myBiopolymerClassContainer.diff_angle[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];

     if((int)myBiopolymerClassContainer.counter_gen[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]%((int)(myBiopolymerClassContainer.tau_1_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]) == 0)) { 

	     if(rand()%100/100.0 <= 1.0/(1.0 + 
				     exp(-myBiopolymerClassContainer.dL_angle1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]/
					     myBiopolymerClassContainer.counter_threebody_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count])) || //(int)(myBiopolymerClassContainer.counter_threebody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count])%100 == 0 ||
			     myBiopolymerClassContainer.counter_threebody_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] == 1) {    

		     myBiopolymerClassContainer.stamp_angle[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = myBiopolymerClassContainer.diff_angle[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];

	     }; 
     };

     double Corr_func_angle_t = 0.0;

     if((sqrt(pow(myBiopolymerClassContainer.stamp_angle[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] -
					     myBiopolymerClassContainer.av_angle[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
					     /myBiopolymerClassContainer.counter_threebody_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count],2))
			     *sqrt(pow(myBiopolymerClassContainer.diff_angle[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] - myBiopolymerClassContainer.av_angle[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]/myBiopolymerClassContainer.counter_threebody_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count],2))) > 0.0)

     {

	     myBiopolymerClassContainer.dL_angle1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] += (myBiopolymerClassContainer.stamp_angle[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] -
			     myBiopolymerClassContainer.av_angle[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
			     /myBiopolymerClassContainer.counter_threebody_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]    
			     )*(myBiopolymerClassContainer.diff_angle[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] - myBiopolymerClassContainer.av_angle[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
				     /myBiopolymerClassContainer.counter_threebody_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count])
			     /
			     (sqrt(pow(myBiopolymerClassContainer.stamp_angle[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] -
				       myBiopolymerClassContainer.av_angle[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
				       /myBiopolymerClassContainer.counter_threebody_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count],2))
			      *sqrt(pow(myBiopolymerClassContainer.diff_angle[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] - myBiopolymerClassContainer.av_angle[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]/myBiopolymerClassContainer.counter_threebody_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count],2)));      

     }

     // FINISH ANGLES

     // BEGIN D_12   
     if(isnan(myBiopolymerClassContainer.diff_x12[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]) == 0)
	     myBiopolymerClassContainer.av_d12_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] += myBiopolymerClassContainer.diff_x12[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];

     if(isnan(myBiopolymerClassContainer.diff_y12[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]) == 0)      
	     myBiopolymerClassContainer.av_d12_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] += myBiopolymerClassContainer.diff_y12[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];

     if(isnan(myBiopolymerClassContainer.diff_z12[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]) == 0)      
	     myBiopolymerClassContainer.av_d12_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] += myBiopolymerClassContainer.diff_z12[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];

     if((int)myBiopolymerClassContainer.counter_gen[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]%((int)(myBiopolymerClassContainer.tau_1_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]) == 0)) {   

	     if(rand()%100/100.0 <= 1.0/(1.0 + 
				     exp(-myBiopolymerClassContainer.dL12_x1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]/
					     myBiopolymerClassContainer.counter_threebody_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count])) ||
			     myBiopolymerClassContainer.counter_threebody_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] == 1) {         

		     myBiopolymerClassContainer.stamp_d12_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = myBiopolymerClassContainer.diff_x12[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];

	     };

	     if(rand()%100/100.0 <= 1.0/(1.0 + 
				     exp(-myBiopolymerClassContainer.dL12_y1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]/
					     myBiopolymerClassContainer.counter_threebody_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count])) ||
			     myBiopolymerClassContainer.counter_threebody_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] == 1) {        

		     myBiopolymerClassContainer.stamp_d12_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = myBiopolymerClassContainer.diff_y12[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];

	     };      

	     if(rand()%100/100.0 <= 1.0/(1.0 + 
				     exp(-myBiopolymerClassContainer.dL12_z1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]/
					     myBiopolymerClassContainer.counter_threebody_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count])) ||
			     myBiopolymerClassContainer.counter_threebody_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] == 1) {        

		     myBiopolymerClassContainer.stamp_d12_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = myBiopolymerClassContainer.diff_z12[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];        

	     };      

	     if(rand()%100/100.0 <= 1.0/(1.0 + 
				     exp(-myBiopolymerClassContainer.dL23_x1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]/
					     myBiopolymerClassContainer.counter_threebody_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count])) ||
			     myBiopolymerClassContainer.counter_threebody_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] == 1) {        

		     myBiopolymerClassContainer.stamp_d23_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = myBiopolymerClassContainer.diff_x23[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];

	     };      

	     if(rand()%100/100.0 <= 1.0/(1.0 + 
				     exp(-myBiopolymerClassContainer.dL23_y1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]/
					     myBiopolymerClassContainer.counter_threebody_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count])) ||
			     myBiopolymerClassContainer.counter_threebody_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] == 1) {        

		     myBiopolymerClassContainer.stamp_d23_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = myBiopolymerClassContainer.diff_y23[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];

	     };       

	     if(rand()%100/100.0 <= 1.0/(1.0 + 
				     exp(-myBiopolymerClassContainer.dL23_z1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]/
					     myBiopolymerClassContainer.counter_threebody_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count])) ||
			     myBiopolymerClassContainer.counter_threebody_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] == 1) {        

		     myBiopolymerClassContainer.stamp_d23_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = myBiopolymerClassContainer.diff_z23[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];        

	     };       

	     if(rand()%100/100.0 <= 1.0/(1.0 + 
				     exp(-myBiopolymerClassContainer.dL13_x1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]/
					     myBiopolymerClassContainer.counter_threebody_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count])) ||
			     myBiopolymerClassContainer.counter_threebody_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] == 1) {        

		     myBiopolymerClassContainer.stamp_d13_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = myBiopolymerClassContainer.diff_x13[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];

	     };         

	     if(rand()%100/100.0 <= 1.0/(1.0 + 
				     exp(-myBiopolymerClassContainer.dL13_y1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]/
					     myBiopolymerClassContainer.counter_threebody_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count])) ||
			     myBiopolymerClassContainer.counter_threebody_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] == 1) {        

		     myBiopolymerClassContainer.stamp_d13_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = myBiopolymerClassContainer.diff_y13[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];

	     };         

	     if(rand()%100/100.0 <= 1.0/(1.0 + 
				     exp(-myBiopolymerClassContainer.dL13_z1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]/
					     myBiopolymerClassContainer.counter_threebody_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count])) ||
			     myBiopolymerClassContainer.counter_threebody_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] == 1) {       

		     myBiopolymerClassContainer.stamp_d13_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = myBiopolymerClassContainer.diff_z13[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];        

	     };        
     };

     d_x  = 0.0;     
     d_x2 = 0.0;
     d_x3 = 0.0;
     d_y  = 0.0;     
     d_y2 = 0.0;
     d_y3 = 0.0;     
     d_z  = 0.0;     
     d_z2 = 0.0;
     d_z3 = 0.0;          
     
     double Corr_func_12_x=0.0;
     if((
			     sqrt(pow(myBiopolymerClassContainer.stamp_d12_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] -
					     myBiopolymerClassContainer.av_d12_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
					     /myBiopolymerClassContainer.counter_threebody_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count],2))    
			     *sqrt(pow(myBiopolymerClassContainer.diff_x12[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] -
					     myBiopolymerClassContainer.av_d12_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]/myBiopolymerClassContainer.counter_threebody_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count],2))) > 0.0)

     {

         d_x = (myBiopolymerClassContainer.stamp_d12_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] -
			     myBiopolymerClassContainer.av_d12_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
			     /myBiopolymerClassContainer.counter_threebody_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]    
			     )*(myBiopolymerClassContainer.diff_x12[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] - myBiopolymerClassContainer.av_d12_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
				     /myBiopolymerClassContainer.counter_threebody_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
			       )/(
				       sqrt(pow(myBiopolymerClassContainer.stamp_d12_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] -
						       myBiopolymerClassContainer.av_d12_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
						       /myBiopolymerClassContainer.counter_threebody_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count],2))    
				       *sqrt(pow(myBiopolymerClassContainer.diff_x12[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] -
						       myBiopolymerClassContainer.av_d12_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]/myBiopolymerClassContainer.counter_threebody_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count],2)));
         
	     myBiopolymerClassContainer.dL12_x1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] += d_x;

     };

     double Corr_func_12_y=0.0;
     if((
			     sqrt(pow(myBiopolymerClassContainer.stamp_d12_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] -
					     myBiopolymerClassContainer.av_d12_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
					     /myBiopolymerClassContainer.counter_threebody_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count],2))    
			     *sqrt(pow(myBiopolymerClassContainer.diff_y12[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] -
					     myBiopolymerClassContainer.av_d12_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]/myBiopolymerClassContainer.counter_threebody_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count],2))) > 0.0)

     {

         d_y = (myBiopolymerClassContainer.stamp_d12_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] -
			     myBiopolymerClassContainer.av_d12_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
			     /myBiopolymerClassContainer.counter_threebody_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]    
			     )*(myBiopolymerClassContainer.diff_y12[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] - myBiopolymerClassContainer.av_d12_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
				     /myBiopolymerClassContainer.counter_threebody_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
			       )/(
				       sqrt(pow(myBiopolymerClassContainer.stamp_d12_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] -
						       myBiopolymerClassContainer.av_d12_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
						       /myBiopolymerClassContainer.counter_threebody_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count],2))    
				       *sqrt(pow(myBiopolymerClassContainer.diff_y12[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] -
						       myBiopolymerClassContainer.av_d12_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]/myBiopolymerClassContainer.counter_threebody_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count],2)));
         
	     myBiopolymerClassContainer.dL12_y1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] += d_y;

     };

     double Corr_func_12_z=0.0;
     if((
			     sqrt(pow(myBiopolymerClassContainer.stamp_d12_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] -
					     myBiopolymerClassContainer.av_d12_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
					     /myBiopolymerClassContainer.counter_threebody_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count],2))    
			     *sqrt(pow(myBiopolymerClassContainer.diff_z12[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] -
					     myBiopolymerClassContainer.av_d12_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]/myBiopolymerClassContainer.counter_threebody_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count],2))) > 0.0)

     {

         d_z = (myBiopolymerClassContainer.stamp_d12_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] -
			     myBiopolymerClassContainer.av_d12_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
			     /myBiopolymerClassContainer.counter_threebody_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]    
			     )*(myBiopolymerClassContainer.diff_z12[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] - myBiopolymerClassContainer.av_d12_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
				     /myBiopolymerClassContainer.counter_threebody_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
			       )/(
				       sqrt(pow(myBiopolymerClassContainer.stamp_d12_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] -
						       myBiopolymerClassContainer.av_d12_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
						       /myBiopolymerClassContainer.counter_threebody_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count],2))    
				       *sqrt(pow(myBiopolymerClassContainer.diff_z12[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] -
						       myBiopolymerClassContainer.av_d12_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]/myBiopolymerClassContainer.counter_threebody_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count],2)));
         
	     myBiopolymerClassContainer.dL12_z1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] += d_z; 

     };

     // END 12
     if(isnan(myBiopolymerClassContainer.diff_x13[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]) == 0)      
	     myBiopolymerClassContainer.av_d13_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] += myBiopolymerClassContainer.diff_x13[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];

     if(isnan(myBiopolymerClassContainer.diff_y13[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]) == 0)      
	     myBiopolymerClassContainer.av_d13_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] += myBiopolymerClassContainer.diff_y13[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];

     if(isnan(myBiopolymerClassContainer.diff_z13[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]) == 0)      
	     myBiopolymerClassContainer.av_d13_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] += myBiopolymerClassContainer.diff_z13[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];

     double Corr_func_13_x=0.0;
     if((
			     sqrt(pow(myBiopolymerClassContainer.stamp_d13_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] -
					     myBiopolymerClassContainer.av_d13_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
					     /myBiopolymerClassContainer.counter_threebody_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count],2))    
			     *sqrt(pow(myBiopolymerClassContainer.diff_x13[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] -
					     myBiopolymerClassContainer.av_d13_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]/myBiopolymerClassContainer.counter_threebody_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count],2))) > 0.0) 

     {

         d_x2 = (myBiopolymerClassContainer.stamp_d13_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] -
			     myBiopolymerClassContainer.av_d13_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
			     /myBiopolymerClassContainer.counter_threebody_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]    
			     )*(myBiopolymerClassContainer.diff_x13[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] - myBiopolymerClassContainer.av_d13_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
				     /myBiopolymerClassContainer.counter_threebody_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
			       )/(
				       sqrt(pow(myBiopolymerClassContainer.stamp_d13_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] -
						       myBiopolymerClassContainer.av_d13_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
						       /myBiopolymerClassContainer.counter_threebody_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count],2))    
				       *sqrt(pow(myBiopolymerClassContainer.diff_x13[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] -
						       myBiopolymerClassContainer.av_d13_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]/myBiopolymerClassContainer.counter_threebody_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count],2)));
         
	     myBiopolymerClassContainer.dL13_x1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] += d_x2;

     };

     double Corr_func_13_y=0.0; 
     if((
			     sqrt(pow(myBiopolymerClassContainer.stamp_d13_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] -
					     myBiopolymerClassContainer.av_d13_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
					     /myBiopolymerClassContainer.counter_threebody_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count],2))    
			     *sqrt(pow(myBiopolymerClassContainer.diff_y13[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] -
					     myBiopolymerClassContainer.av_d13_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]/myBiopolymerClassContainer.counter_threebody_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count],2))) > 0.0)

     {

         d_y2 = (myBiopolymerClassContainer.stamp_d13_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] -
			     myBiopolymerClassContainer.av_d13_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
			     /myBiopolymerClassContainer.counter_threebody_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]    
			     )*(myBiopolymerClassContainer.diff_y13[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] - myBiopolymerClassContainer.av_d13_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
				     /myBiopolymerClassContainer.counter_threebody_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
			       )/(
				       sqrt(pow(myBiopolymerClassContainer.stamp_d13_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] -
						       myBiopolymerClassContainer.av_d13_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
						       /myBiopolymerClassContainer.counter_threebody_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count],2))    
				       *sqrt(pow(myBiopolymerClassContainer.diff_y13[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] -
						       myBiopolymerClassContainer.av_d13_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]/myBiopolymerClassContainer.counter_threebody_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count],2)));
         
	     myBiopolymerClassContainer.dL13_y1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] += d_y2;

     };

     double Corr_func_13_z=0.0; 
     if((
			     sqrt(pow(myBiopolymerClassContainer.stamp_d13_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] -
					     myBiopolymerClassContainer.av_d13_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
					     /myBiopolymerClassContainer.counter_threebody_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count],2))    
			     *sqrt(pow(myBiopolymerClassContainer.diff_z13[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] -
					     myBiopolymerClassContainer.av_d13_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]/myBiopolymerClassContainer.counter_threebody_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count],2))) > 0.0)

     {

         d_z2 = (myBiopolymerClassContainer.stamp_d13_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] -
			     myBiopolymerClassContainer.av_d13_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
			     /myBiopolymerClassContainer.counter_threebody_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]    
			     )*(myBiopolymerClassContainer.diff_z13[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] - myBiopolymerClassContainer.av_d13_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
				     /myBiopolymerClassContainer.counter_threebody_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
			       )/(
				       sqrt(pow(myBiopolymerClassContainer.stamp_d13_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] -
						       myBiopolymerClassContainer.av_d13_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
						       /myBiopolymerClassContainer.counter_threebody_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count],2))    
				       *sqrt(pow(myBiopolymerClassContainer.diff_z13[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] -
						       myBiopolymerClassContainer.av_d13_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]/myBiopolymerClassContainer.counter_threebody_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count],2))); 
         
	     myBiopolymerClassContainer.dL13_z1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] += d_z2;

     };

     // END 13

     // BEGIN 23

     if(isnan(myBiopolymerClassContainer.diff_x23[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]) == 0)     
	     myBiopolymerClassContainer.av_d23_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] += myBiopolymerClassContainer.diff_x23[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];

     if(isnan(myBiopolymerClassContainer.diff_y23[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]) == 0)      
	     myBiopolymerClassContainer.av_d23_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] += myBiopolymerClassContainer.diff_y23[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];

     if(isnan(myBiopolymerClassContainer.diff_z23[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]) == 0)      
	     myBiopolymerClassContainer.av_d23_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] += myBiopolymerClassContainer.diff_z23[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];

     //  if((int)(myBiopolymerClassContainer.counter_threebody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count])%100 == 0 ||
     //     myBiopolymerClassContainer.counter_threebody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] == 1) {         

     double Corr_func_23_x=0.0; 
     if((
			     sqrt(pow(myBiopolymerClassContainer.stamp_d23_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] -
					     myBiopolymerClassContainer.av_d23_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
					     /myBiopolymerClassContainer.counter_threebody_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count],2))    
			     *sqrt(pow(myBiopolymerClassContainer.diff_x23[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] -
					     myBiopolymerClassContainer.av_d23_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]/myBiopolymerClassContainer.counter_threebody_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count],2))) > 0.0) 

     {

         d_x3 = (myBiopolymerClassContainer.stamp_d23_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] -
			     myBiopolymerClassContainer.av_d23_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
			     /myBiopolymerClassContainer.counter_threebody_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]    
			     )*(myBiopolymerClassContainer.diff_x23[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] - myBiopolymerClassContainer.av_d23_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
				     /myBiopolymerClassContainer.counter_threebody_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
			       )/(
				       sqrt(pow(myBiopolymerClassContainer.stamp_d23_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] -
						       myBiopolymerClassContainer.av_d23_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
						       /myBiopolymerClassContainer.counter_threebody_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count],2))    
				       *sqrt(pow(myBiopolymerClassContainer.diff_x23[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] -
						       myBiopolymerClassContainer.av_d23_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]/myBiopolymerClassContainer.counter_threebody_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count],2)));
         
	     myBiopolymerClassContainer.dL23_x1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] += d_x3;

     };

     double Corr_func_23_y=0.0;
     if((
			     sqrt(pow(myBiopolymerClassContainer.stamp_d23_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] -
					     myBiopolymerClassContainer.av_d23_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
					     /myBiopolymerClassContainer.counter_threebody_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count],2))    
			     *sqrt(pow(myBiopolymerClassContainer.diff_y23[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] -
					     myBiopolymerClassContainer.av_d23_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]/myBiopolymerClassContainer.counter_threebody_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count],2))) > 0.0)    

     {

         d_y3 = (myBiopolymerClassContainer.stamp_d23_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] -
			     myBiopolymerClassContainer.av_d23_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
			     /myBiopolymerClassContainer.counter_threebody_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]    
			     )*(myBiopolymerClassContainer.diff_y23[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] - myBiopolymerClassContainer.av_d23_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
				     /myBiopolymerClassContainer.counter_threebody_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
			       )/(
				       sqrt(pow(myBiopolymerClassContainer.stamp_d23_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] -
						       myBiopolymerClassContainer.av_d23_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
						       /myBiopolymerClassContainer.counter_threebody_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count],2))    
				       *sqrt(pow(myBiopolymerClassContainer.diff_y23[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] -
						       myBiopolymerClassContainer.av_d23_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]/myBiopolymerClassContainer.counter_threebody_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count],2)));
         
	     myBiopolymerClassContainer.dL23_y1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] += d_y3;

     };

     double Corr_func_23_z=0.0; 
     if((
			     sqrt(pow(myBiopolymerClassContainer.stamp_d23_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] -
					     myBiopolymerClassContainer.av_d23_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
					     /myBiopolymerClassContainer.counter_threebody_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count],2))    
			     *sqrt(pow(myBiopolymerClassContainer.diff_z23[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] -
					     myBiopolymerClassContainer.av_d23_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]/myBiopolymerClassContainer.counter_threebody_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count],2))) > 0.0) 

     {

         d_z3 = (myBiopolymerClassContainer.stamp_d23_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] -
			     myBiopolymerClassContainer.av_d23_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
			     /myBiopolymerClassContainer.counter_threebody_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]    
			     )*(myBiopolymerClassContainer.diff_z23[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] - myBiopolymerClassContainer.av_d23_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
				     /myBiopolymerClassContainer.counter_threebody_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
			       )/(
				       sqrt(pow(myBiopolymerClassContainer.stamp_d23_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] -
						       myBiopolymerClassContainer.av_d23_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
						       /myBiopolymerClassContainer.counter_threebody_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count],2))    
				       *sqrt(pow(myBiopolymerClassContainer.diff_z23[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] -
						       myBiopolymerClassContainer.av_d23_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]/myBiopolymerClassContainer.counter_threebody_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count],2)));
         
	     myBiopolymerClassContainer.dL23_z1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] += d_z3;  

      }
     };
    };
    };     
  };

  int r22;

  r22 = 0;

  for (int r=0;r<myParameterReader.ntc_class_container.numNTC_Torsions();r++) { 

	  if((myParameterReader.ntc_class_container.getNTC_Class(r)).twobody == 1) {

		  // if((myParameterReader.ntc_class_container.getNTC_Class(r)).count != (myParameterReader.ntc_class_container.getNTC_Class(r+1)).count) {    

		  for(n22=0; n22 <= 0; n22 ++) {     

			  String chainId1=(myParameterReader.ntc_class_container.getNTC_Class(r)).NtC_FirstBPChain; 

			  ResidueID residueNumber1=(myParameterReader.ntc_class_container.getNTC_Class(r)).FirstBPResidue;      
			  ResidueID residueNumber2=(myParameterReader.ntc_class_container.getNTC_Class(r)).sec_res_twobody;

			  NTC_PAR_BondRow myNTC_PAR_BondRow = myNTC_PAR_Class.myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[(myParameterReader.ntc_class_container.getNTC_Class(r)).NTC_PAR_BondRowIndex];

              all_refs_true = 0;  
              
			  if( myBiopolymerClassContainer.updBiopolymerClass(chainId1).hasAtom(residueNumber1,"C2") == true) {      

				  body1 = myBiopolymerClassContainer.updAtomMobilizedBody(matter,chainId1,residueNumber1,"C2");
				  Vec3 stationA = myBiopolymerClassContainer.getAtomLocationInMobilizedBodyFrame(chainId1,residueNumber1,"C2");
				  state_1 = myBiopolymerClassContainer.calcAtomLocationInGroundFrame(state,chainId1,residueNumber1,"C2");       
                  all_refs_true += 1;
                  
			  }

			  String chainId2=(myParameterReader.ntc_class_container.getNTC_Class(r)).sec_res_chain_twobody;

			  if( myBiopolymerClassContainer.updBiopolymerClass(chainId2).hasAtom(residueNumber2,"C2") == true) {      

				  body2 = myBiopolymerClassContainer.updAtomMobilizedBody(matter,chainId2,residueNumber2,"C2");
				  Vec3 stationB = myBiopolymerClassContainer.getAtomLocationInMobilizedBodyFrame(chainId2,residueNumber2,"C2");
				  state_2 = myBiopolymerClassContainer.calcAtomLocationInGroundFrame(state,chainId2,residueNumber2,"C2");     \
                  all_refs_true += 1;
				  
			  }

		if(all_refs_true == 2) {	  
			  

			  Vec3 d12 = state_2 - state_1;

			  d_t_12 = d12.norm();

			  if((myParameterReader.ntc_class_container.getNTC_Class(r)).count != (myParameterReader.ntc_class_container.getNTC_Class(r+1)).count && n22 == 0){   

				  diff_vec12 += d_t_12;

			  }

			  //  d12 = (state_1%state_2);

			  if((myParameterReader.ntc_class_container.getNTC_Class(r)).count != (myParameterReader.ntc_class_container.getNTC_Class(r+1)).count && n22 == 0)
				  myBiopolymerClassContainer.counter_twobody_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] += 1.0;
			  if((myParameterReader.ntc_class_container.getNTC_Class(r)).count != (myParameterReader.ntc_class_container.getNTC_Class(r+1)).count && n22 == 0)
				  myBiopolymerClassContainer.counter_twobody_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] += 1.0;      
			  if((myParameterReader.ntc_class_container.getNTC_Class(r)).count != (myParameterReader.ntc_class_container.getNTC_Class(r+1)).count && n22 == 0)
				  myBiopolymerClassContainer.counter_twobody_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] += 1.0;

			  //    std::cout << myBiopolymerClassContainer.counter_twobody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] << " " << (myParameterReader.ntc_class_container.getNTC_Class(r)).count << " r: " << r << " -- counter twobody \n";

			  if(myBiopolymerClassContainer.counter_gen[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] == 1.0 && r22==0) {

				  myBiopolymerClassContainer.av_d12b_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 0.0;

			  }

			  if(myBiopolymerClassContainer.counter_gen[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] == 1.0 && r22==0) {   

				  myBiopolymerClassContainer.av_d12b_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 0.0;

			  }

			  if(myBiopolymerClassContainer.counter_gen[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] == 1.0 && r22==0) {   

				  myBiopolymerClassContainer.av_d12b_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 0.0;

			  }

			  r22 ++;
              
              if(myBiopolymerClassContainer.counter_gen[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] == 1.0 || 
                 myBiopolymerClassContainer.reset_2body[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] == 1) {
                 
                myBiopolymerClassContainer.reset_2body[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 0;  
                  
                for(k=1;i<=100;i++) {
                
                    myBiopolymerClassContainer.alpha_2_x_twobody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][k] = 0.0;
                    myBiopolymerClassContainer.alpha_2_y_twobody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][k] = 0.0;                         
                    myBiopolymerClassContainer.alpha_2_z_twobody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][k] = 0.0;
                         
                };
                   
                myBiopolymerClassContainer.bias_2body_x_last[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] =
                myBiopolymerClassContainer.dL12b_x1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
                /myBiopolymerClassContainer.counter_twobody_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];
                
                myBiopolymerClassContainer.bias_2body_y_last[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] =
                myBiopolymerClassContainer.dL12b_y1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
                /myBiopolymerClassContainer.counter_twobody_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];

                myBiopolymerClassContainer.bias_2body_z_last[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] =
                myBiopolymerClassContainer.dL12b_z1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
                /myBiopolymerClassContainer.counter_twobody_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];               
              
                myBiopolymerClassContainer.dL12b_x1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 0.0;
                myBiopolymerClassContainer.dL12b_y1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 0.0;                
                myBiopolymerClassContainer.dL12b_z1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 0.0;
                
                myBiopolymerClassContainer.av_d12b_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 0.0;
                myBiopolymerClassContainer.av_d12b_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 0.0;                
                myBiopolymerClassContainer.av_d12b_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 0.0;

                myBiopolymerClassContainer.counter_twobody_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 1.0;
                myBiopolymerClassContainer.counter_twobody_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 1.0;                
                myBiopolymerClassContainer.counter_twobody_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 1.0;
                
             if(myBiopolymerClassContainer.counter_gen[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] == 1.0) {
               
                myBiopolymerClassContainer.bias_2body_x_last[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 1.0;
                myBiopolymerClassContainer.bias_2body_y_last[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 1.0;    
                myBiopolymerClassContainer.bias_2body_z_last[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 1.0;                  
                 
             };   
                
                myBiopolymerClassContainer.bias_2body_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 0.0;
                myBiopolymerClassContainer.bias_2body_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 0.0;    
                myBiopolymerClassContainer.bias_2body_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 0.0;                
                
          for(k=1;k<=100;k++) {
         
                myBiopolymerClassContainer.hist_global_12_2body_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][k] = 0.0;
                myBiopolymerClassContainer.hist_global_12_2body_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][k] = 0.0;
                myBiopolymerClassContainer.hist_global_12_2body_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][k] = 0.0;
          
                myBiopolymerClassContainer.probability_x3_twobody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][k] = 0.0;
                myBiopolymerClassContainer.probability_y3_twobody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][k] = 0.0;       
                myBiopolymerClassContainer.probability_z3_twobody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][k] = 0.0;           
          
                myBiopolymerClassContainer.weight_12_2body_x[k] = 0.0;
                myBiopolymerClassContainer.weight_12_2body_y[k] = 0.0;          
                myBiopolymerClassContainer.weight_12_2body_z[k] = 0.0;
                 
          
                };               
              };

			  myBiopolymerClassContainer.diff_x12b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = (d12[0]);
			  myBiopolymerClassContainer.diff_y12b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = (d12[1]); 
			  myBiopolymerClassContainer.diff_z12b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = (d12[2]); 

			  myBiopolymerClassContainer.delta_x12b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = (d12[0] -
					  myBiopolymerClassContainer.d_12b_last_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]);
			  myBiopolymerClassContainer.delta_y12b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = (d12[1] -
					  myBiopolymerClassContainer.d_12b_last_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]); 
			  myBiopolymerClassContainer.delta_z12b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = (d12[2] -
					  myBiopolymerClassContainer.d_12b_last_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]);       

			  myBiopolymerClassContainer.d_12b_last_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = (d12[0]);
			  myBiopolymerClassContainer.d_12b_last_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = (d12[1]);     
			  myBiopolymerClassContainer.d_12b_last_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = (d12[2]);        

			  if(isnan(myBiopolymerClassContainer.diff_x12b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]) == 0)  myBiopolymerClassContainer.av_d12b_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] += myBiopolymerClassContainer.diff_x12b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];

			  if(isnan(myBiopolymerClassContainer.diff_y12b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]) == 0)      
				  myBiopolymerClassContainer.av_d12b_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] += myBiopolymerClassContainer.diff_y12b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];

			  if(isnan(myBiopolymerClassContainer.diff_z12b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]) == 0)      
				  myBiopolymerClassContainer.av_d12b_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] += myBiopolymerClassContainer.diff_z12b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];

			  if((int)myBiopolymerClassContainer.counter_gen[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]%((int)myBiopolymerClassContainer.tau_2_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]) == 0) {   

				  if(rand()%100/100.0 <= 1.0/(1.0 + 
							  exp(-myBiopolymerClassContainer.dL12b_x1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]/
								  myBiopolymerClassContainer.counter_twobody_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count])) ||
						  myBiopolymerClassContainer.counter_twobody_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] == 1) {        

					  myBiopolymerClassContainer.stamp_d12b_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = myBiopolymerClassContainer.diff_x12b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];

				  };
			  }

			  if((int)myBiopolymerClassContainer.counter_gen[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]%((int)myBiopolymerClassContainer.tau_2_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]) == 0) {        

				  if(rand()%100/100.0 <= 1.0/(1.0 + 
							  exp(-myBiopolymerClassContainer.dL12b_y1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]/
								  myBiopolymerClassContainer.counter_twobody_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count])) ||
						  myBiopolymerClassContainer.counter_twobody_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] == 1) {       

					  myBiopolymerClassContainer.stamp_d12b_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = myBiopolymerClassContainer.diff_y12b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];

				  };  
			  }

			  if((int)myBiopolymerClassContainer.counter_gen[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]%((int)myBiopolymerClassContainer.tau_2_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]) == 0) {      

				  if(rand()%100/100.0 <= 1.0/(1.0 + 
							  exp(-myBiopolymerClassContainer.dL12b_z1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]/
								  myBiopolymerClassContainer.counter_twobody_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count])) ||
						  myBiopolymerClassContainer.counter_twobody_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] == 1) {         

					  myBiopolymerClassContainer.stamp_d12b_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = myBiopolymerClassContainer.diff_z12b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];        

				  };    
			  };

			  double Corr_func_12b_x = 0.0;

              d_x = 0.0;
              d_y = 0.0;
              d_z = 0.0;
              
			  if((
						  sqrt(pow(myBiopolymerClassContainer.stamp_d12b_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] -
								  myBiopolymerClassContainer.av_d12b_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
								  /myBiopolymerClassContainer.counter_twobody_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count],2))    
						  *sqrt(pow(myBiopolymerClassContainer.diff_x12b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] -
								  myBiopolymerClassContainer.av_d12b_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]/myBiopolymerClassContainer.counter_twobody_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count],2))) > 0.0) 

			  {

                  d_x = (myBiopolymerClassContainer.stamp_d12b_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] -
							  myBiopolymerClassContainer.av_d12b_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
							  /myBiopolymerClassContainer.counter_twobody_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]    
					     )*(myBiopolymerClassContainer.diff_x12b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] - myBiopolymerClassContainer.av_d12b_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
						     /myBiopolymerClassContainer.counter_twobody_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
					       )/(
						       sqrt(pow(myBiopolymerClassContainer.stamp_d12b_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] -
								       myBiopolymerClassContainer.av_d12b_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
								       /myBiopolymerClassContainer.counter_twobody_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count],2))    
						       *sqrt(pow(myBiopolymerClassContainer.diff_x12b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] -
								       myBiopolymerClassContainer.av_d12b_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]/myBiopolymerClassContainer.counter_twobody_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count],2)));
                  
				  myBiopolymerClassContainer.dL12b_x1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]   
					  += d_x;


			  };

			  double Corr_func_12b_y = 0.0;

			  if((
						  sqrt(pow(myBiopolymerClassContainer.stamp_d12b_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] -
								  myBiopolymerClassContainer.av_d12b_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
								  /myBiopolymerClassContainer.counter_twobody_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count],2))    
						  *sqrt(pow(myBiopolymerClassContainer.diff_y12b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] -
								  myBiopolymerClassContainer.av_d12b_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]/myBiopolymerClassContainer.counter_twobody_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count],2))) > 0.0)

			  {

                  d_y = (myBiopolymerClassContainer.stamp_d12b_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] -
							  myBiopolymerClassContainer.av_d12b_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
							  /myBiopolymerClassContainer.counter_twobody_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]    
					     )*(myBiopolymerClassContainer.diff_y12b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] - myBiopolymerClassContainer.av_d12b_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
						     /myBiopolymerClassContainer.counter_twobody_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
					       )/(
						       sqrt(pow(myBiopolymerClassContainer.stamp_d12b_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] -
								       myBiopolymerClassContainer.av_d12b_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
								       /myBiopolymerClassContainer.counter_twobody_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count],2))    
						       *sqrt(pow(myBiopolymerClassContainer.diff_y12b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] -
								       myBiopolymerClassContainer.av_d12b_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]/myBiopolymerClassContainer.counter_twobody_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count],2)));
                  
				  myBiopolymerClassContainer.dL12b_y1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]   
					  += d_y;

			  };

			  double Corr_func_12b_z = 0.0;

			  if((
						  sqrt(pow(myBiopolymerClassContainer.stamp_d12b_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] -
								  myBiopolymerClassContainer.av_d12b_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
								  /myBiopolymerClassContainer.counter_twobody_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count],2))    
						  *sqrt(pow(myBiopolymerClassContainer.diff_z12b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] -
								  myBiopolymerClassContainer.av_d12b_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]/myBiopolymerClassContainer.counter_twobody_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count],2))) > 0.0)

			  {

                  d_z = (myBiopolymerClassContainer.stamp_d12b_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] -
							  myBiopolymerClassContainer.av_d12b_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
							  /myBiopolymerClassContainer.counter_twobody_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]    
					     )*(myBiopolymerClassContainer.diff_z12b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] - myBiopolymerClassContainer.av_d12b_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
						     /myBiopolymerClassContainer.counter_twobody_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
					       )/(
						       sqrt(pow(myBiopolymerClassContainer.stamp_d12b_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] -
								       myBiopolymerClassContainer.av_d12b_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
								       /myBiopolymerClassContainer.counter_twobody_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count],2))    
						       *sqrt(pow(myBiopolymerClassContainer.diff_z12b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] -
								       myBiopolymerClassContainer.av_d12b_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]/myBiopolymerClassContainer.counter_twobody_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count],2)));
                  
				  myBiopolymerClassContainer.dL12b_z1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]   
					  += d_z; 

			  };
            };
		  };
	  };
	  };

	  int res_cycle; 

	  Vec3 state_3;

	  for (int r=0;r<myParameterReader.ntc_class_container.numNTC_Torsions();r++) { 

		  if((myParameterReader.ntc_class_container.getNTC_Class(r)).threebody == 1) {

			  //  if((myParameterReader.ntc_class_container.getNTC_Class(r)).count != (myParameterReader.ntc_class_container.getNTC_Class(r+1)).count) {  

			  n22 = 0;

			  Vec3 d12,d13,d23;
			  Vec3 f_bias;
       
      String chainId1=(myParameterReader.ntc_class_container.getNTC_Class(r)).NtC_FirstBPChain; 
      
      ResidueID residueNumber1=(myParameterReader.ntc_class_container.getNTC_Class(r)).FirstBPResidue;      
      ResidueID residueNumber2=(myParameterReader.ntc_class_container.getNTC_Class(r)).sec_res_threebody;
      ResidueID residueNumber3=(myParameterReader.ntc_class_container.getNTC_Class(r)).third_res_threebody;
      
      NTC_PAR_BondRow myNTC_PAR_BondRow = myNTC_PAR_Class.myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[(myParameterReader.ntc_class_container.getNTC_Class(r)).NTC_PAR_BondRowIndex];
      
      all_refs_true = 0;
      
    if( myBiopolymerClassContainer.updBiopolymerClass(chainId1).hasAtom(residueNumber1,"C2") == true) {      
      
      body1 = myBiopolymerClassContainer.updAtomMobilizedBody(matter,chainId1,residueNumber1,"C2");
      Vec3 stationA = myBiopolymerClassContainer.getAtomLocationInMobilizedBodyFrame(chainId1,residueNumber1,"C2");
      state_1 = myBiopolymerClassContainer.calcAtomLocationInGroundFrame(state,chainId1,residueNumber1,"C2");       
      all_refs_true += 1;
      
    }
      
      String chainId2=(myParameterReader.ntc_class_container.getNTC_Class(r)).sec_res_chain_threebody;
     
    if( myBiopolymerClassContainer.updBiopolymerClass(chainId2).hasAtom(residueNumber2,"C2") == true) {      
      
      body2 = myBiopolymerClassContainer.updAtomMobilizedBody(matter,chainId2,residueNumber2,"C2");
      Vec3 stationB = myBiopolymerClassContainer.getAtomLocationInMobilizedBodyFrame(chainId2,residueNumber2,"C2");
      state_2 = myBiopolymerClassContainer.calcAtomLocationInGroundFrame(state,chainId2,residueNumber2,"C2");     
      all_refs_true += 1;
      
    }
      
      String chainId3=(myParameterReader.ntc_class_container.getNTC_Class(r)).third_res_chain_threebody;
     
    if( myBiopolymerClassContainer.updBiopolymerClass(chainId3).hasAtom(residueNumber3,"C2") == true) {      
      
      body3 = myBiopolymerClassContainer.updAtomMobilizedBody(matter,chainId3,residueNumber3,"C2");
      Vec3 stationC = myBiopolymerClassContainer.getAtomLocationInMobilizedBodyFrame(chainId3,residueNumber3,"C2");
      state_3 = myBiopolymerClassContainer.calcAtomLocationInGroundFrame(state,chainId3,residueNumber3,"C2");      
      all_refs_true += 1;     
      
    }
      
    double d_1_x, d_1_y, d_1_z;  
     
   if(all_refs_true == 3) { 
    
     myBiopolymerClassContainer.bias_angle[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = myBiopolymerClassContainer.dL_angle1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     /myBiopolymerClassContainer.counter_threebody_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];  
     
     myBiopolymerClassContainer.bias_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]  = 
     (myBiopolymerClassContainer.dL12_x1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     /myBiopolymerClassContainer.counter_threebody_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     *myBiopolymerClassContainer.dL13_x1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     /myBiopolymerClassContainer.counter_threebody_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     *myBiopolymerClassContainer.dL23_x1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     /myBiopolymerClassContainer.counter_threebody_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count])
     *myBiopolymerClassContainer.bias_angle[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]*
     myBiopolymerClassContainer.bias_3body_x_last[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];     
     
     myBiopolymerClassContainer.bias_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]  = 
     (myBiopolymerClassContainer.dL12_y1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     /myBiopolymerClassContainer.counter_threebody_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     *myBiopolymerClassContainer.dL13_y1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     /myBiopolymerClassContainer.counter_threebody_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     *myBiopolymerClassContainer.dL23_y1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     /myBiopolymerClassContainer.counter_threebody_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count])
     *myBiopolymerClassContainer.bias_angle[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]*
     myBiopolymerClassContainer.bias_3body_y_last[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];
     
     myBiopolymerClassContainer.bias_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]  = 
     (myBiopolymerClassContainer.dL12_z1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     /myBiopolymerClassContainer.counter_threebody_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     *myBiopolymerClassContainer.dL13_z1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     /myBiopolymerClassContainer.counter_threebody_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     *myBiopolymerClassContainer.dL23_z1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     /myBiopolymerClassContainer.counter_threebody_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count])
     *myBiopolymerClassContainer.bias_angle[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]*
     myBiopolymerClassContainer.bias_3body_z_last[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]; 
      
      double d_b1_x, d_b1_y, d_b1_z,d_b2_x,d_b2_y,d_b2_z,d_b3_x,d_b3_y,d_b3_z;
           
      d_t_12 = d12.norm();
      d_t_13 = d13.norm();
      d_t_23 = d23.norm();
   
      
     if((int)myBiopolymerClassContainer.counter_gen[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]%
        (int)myBiopolymerClassContainer.tau_1_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] == 0) { 
      
      range_x = -1.0;
      range_y = -1.0;
      range_z = -1.0;
      
  for(k=1;k<=100;k++) {    
      
 //    if((int)myBiopolymerClassContainer.counter_gen[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]%
 //       (int)myBiopolymerClassContainer.tau_1_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] == 0) {  
      
        p_3body_x = 0.0;
        p_3body_y = 0.0;
        p_3body_z = 0.0;
      
        myBiopolymerClassContainer.hist_global_12_3body_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][k] +=
        exp(-pow(myBiopolymerClassContainer.bias_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] - range_x,2)/0.02);          
     
        myBiopolymerClassContainer.weight_12_3body_x[k] +=
        exp(-pow(myBiopolymerClassContainer.bias_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] - range_x,2)/0.02);     
     
        if(myBiopolymerClassContainer.weight_12_3body_x[k] != 0.0) {
          
           if(myBiopolymerClassContainer.hist_global_12_3body_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][k]/myBiopolymerClassContainer.weight_12_3body_x[k] > 1E-8 &&
              myBiopolymerClassContainer.hist_global_12_3body_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][k]/myBiopolymerClassContainer.weight_12_3body_x[k] < 
              p_3body_x) {
               
              p_3body_x =
              myBiopolymerClassContainer.hist_global_12_3body_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][k]/myBiopolymerClassContainer.weight_12_3body_x[k]; 
               
            }; 
        };
 //    }
     
 //    if((int)myBiopolymerClassContainer.counter_gen[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]%
 //       (int)myBiopolymerClassContainer.tau_1_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] == 0) {  
      
        myBiopolymerClassContainer.hist_global_12_3body_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][k] +=
        exp(-pow(myBiopolymerClassContainer.bias_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] - range_y,2)/0.02);          
     
        myBiopolymerClassContainer.weight_12_3body_y[k] +=
        exp(-pow(myBiopolymerClassContainer.bias_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] - range_y,2)/0.02);     
     
        if(myBiopolymerClassContainer.weight_12_3body_y[k] != 0.0) {
          
           if(myBiopolymerClassContainer.hist_global_12_3body_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][k]/myBiopolymerClassContainer.weight_12_3body_y[k] > 1E-8 &&
              myBiopolymerClassContainer.hist_global_12_3body_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][k]/myBiopolymerClassContainer.weight_12_3body_y[k] < 
              p_3body_y) {
               
              p_3body_y =
              myBiopolymerClassContainer.hist_global_12_3body_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][k]/myBiopolymerClassContainer.weight_12_3body_y[k]; 
               
            }; 
        };     
     
 //    }
     
 //    if((int)myBiopolymerClassContainer.counter_gen[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]%
 //       (int)myBiopolymerClassContainer.tau_1_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] == 0) {  
      
        myBiopolymerClassContainer.hist_global_12_3body_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][k] +=
        exp(-pow(myBiopolymerClassContainer.bias_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] - range_z,2)/0.02);          
     
        myBiopolymerClassContainer.weight_12_3body_z[k] +=
        exp(-pow(myBiopolymerClassContainer.bias_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] - range_z,2)/0.02);     
     
        if(myBiopolymerClassContainer.weight_12_3body_z[k] != 0.0) {
          
           if(myBiopolymerClassContainer.hist_global_12_3body_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][k]/myBiopolymerClassContainer.weight_12_3body_z[k] > 1E-8 &&
              myBiopolymerClassContainer.hist_global_12_3body_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][k]/myBiopolymerClassContainer.weight_12_3body_z[k] < 
              p_3body_z) {
               
              p_3body_z =
              myBiopolymerClassContainer.hist_global_12_3body_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][k]/myBiopolymerClassContainer.weight_12_3body_z[k]; 
               
            }; 
        };     
     
  //   }    
        
         ran = rand()%100/100.0;
        
         if(ran < (p_3body_x + p_3body_y + p_3body_z)/3.0) {
              
                   myBiopolymerClassContainer.reset_3body[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 1;
                  
              } else {
                
                   myBiopolymerClassContainer.reset_3body[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 0;
                  
         };        
       
        if(sqrt(pow(myBiopolymerClassContainer.bias_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] - range_x,2)) < 0.02) {
            myBiopolymerClassContainer.meas_12_3body_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = k;         
        };
        if(sqrt(pow(myBiopolymerClassContainer.bias_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] - range_y,2)) < 0.02) {
            myBiopolymerClassContainer.meas_12_3body_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = k;         
        };
        if(sqrt(pow(myBiopolymerClassContainer.bias_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] - range_z,2)) < 0.02) {
            myBiopolymerClassContainer.meas_12_3body_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = k;            
        };
        
        range_x += 0.02;
        range_y += 0.02;
        range_z += 0.02;
     };

     /*
     
if((int)myBiopolymerClassContainer.counter_gen[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]%
  ((int)myBiopolymerClassContainer.tau_1_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]) == 0) {     
            
     ran = rand()%100/100.0;
     
        if( ran <= 1.0/(1.0 + exp(-myBiopolymerClassContainer.bias_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]))) {            
            
             myBiopolymerClassContainer.tau_1_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = myBiopolymerClassContainer.counter_threebody_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];
             myBiopolymerClassContainer.counter_threebody_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 1.0;
             
             myBiopolymerClassContainer.dL12_x1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 0.0;             
             myBiopolymerClassContainer.dL13_x1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 0.0;
             myBiopolymerClassContainer.dL23_x1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 0.0;         
                   
      for(k=1;k<=100;k++) {
             myBiopolymerClassContainer.av_d12_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]    = 0.0;             
             myBiopolymerClassContainer.av_d13_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]    = 0.0;              
             myBiopolymerClassContainer.av_d23_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]    = 0.0;               
             
             myBiopolymerClassContainer.beta_2_x_threebody_last[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 0.0;
             
             int n33;
             
             for(n33=1;n33<=100;n33++) {
                
                 myBiopolymerClassContainer.hist_global_12_3body_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][n33] = 0.0;
                 myBiopolymerClassContainer.weight_12_3body_x[n33] = 0.0;              
                 myBiopolymerClassContainer.alpha_2_x_threebody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][n33] = 0.0;
                 myBiopolymerClassContainer.probability_x3_threebody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][n33] = 0.0;
              };
      };
     };
     
if((int)myBiopolymerClassContainer.counter_gen[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]%
  ((int)myBiopolymerClassContainer.tau_1_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]) == 0) {     
          
     ran = rand()%100/100.0;
     
        if( ran <= 1.0/(1.0 + exp(-myBiopolymerClassContainer.bias_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]))) {
            
             myBiopolymerClassContainer.tau_1_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = myBiopolymerClassContainer.counter_threebody_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];
             myBiopolymerClassContainer.counter_threebody_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 1.0;
             
             myBiopolymerClassContainer.dL12_y1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 0.0;             
             myBiopolymerClassContainer.dL13_y1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 0.0;
             myBiopolymerClassContainer.dL23_y1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 0.0;         
             
             myBiopolymerClassContainer.av_d12_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]    = 0.0;             
             myBiopolymerClassContainer.av_d13_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]    = 0.0;              
             myBiopolymerClassContainer.av_d23_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]    = 0.0;               
             
             myBiopolymerClassContainer.beta_2_y_threebody_last[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 0.0;
             
             int n33;
             
             for(n33=1;n33<=100;n33++) {
                
                 myBiopolymerClassContainer.hist_global_12_3body_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][n33] = 0.0;
                 myBiopolymerClassContainer.weight_12_3body_y[n33] = 0.0;              
                 myBiopolymerClassContainer.alpha_2_y_threebody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][n33] = 0.0;
                 myBiopolymerClassContainer.probability_y3_threebody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][n33] = 0.0;
                 
              };
      };
     };     
     
if((int)myBiopolymerClassContainer.counter_gen[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]%
  ((int)myBiopolymerClassContainer.tau_1_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]) == 0) {     
     
     ran = rand()%100/100.0;
     
        if( ran <= 1.0/(1.0 + exp(-myBiopolymerClassContainer.bias_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]))) {
            
             myBiopolymerClassContainer.tau_1_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = myBiopolymerClassContainer.counter_threebody_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];
             myBiopolymerClassContainer.counter_threebody_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 1.0;
             
             myBiopolymerClassContainer.dL12_z1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 0.0;             
             myBiopolymerClassContainer.dL13_z1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 0.0;
             myBiopolymerClassContainer.dL23_z1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 0.0;         
             
             myBiopolymerClassContainer.av_d12_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]    = 0.0;             
             myBiopolymerClassContainer.av_d13_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]    = 0.0;              
             myBiopolymerClassContainer.av_d23_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]    = 0.0;               
             
             myBiopolymerClassContainer.beta_2_z_threebody_last[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 0.0;
             
             int n33;
             
             for(n33=1;n33<=100;n33++) {
                
                 myBiopolymerClassContainer.hist_global_12_3body_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][n33] = 0.0;
                 myBiopolymerClassContainer.weight_12_3body_z[n33] = 0.0;              
                 myBiopolymerClassContainer.alpha_2_z_threebody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][n33] = 0.0;
                 myBiopolymerClassContainer.probability_z3_threebody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][n33] = 0.0;
                 
              };
      }; */
      };
     };      
    };
   };
   
for (int r=0;r<myParameterReader.ntc_class_container.numNTC_Torsions();r++) { 
      
  if((myParameterReader.ntc_class_container.getNTC_Class(r)).twobody == 1) {
      
 // if((myParameterReader.ntc_class_container.getNTC_Class(r)).count != (myParameterReader.ntc_class_container.getNTC_Class(r+1)).count) {
      
      Vec3 d12;
      Vec3 f_bias; 
      
      n22 = 0;
      
      String chainId1=(myParameterReader.ntc_class_container.getNTC_Class(r)).NtC_FirstBPChain; 
      
      NTC_PAR_BondRow myNTC_PAR_BondRow = myNTC_PAR_Class.myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[(myParameterReader.ntc_class_container.getNTC_Class(r)).NTC_PAR_BondRowIndex];
      
      ResidueID residueNumber1=(myParameterReader.ntc_class_container.getNTC_Class(r)).FirstBPResidue;      
      ResidueID residueNumber2=(myParameterReader.ntc_class_container.getNTC_Class(r)).sec_res_twobody;
      
      all_refs_true = 0;
      
     if( myBiopolymerClassContainer.updBiopolymerClass(chainId1).hasAtom(residueNumber1,"C2") == true) {     
      
      body1 = myBiopolymerClassContainer.updAtomMobilizedBody(matter,chainId1,residueNumber1,"C2");
      Vec3 stationA = myBiopolymerClassContainer.getAtomLocationInMobilizedBodyFrame(chainId1,residueNumber1,"C2");
      state_1 = myBiopolymerClassContainer.calcAtomLocationInGroundFrame(state,chainId1,residueNumber1,"C2");       
      all_refs_true += 1;
      
     }
      
      String chainId2=(myParameterReader.ntc_class_container.getNTC_Class(r)).sec_res_chain_twobody;
      
    if( myBiopolymerClassContainer.updBiopolymerClass(chainId2).hasAtom(residueNumber2,"C2") == true) {      
      
      body2 = myBiopolymerClassContainer.updAtomMobilizedBody(matter,chainId2,residueNumber2,"C2");
      Vec3 stationB = myBiopolymerClassContainer.getAtomLocationInMobilizedBodyFrame(chainId2,residueNumber2,"C2");
      state_2 = myBiopolymerClassContainer.calcAtomLocationInGroundFrame(state,chainId2,residueNumber2,"C2");       
      all_refs_true += 1;      
      
    }
      
   //   Vec3 d12 = state_1%state_2;      
   
   if(all_refs_true == 2) {
   
      
     myBiopolymerClassContainer.bias_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]  = myBiopolymerClassContainer.dL12b_x1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     /myBiopolymerClassContainer.counter_twobody_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]*
     myBiopolymerClassContainer.bias_2body_x_last[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];
     myBiopolymerClassContainer.bias_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]  = myBiopolymerClassContainer.dL12b_y1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     /myBiopolymerClassContainer.counter_twobody_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]*
     myBiopolymerClassContainer.bias_2body_y_last[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];
     myBiopolymerClassContainer.bias_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]  = myBiopolymerClassContainer.dL12b_z1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     /myBiopolymerClassContainer.counter_twobody_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]*
     myBiopolymerClassContainer.bias_2body_z_last[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];      
     
      double d_b1_x, d_b1_y, d_b1_z,d_b2_x,d_b2_y,d_b2_z,d_b3_x,d_b3_y,d_b3_z;
     
      d_t_12 = d12.norm();
      
    if( (int)myBiopolymerClassContainer.counter_gen[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]%(int)myBiopolymerClassContainer.tau_1_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] == 0) {   
        
      range_x = -1.0;
      range_y = -1.0;
      range_z = -1.0;
      
      for(k=1;k<=100;k++) {
          
        p_2body_x = 0.0;
        p_2body_y = 0.0;
        p_2body_z = 0.0;
     
 //    if((int)myBiopolymerClassContainer.counter_gen[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]%
 //       (int)myBiopolymerClassContainer.tau_2_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] == 0) {  
      
        myBiopolymerClassContainer.hist_global_12_2body_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][k] +=
        exp(-pow(myBiopolymerClassContainer.bias_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] - range_x,2)/0.02);          
     
        myBiopolymerClassContainer.weight_12_2body_x[k] +=
        exp(-pow(myBiopolymerClassContainer.bias_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] - range_x,2)/0.02); 
     
        if(myBiopolymerClassContainer.weight_12_2body_x[k] != 0.0) {
          
           if(myBiopolymerClassContainer.hist_global_12_2body_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][k]/myBiopolymerClassContainer.weight_12_2body_x[k] > 1E-8 &&
              myBiopolymerClassContainer.hist_global_12_2body_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][k]/myBiopolymerClassContainer.weight_12_2body_x[k] < 
              p_2body_x) {
               
              p_2body_x =
              myBiopolymerClassContainer.hist_global_12_2body_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][k]/myBiopolymerClassContainer.weight_12_2body_x[k]; 
               
            }; 
        };     
     
 //    }
     
 //    if((int)myBiopolymerClassContainer.counter_gen[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]%
 //       (int)myBiopolymerClassContainer.tau_2_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] == 0) {  
      
        myBiopolymerClassContainer.hist_global_12_2body_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][k] +=
        exp(-pow(myBiopolymerClassContainer.bias_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] - range_y,2)/0.02);          
     
        myBiopolymerClassContainer.weight_12_2body_y[k] +=
        exp(-pow(myBiopolymerClassContainer.bias_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] - range_y,2)/0.02);     
     
        if(myBiopolymerClassContainer.weight_12_2body_y[k] != 0.0) {
          
           if(myBiopolymerClassContainer.hist_global_12_2body_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][k]/myBiopolymerClassContainer.weight_12_2body_y[k] > 1E-8 &&
              myBiopolymerClassContainer.hist_global_12_2body_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][k]/myBiopolymerClassContainer.weight_12_2body_y[k] < 
              p_2body_y) {
               
              p_2body_y =
              myBiopolymerClassContainer.hist_global_12_2body_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][k]/myBiopolymerClassContainer.weight_12_2body_y[k]; 
               
            }; 
        };       
     
  //   }
     
  //   if((int)myBiopolymerClassContainer.counter_gen[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]%
  //      (int)myBiopolymerClassContainer.tau_2_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] == 0) {  
      
        myBiopolymerClassContainer.hist_global_12_2body_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][k] +=
        exp(-pow(myBiopolymerClassContainer.bias_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] - range_z,2)/0.02);          
     
        myBiopolymerClassContainer.weight_12_2body_z[k] +=
        exp(-pow(myBiopolymerClassContainer.bias_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] - range_z,2)/0.02);     
     
        if(myBiopolymerClassContainer.weight_12_2body_z[k] != 0.0) {
          
           if(myBiopolymerClassContainer.hist_global_12_2body_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][k]/myBiopolymerClassContainer.weight_12_2body_z[k] > 1E-8 &&
              myBiopolymerClassContainer.hist_global_12_2body_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][k]/myBiopolymerClassContainer.weight_12_2body_z[k] < 
              p_2body_z) {
               
              p_2body_z =
              myBiopolymerClassContainer.hist_global_12_2body_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][k]/myBiopolymerClassContainer.weight_12_2body_z[k]; 
               
            }; 
        };       
     
         ran = rand()%100/100.0;
        
         if(ran < (p_2body_x + p_2body_y + p_2body_z)/3.0) {
              
                   myBiopolymerClassContainer.reset_2body[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 1;
                  
              } else {
                
                   myBiopolymerClassContainer.reset_2body[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 0;
                  
         };          
        
   //  }        
      
        if(sqrt(pow(myBiopolymerClassContainer.bias_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] - range_x,2)) < 0.02) {
            myBiopolymerClassContainer.meas_12_2body_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = k;        
        };
        if(sqrt(pow(myBiopolymerClassContainer.bias_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] - range_y,2)) < 0.02) {
            myBiopolymerClassContainer.meas_12_2body_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = k;           
        };
        if(sqrt(pow(myBiopolymerClassContainer.bias_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] - range_z,2)) < 0.02) {
            myBiopolymerClassContainer.meas_12_2body_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = k;             
        };
        
        range_x += 0.02;
        range_y += 0.02;
        range_z += 0.02;
        
  //   };        
   };


/*   
   
if((int)myBiopolymerClassContainer.counter_gen[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]%
  ((int)myBiopolymerClassContainer.tau_2_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]) == 0) {     
     
     ran = rand()%100/100.0;
     
        if( ran <= 1.0/(1.0 + exp(-myBiopolymerClassContainer.bias_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]))) {             
            
             myBiopolymerClassContainer.tau_2_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = myBiopolymerClassContainer.counter_twobody_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];
             myBiopolymerClassContainer.counter_twobody_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 1.0;
             
             myBiopolymerClassContainer.dL12b_x1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 0.0;                   
             myBiopolymerClassContainer.av_d12b_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]    = 0.0;                          
             
             myBiopolymerClassContainer.beta_2_x_twobody_last[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 0.0;
             
             int n33;
             
             for(n33=1;n33<=100;n33++) {
                
                 myBiopolymerClassContainer.hist_global_12_2body_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][n33] = 0.0;
                 myBiopolymerClassContainer.weight_12_2body_x[n33] = 0.0;              
                 myBiopolymerClassContainer.alpha_2_x_twobody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][n33] = 0.0;
                 myBiopolymerClassContainer.probability_x3_twobody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][n33] = 0.0;
                 
              };
            };     
           };
   
   
if((int)myBiopolymerClassContainer.counter_gen[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]%
  ((int)myBiopolymerClassContainer.tau_2_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]) == 0) {     
     
     ran = rand()%100/100.0;
     
        if( ran <= 1.0/(1.0 + exp(-myBiopolymerClassContainer.bias_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]))) {        
            
             myBiopolymerClassContainer.tau_2_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = myBiopolymerClassContainer.counter_twobody_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];
             myBiopolymerClassContainer.counter_twobody_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 1.0;
             
             myBiopolymerClassContainer.dL12b_y1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 0.0;                   
             myBiopolymerClassContainer.av_d12b_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]    = 0.0;                          
             
             myBiopolymerClassContainer.beta_2_y_twobody_last[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 0.0;
             
             int n33;
             
             for(n33=1;n33<=100;n33++) {
                
                 myBiopolymerClassContainer.hist_global_12_2body_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][n33] = 0.0;
                 myBiopolymerClassContainer.weight_12_2body_y[n33] = 0.0;              
                 myBiopolymerClassContainer.alpha_2_y_twobody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][n33] = 0.0;
                 myBiopolymerClassContainer.probability_y3_twobody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][n33] = 0.0;
                 
              };
      };
     };     
   
     
if((int)myBiopolymerClassContainer.counter_gen[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]%
  ((int)myBiopolymerClassContainer.tau_2_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]) == 0) {     
     
     ran = rand()%100/100.0;
     
        if( ran <= 1.0/(1.0 + exp(-myBiopolymerClassContainer.bias_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]))) {
            
             myBiopolymerClassContainer.tau_2_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = myBiopolymerClassContainer.counter_twobody_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];
             myBiopolymerClassContainer.counter_twobody_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 1.0;
             
             myBiopolymerClassContainer.dL12b_z1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 0.0;                   
             myBiopolymerClassContainer.av_d12b_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 0.0;                          
             myBiopolymerClassContainer.beta_2_z_twobody_last[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 0.0;
             
             int n33;
             
             for(n33=1;n33<=100;n33++) {
                
                 myBiopolymerClassContainer.hist_global_12_2body_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][n33] = 0.0;
                 myBiopolymerClassContainer.weight_12_2body_z[n33] = 0.0;              
                 myBiopolymerClassContainer.alpha_2_z_twobody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][n33] = 0.0;
                 myBiopolymerClassContainer.probability_z3_twobody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][n33] = 0.0;
                 
              }; */
          }; 
        };  
    };
};
    
double delta_l, delta_b, delta_l2, delta_b2;

for (int r=0;r<myParameterReader.ntc_class_container.numNTC_Torsions();r++) { 

          
   if((myParameterReader.ntc_class_container.getNTC_Class(r)).threebody == 1) {

     myBiopolymerClassContainer.bias_angle[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = myBiopolymerClassContainer.dL_angle1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     /myBiopolymerClassContainer.counter_threebody_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];  
     
     myBiopolymerClassContainer.bias_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]  = 
     (myBiopolymerClassContainer.dL12_x1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     /myBiopolymerClassContainer.counter_threebody_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     *myBiopolymerClassContainer.dL13_x1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     /myBiopolymerClassContainer.counter_threebody_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     *myBiopolymerClassContainer.dL23_x1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     /myBiopolymerClassContainer.counter_threebody_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count])
     *myBiopolymerClassContainer.bias_angle[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]*
     myBiopolymerClassContainer.bias_3body_x_last[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];     
     
     myBiopolymerClassContainer.bias_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]  = 
     (myBiopolymerClassContainer.dL12_y1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     /myBiopolymerClassContainer.counter_threebody_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     *myBiopolymerClassContainer.dL13_y1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     /myBiopolymerClassContainer.counter_threebody_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     *myBiopolymerClassContainer.dL23_y1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     /myBiopolymerClassContainer.counter_threebody_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count])
     *myBiopolymerClassContainer.bias_angle[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]*
     myBiopolymerClassContainer.bias_3body_y_last[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];
     
     myBiopolymerClassContainer.bias_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]  = 
     (myBiopolymerClassContainer.dL12_z1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     /myBiopolymerClassContainer.counter_threebody_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     *myBiopolymerClassContainer.dL13_z1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     /myBiopolymerClassContainer.counter_threebody_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     *myBiopolymerClassContainer.dL23_z1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     /myBiopolymerClassContainer.counter_threebody_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count])
     *myBiopolymerClassContainer.bias_angle[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]*
     myBiopolymerClassContainer.bias_3body_z_last[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]; 
     
   };   
    
   if((myParameterReader.ntc_class_container.getNTC_Class(r)).twobody == 1) {

     myBiopolymerClassContainer.bias_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]  = myBiopolymerClassContainer.dL12b_x1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     /myBiopolymerClassContainer.counter_twobody_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]*
     myBiopolymerClassContainer.bias_2body_x_last[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];
     myBiopolymerClassContainer.bias_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]  = myBiopolymerClassContainer.dL12b_y1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     /myBiopolymerClassContainer.counter_twobody_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]*
     myBiopolymerClassContainer.bias_2body_y_last[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];
     myBiopolymerClassContainer.bias_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]  = myBiopolymerClassContainer.dL12b_z1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     /myBiopolymerClassContainer.counter_twobody_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]*
     myBiopolymerClassContainer.bias_2body_z_last[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];      
     
  };    
};

  
for (int r=0;r<myParameterReader.ntc_class_container.numNTC_Torsions();r++) { 
    
 //if((myParameterReader.ntc_class_container.getNTC_Class(r)).count != (myParameterReader.ntc_class_container.getNTC_Class(r+1)).count) {   
 
  if((myParameterReader.ntc_class_container.getNTC_Class(r)).threebody == 1) {  
  
   for(n22 = 0; n22 <= 0; n22 ++) {   
      
      Vec3 d12,d13,d23;
      Vec3 f_bias;      
      
      String chainId1=(myParameterReader.ntc_class_container.getNTC_Class(r)).NtC_FirstBPChain; 
      
      ResidueID residueNumber1=(myParameterReader.ntc_class_container.getNTC_Class(r)).FirstBPResidue;      
      ResidueID residueNumber2=(myParameterReader.ntc_class_container.getNTC_Class(r)).sec_res_threebody; 
      ResidueID residueNumber3=(myParameterReader.ntc_class_container.getNTC_Class(r)).third_res_threebody; 
      
      NTC_PAR_BondRow myNTC_PAR_BondRow = myNTC_PAR_Class.myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[(myParameterReader.ntc_class_container.getNTC_Class(r)).NTC_PAR_BondRowIndex];
      
      all_refs_true = 0;
      
    if( myBiopolymerClassContainer.updBiopolymerClass(chainId1).hasAtom(residueNumber1,"C2") == true) { 
      
      body1 = myBiopolymerClassContainer.updAtomMobilizedBody(matter,chainId1,residueNumber1,"C2");
      Vec3 stationA = myBiopolymerClassContainer.getAtomLocationInMobilizedBodyFrame(chainId1,residueNumber1,"C2");
      state_1 = myBiopolymerClassContainer.calcAtomLocationInGroundFrame(state,chainId1,residueNumber1,"C2");       
      all_refs_true += 1;
      
    }
      
      String chainId2=(myParameterReader.ntc_class_container.getNTC_Class(r)).sec_res_chain_threebody;
     
    if( myBiopolymerClassContainer.updBiopolymerClass(chainId2).hasAtom(residueNumber2,"C2") == true) {      
      
      body2 = myBiopolymerClassContainer.updAtomMobilizedBody(matter,chainId2,residueNumber2,"C2");
      Vec3 stationB = myBiopolymerClassContainer.getAtomLocationInMobilizedBodyFrame(chainId2,residueNumber2,"C2");
      state_2 = myBiopolymerClassContainer.calcAtomLocationInGroundFrame(state,chainId2,residueNumber2,"C2");     
      all_refs_true += 1;
      
    }
      
      String chainId3=(myParameterReader.ntc_class_container.getNTC_Class(r)).third_res_chain_threebody;
     
    if( myBiopolymerClassContainer.updBiopolymerClass(chainId3).hasAtom(residueNumber3,"C2") == true) {      
      
      body3 = myBiopolymerClassContainer.updAtomMobilizedBody(matter,chainId3,residueNumber3,"C2");
      Vec3 stationC = myBiopolymerClassContainer.getAtomLocationInMobilizedBodyFrame(chainId3,residueNumber3,"C2");
      state_3 = myBiopolymerClassContainer.calcAtomLocationInGroundFrame(state,chainId3,residueNumber3,"C2");          
      all_refs_true += 1;
            
    }
      
 // 12
 
  myBiopolymerClassContainer.h_div_w_3body_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 0.0;
  myBiopolymerClassContainer.h_div_w_3body_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 0.0;  
  myBiopolymerClassContainer.h_div_w_3body_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 0.0; 
 
 if(all_refs_true == 3) {
 
     myBiopolymerClassContainer.bias_angle[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = myBiopolymerClassContainer.dL_angle1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     /myBiopolymerClassContainer.counter_threebody_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];  
     
     myBiopolymerClassContainer.bias_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]  = 
     (myBiopolymerClassContainer.dL12_x1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     /myBiopolymerClassContainer.counter_threebody_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     *myBiopolymerClassContainer.dL13_x1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     /myBiopolymerClassContainer.counter_threebody_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     *myBiopolymerClassContainer.dL23_x1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     /myBiopolymerClassContainer.counter_threebody_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count])
     *myBiopolymerClassContainer.bias_angle[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]*
     myBiopolymerClassContainer.bias_3body_x_last[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];     
     
     myBiopolymerClassContainer.bias_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]  = 
     (myBiopolymerClassContainer.dL12_y1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     /myBiopolymerClassContainer.counter_threebody_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     *myBiopolymerClassContainer.dL13_y1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     /myBiopolymerClassContainer.counter_threebody_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     *myBiopolymerClassContainer.dL23_y1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     /myBiopolymerClassContainer.counter_threebody_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count])
     *myBiopolymerClassContainer.bias_angle[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]*
     myBiopolymerClassContainer.bias_3body_y_last[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];
     
     myBiopolymerClassContainer.bias_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]  = 
     (myBiopolymerClassContainer.dL12_z1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     /myBiopolymerClassContainer.counter_threebody_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     *myBiopolymerClassContainer.dL13_z1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     /myBiopolymerClassContainer.counter_threebody_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     *myBiopolymerClassContainer.dL23_z1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     /myBiopolymerClassContainer.counter_threebody_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count])
     *myBiopolymerClassContainer.bias_angle[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]*
     myBiopolymerClassContainer.bias_3body_z_last[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];  
      
 if(myBiopolymerClassContainer.meas_12_3body_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] > 1 &&
    myBiopolymerClassContainer.meas_12_3body_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] <= 100 &&
    myBiopolymerClassContainer.meas_12_3body_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] > 1 &&
    myBiopolymerClassContainer.meas_12_3body_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] <= 100 && 
    myBiopolymerClassContainer.meas_12_3body_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] > 1 &&
    myBiopolymerClassContainer.meas_12_3body_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] <= 100
) {
     
   f_bias[0] = 0.0;
   f_bias[1] = 0.0;
   f_bias[2] = 0.0;
   
   Vec3 random_force;
      
   if( myBiopolymerClassContainer.hist_global_12_3body_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][myBiopolymerClassContainer.meas_12_3body_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]] != 0.0 &&
       myBiopolymerClassContainer.weight_12_3body_x[myBiopolymerClassContainer.meas_12_3body_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]] != 0.0 &&
       myBiopolymerClassContainer.delta_x12[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] != 0.0 && p_3body_x != 0.0 &&
       myBiopolymerClassContainer.hist_global_12_3body_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][myBiopolymerClassContainer.meas_12_3body_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]-1] != 0.0 &&
       myBiopolymerClassContainer.weight_12_3body_x[myBiopolymerClassContainer.meas_12_3body_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]-1] != 0.0)
       
   {
       
       d_x = (-log(((myBiopolymerClassContainer.hist_global_12_3body_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][myBiopolymerClassContainer.meas_12_3body_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]]/
           myBiopolymerClassContainer.weight_12_3body_x[myBiopolymerClassContainer.meas_12_3body_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]]
                    ))/p_3body_x) - myBiopolymerClassContainer.hist_global_12_3body_x_last[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]);
       
       myBiopolymerClassContainer.h_div_w_3body_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = myBiopolymerClassContainer.hist_global_12_3body_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][myBiopolymerClassContainer.meas_12_3body_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]]/
           myBiopolymerClassContainer.weight_12_3body_x[myBiopolymerClassContainer.meas_12_3body_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]];
       
  myBiopolymerClassContainer.l_max_t[0]+= myBiopolymerClassContainer.hist_global_12_3body_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][myBiopolymerClassContainer.meas_12_3body_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]]/myBiopolymerClassContainer.weight_12_3body_x[myBiopolymerClassContainer.meas_12_3body_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]];       
       
   myBiopolymerClassContainer.counter_l_max_t[0]+= 1.0;       
   
       if(d_x != 0.0) {
           
           myBiopolymerClassContainer.bias_3body_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] += myBiopolymerClassContainer.delta_x12[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]/d_x; 
        
           myBiopolymerClassContainer.f_bias_3body_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = myBiopolymerClassContainer.delta_x12[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]/d_x;
           
       }
    
       myBiopolymerClassContainer.hist_global_12_3body_x_last[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] =
       -log(((myBiopolymerClassContainer.hist_global_12_3body_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][myBiopolymerClassContainer.meas_12_3body_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]]/
           myBiopolymerClassContainer.weight_12_3body_x[myBiopolymerClassContainer.meas_12_3body_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]]
                    ))/p_3body_x); 
       
   };
           
   if( myBiopolymerClassContainer.hist_global_12_3body_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][myBiopolymerClassContainer.meas_12_3body_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]] != 0.0 &&
       myBiopolymerClassContainer.weight_12_3body_y[myBiopolymerClassContainer.meas_12_3body_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]] != 0.0 &&
       myBiopolymerClassContainer.delta_y12[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] != 0.0 && p_3body_y != 0.0 &&
       myBiopolymerClassContainer.hist_global_12_3body_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][myBiopolymerClassContainer.meas_12_3body_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]-1] != 0.0 &&
       myBiopolymerClassContainer.weight_12_3body_y[myBiopolymerClassContainer.meas_12_3body_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]-1] != 0.0)
       
   {
       
       d_y = (-log(((myBiopolymerClassContainer.hist_global_12_3body_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][myBiopolymerClassContainer.meas_12_3body_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]]/
           myBiopolymerClassContainer.weight_12_3body_y[myBiopolymerClassContainer.meas_12_3body_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]]
                    ))/p_3body_y) - myBiopolymerClassContainer.hist_global_12_3body_y_last[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]);
       
       myBiopolymerClassContainer.h_div_w_3body_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = myBiopolymerClassContainer.hist_global_12_3body_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][myBiopolymerClassContainer.meas_12_3body_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]]/
           myBiopolymerClassContainer.weight_12_3body_y[myBiopolymerClassContainer.meas_12_3body_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]];       
       
  myBiopolymerClassContainer.l_max_t[0]+= myBiopolymerClassContainer.hist_global_12_3body_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][myBiopolymerClassContainer.meas_12_3body_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]]/myBiopolymerClassContainer.weight_12_3body_y[myBiopolymerClassContainer.meas_12_3body_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]];       
       
   myBiopolymerClassContainer.counter_l_max_t[0]+= 1.0;       
   
       if(d_y != 0.0) {
           
           myBiopolymerClassContainer.bias_3body_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] +=
           myBiopolymerClassContainer.delta_y12[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]/d_y;
           
           myBiopolymerClassContainer.f_bias_3body_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = myBiopolymerClassContainer.delta_y12[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]/d_y;     
    
       }
           
       myBiopolymerClassContainer.hist_global_12_3body_y_last[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] =
       -log(((myBiopolymerClassContainer.hist_global_12_3body_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][myBiopolymerClassContainer.meas_12_3body_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]]/
           myBiopolymerClassContainer.weight_12_3body_y[myBiopolymerClassContainer.meas_12_3body_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]]
                    ))/p_3body_y);   
       
   };   
                

   if( myBiopolymerClassContainer.hist_global_12_3body_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][myBiopolymerClassContainer.meas_12_3body_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]] != 0.0 &&
       myBiopolymerClassContainer.weight_12_3body_z[myBiopolymerClassContainer.meas_12_3body_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]] != 0.0 &&
       myBiopolymerClassContainer.delta_z12[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] != 0.0 && p_3body_z != 0.0 &&
       myBiopolymerClassContainer.hist_global_12_3body_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][myBiopolymerClassContainer.meas_12_3body_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]-1] != 0.0 &&
       myBiopolymerClassContainer.weight_12_3body_z[myBiopolymerClassContainer.meas_12_3body_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]-1] != 0.0)
       
   {
       
       d_z = (-log(((myBiopolymerClassContainer.hist_global_12_3body_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][myBiopolymerClassContainer.meas_12_3body_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]]/
           myBiopolymerClassContainer.weight_12_3body_z[myBiopolymerClassContainer.meas_12_3body_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]]
                    ))/p_3body_z) - myBiopolymerClassContainer.hist_global_12_3body_z_last[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]);
       
       myBiopolymerClassContainer.h_div_w_3body_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = myBiopolymerClassContainer.hist_global_12_3body_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][myBiopolymerClassContainer.meas_12_3body_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]]/
           myBiopolymerClassContainer.weight_12_3body_z[myBiopolymerClassContainer.meas_12_3body_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]];       
       
  myBiopolymerClassContainer.l_max_t[0]+= myBiopolymerClassContainer.hist_global_12_3body_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][myBiopolymerClassContainer.meas_12_3body_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]]/myBiopolymerClassContainer.weight_12_3body_z[myBiopolymerClassContainer.meas_12_3body_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]];       
     
   myBiopolymerClassContainer.counter_l_max_t[0]+= 1.0;   
   
       if(d_z != 0.0) {
        
           myBiopolymerClassContainer.bias_3body_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] +=
           myBiopolymerClassContainer.delta_z12[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]/d_z;           
           
           myBiopolymerClassContainer.f_bias_3body_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = myBiopolymerClassContainer.delta_z12[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]/d_z;     
    
       }
           
       myBiopolymerClassContainer.hist_global_12_3body_z_last[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] =
       -log(((myBiopolymerClassContainer.hist_global_12_3body_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][myBiopolymerClassContainer.meas_12_3body_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]]/
           myBiopolymerClassContainer.weight_12_3body_z[myBiopolymerClassContainer.meas_12_3body_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]]
                    ))/p_3body_z); 
       
   };   
   
   myBiopolymerClassContainer.bias_3body_norm[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 0.0;
   
   d_tt =
   sqrt(pow(myBiopolymerClassContainer.bias_3body_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count],2) +
        pow(myBiopolymerClassContainer.bias_3body_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count],2) +
        pow(myBiopolymerClassContainer.bias_3body_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count],2));
   
    if(d_tt > myBiopolymerClassContainer.bias_3body_norm[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]) {
       
       myBiopolymerClassContainer.bias_3body_norm[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = d_tt; 
        
     };
    }
   }
  }
 }
}
   
for (int r=0;r<myParameterReader.ntc_class_container.numNTC_Torsions();r++) { 
    
 //if((myParameterReader.ntc_class_container.getNTC_Class(r)).count != (myParameterReader.ntc_class_container.getNTC_Class(r+1)).count) {   
 
  if((myParameterReader.ntc_class_container.getNTC_Class(r)).threebody == 1) {  
  
   for(n22 = 0; n22 <= 0; n22 ++) {   
      
      Vec3 d12,d13,d23;
      Vec3 f_bias;
      
      String chainId1=(myParameterReader.ntc_class_container.getNTC_Class(r)).NtC_FirstBPChain; 
      
      ResidueID residueNumber1=(myParameterReader.ntc_class_container.getNTC_Class(r)).FirstBPResidue;      
      ResidueID residueNumber2=(myParameterReader.ntc_class_container.getNTC_Class(r)).sec_res_threebody; 
      ResidueID residueNumber3=(myParameterReader.ntc_class_container.getNTC_Class(r)).third_res_threebody; 
      
      NTC_PAR_BondRow myNTC_PAR_BondRow = myNTC_PAR_Class.myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[(myParameterReader.ntc_class_container.getNTC_Class(r)).NTC_PAR_BondRowIndex];
      
      all_refs_true = 0;
      
    if( myBiopolymerClassContainer.updBiopolymerClass(chainId1).hasAtom(residueNumber1,"C2") == true) {
      
      body1 = myBiopolymerClassContainer.updAtomMobilizedBody(matter,chainId1,residueNumber1,"C2");
      Vec3 stationA = myBiopolymerClassContainer.getAtomLocationInMobilizedBodyFrame(chainId1,residueNumber1,"C2");
      state_1 = myBiopolymerClassContainer.calcAtomLocationInGroundFrame(state,chainId1,residueNumber1,"C2");       
      all_refs_true += 1;
      }
      
      String chainId2=(myParameterReader.ntc_class_container.getNTC_Class(r)).sec_res_chain_threebody;
      
    if( myBiopolymerClassContainer.updBiopolymerClass(chainId2).hasAtom(residueNumber2,"C2") == true) {      
      
      body2 = myBiopolymerClassContainer.updAtomMobilizedBody(matter,chainId2,residueNumber2,"C2");
      Vec3 stationB = myBiopolymerClassContainer.getAtomLocationInMobilizedBodyFrame(chainId2,residueNumber2,"C2");
      state_2 = myBiopolymerClassContainer.calcAtomLocationInGroundFrame(state,chainId2,residueNumber2,"C2");     
      all_refs_true += 1;
    }
      
      String chainId3=(myParameterReader.ntc_class_container.getNTC_Class(r)).third_res_chain_threebody;
      
    if( myBiopolymerClassContainer.updBiopolymerClass(chainId3).hasAtom(residueNumber3,"C2") == true) {    
      
      body3 = myBiopolymerClassContainer.updAtomMobilizedBody(matter,chainId3,residueNumber3,"C2");
      Vec3 stationC = myBiopolymerClassContainer.getAtomLocationInMobilizedBodyFrame(chainId3,residueNumber3,"C2");
      state_3 = myBiopolymerClassContainer.calcAtomLocationInGroundFrame(state,chainId3,residueNumber3,"C2");          
      all_refs_true += 1;
    }
      
 // 12
 
  if(all_refs_true == 3) {
 
     myBiopolymerClassContainer.bias_angle[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = myBiopolymerClassContainer.dL_angle1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     /myBiopolymerClassContainer.counter_threebody_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];  
     
     myBiopolymerClassContainer.bias_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]  = 
     (myBiopolymerClassContainer.dL12_x1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     /myBiopolymerClassContainer.counter_threebody_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     *myBiopolymerClassContainer.dL13_x1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     /myBiopolymerClassContainer.counter_threebody_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     *myBiopolymerClassContainer.dL23_x1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     /myBiopolymerClassContainer.counter_threebody_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count])
     *myBiopolymerClassContainer.bias_angle[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]*
     myBiopolymerClassContainer.bias_3body_x_last[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];     
     
     myBiopolymerClassContainer.bias_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]  = 
     (myBiopolymerClassContainer.dL12_y1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     /myBiopolymerClassContainer.counter_threebody_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     *myBiopolymerClassContainer.dL13_y1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     /myBiopolymerClassContainer.counter_threebody_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     *myBiopolymerClassContainer.dL23_y1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     /myBiopolymerClassContainer.counter_threebody_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count])
     *myBiopolymerClassContainer.bias_angle[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]*
     myBiopolymerClassContainer.bias_3body_y_last[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];
     
     myBiopolymerClassContainer.bias_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]  = 
     (myBiopolymerClassContainer.dL12_z1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     /myBiopolymerClassContainer.counter_threebody_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     *myBiopolymerClassContainer.dL13_z1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     /myBiopolymerClassContainer.counter_threebody_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     *myBiopolymerClassContainer.dL23_z1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     /myBiopolymerClassContainer.counter_threebody_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count])
     *myBiopolymerClassContainer.bias_angle[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]*
     myBiopolymerClassContainer.bias_3body_z_last[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];  
      
 if(myBiopolymerClassContainer.meas_12_3body_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] > 1 &&
    myBiopolymerClassContainer.meas_12_3body_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] <= 100 &&
    myBiopolymerClassContainer.meas_12_3body_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] > 1 &&
    myBiopolymerClassContainer.meas_12_3body_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] <= 100 && 
    myBiopolymerClassContainer.meas_12_3body_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] > 1 &&
    myBiopolymerClassContainer.meas_12_3body_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] <= 100
) {
     
  f_bias[0] = myBiopolymerClassContainer.f_bias_3body_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];
  f_bias[1] = myBiopolymerClassContainer.f_bias_3body_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];  
  f_bias[2] = myBiopolymerClassContainer.f_bias_3body_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];
  
  d_tt = sqrt(pow(f_bias[0],2) + pow(f_bias[1],2) + pow(f_bias[2],2));
  
  if(d_tt > 0.0) {
   
     f_bias[0] /= d_tt;
     f_bias[1] /= d_tt;
     f_bias[2] /= d_tt;
      
  }
  
 if(myBiopolymerClassContainer.bias_3body_norm[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] > 0.0) {  
   
     sigma_x = myBiopolymerClassContainer.bias_3body_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]/myBiopolymerClassContainer.bias_3body_norm[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];
     sigma_y = myBiopolymerClassContainer.bias_3body_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]/myBiopolymerClassContainer.bias_3body_norm[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];
     sigma_z = myBiopolymerClassContainer.bias_3body_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]/myBiopolymerClassContainer.bias_3body_norm[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];      
        
   if( (int)myBiopolymerClassContainer.counter_gen[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]%(int)myBiopolymerClassContainer.tau_1_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] == 0) {   
      
      range_x = -1E+0;
      range_y = -1E+0;
      range_z = -1E+0;

      for(k=1;k<=100;k++) {

       if(sigma_x != 0.0) {     
          
          myBiopolymerClassContainer.alpha_2_x_threebody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][k] -=
          (myParameterReader.ntc_class_container.getNTC_Class(r)).alpha*exp(-pow(
          (sigma_x - range_x)/(0.1),2))*exp(myBiopolymerClassContainer.alpha_2_x_threebody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][k]/1000.0);
          
          range_x += 0.02;
    
       };

      if(sigma_y != 0.0) {  
        
          myBiopolymerClassContainer.alpha_2_y_threebody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][k] -=
          (myParameterReader.ntc_class_container.getNTC_Class(r)).alpha*exp(-pow(
          (sigma_y - range_y)/(0.1),2))*exp(myBiopolymerClassContainer.alpha_2_y_threebody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][k]/1000.0);
          
          range_y += 0.02;
    
       };
            
      if(sigma_z != 0.0) {  
        
          myBiopolymerClassContainer.alpha_2_z_threebody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][k] -=
          (myParameterReader.ntc_class_container.getNTC_Class(r)).alpha*exp(-pow(
          (sigma_z - range_z)/(0.1),2))*exp(myBiopolymerClassContainer.alpha_2_z_threebody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][k]/1000.0);
          
          range_z += 0.02;
    
       };
      };      
     };     
   };
       
      pos_x = (int)(sigma_x+ 1E+0)*50.0;
      pos_y = (int)(sigma_y+ 1E+0)*50.0;      
      pos_z = (int)(sigma_z+ 1E+0)*50.0;            
      
      double abbassin_x;
      double abbassin_to_x,abbassin_to_y,abbassin_to_z;
      
  if(pos_x > 1 && pos_x <= 100 && pos_y > 1 && pos_y <= 100 && pos_z > 1 && pos_z <= 100) {    
      
      double p_min_x,p_min_y,p_min_z;
      
      p_min_x = 2E+0*myBiopolymerClassContainer.bias_3body_norm[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];

      abbassin_to_x = 0.0;
      abbassin_to_y = 0.0;
      abbassin_to_z = 0.0;      
      cosgamma      = 0.0;
      
      if(p_min_x != 0.0) {
                                                            
      abbassin_to_x = ((f_bias[0]*(myBiopolymerClassContainer.alpha_2_x_threebody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][pos_x] 
                                                                - myBiopolymerClassContainer.alpha_2_x_threebody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][pos_x-1])*1E+2/p_min_x));                                                            

      }
      
      p_min_y = 2E+0*myBiopolymerClassContainer.bias_3body_norm[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];
      
      if(p_min_y != 0.0) {      
      
                                                                
       abbassin_to_y = ((f_bias[1]*(myBiopolymerClassContainer.alpha_2_y_threebody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][pos_y] 
                                                                - myBiopolymerClassContainer.alpha_2_y_threebody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][pos_y-1])*1E+2/p_min_y));                                                         
                                                            
      
      }
      
      p_min_z = 2E+0*myBiopolymerClassContainer.bias_3body_norm[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];
      
      if(p_min_z != 0.0) {        
      

      abbassin_to_z = ((f_bias[2]*(myBiopolymerClassContainer.alpha_2_z_threebody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][pos_z] 
                                                                - myBiopolymerClassContainer.alpha_2_z_threebody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][pos_z-1])*1E+2/p_min_z));                                                      
                                                            
      }
      
      abbassin_x = sqrt(pow(abbassin_to_x,2) + pow(abbassin_to_y,2) + pow(abbassin_to_z,2));
      
      d_tt = sqrt(pow(bodyForces[body1.getMobilizedBodyIndex()][1][0] + bodyForces[body1.getMobilizedBodyIndex()][1][1] + bodyForces[body1.getMobilizedBodyIndex()][1][2],2));      
      
    if( (abbassin_x) > 0.0 && d_tt > 0.0) {  
      
      cosgamma = (abbassin_to_x*bodyForces[body1.getMobilizedBodyIndex()][1][0] + abbassin_to_y*bodyForces[body1.getMobilizedBodyIndex()][1][1] +
                  abbassin_to_z*bodyForces[body1.getMobilizedBodyIndex()][1][2])/(abbassin_x*d_tt);
             
      
    }
      
     myBiopolymerClassContainer.l_max_tt[0] += cosgamma;
      myBiopolymerClassContainer.counter_l_max_tt[0] += 1.0;
      
      bodyForces[body1.getMobilizedBodyIndex()][1][0] += abbassin_to_x*cosgamma;
      bodyForces[body1.getMobilizedBodyIndex()][1][1] += abbassin_to_y*cosgamma;       
      bodyForces[body1.getMobilizedBodyIndex()][1][2] += abbassin_to_z*cosgamma;        
      
  if(rand()%100/100.0 < 0.5) {    
      
      abbassin_to_x = 0.0;
      abbassin_to_y = 0.0;
      abbassin_to_z = 0.0;
      cosgamma      = 0.0;
      
      if(p_min_x != 0.0) {   
        

      abbassin_to_x = -((f_bias[0]*(myBiopolymerClassContainer.alpha_2_x_threebody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][pos_x] 
                                                                - myBiopolymerClassContainer.alpha_2_x_threebody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][pos_x-1])*1E+2/p_min_x));                                                              
      
      }
      
      
      if(p_min_y != 0.0) { 
          
                                                            
       abbassin_to_y = -((f_bias[1]*(myBiopolymerClassContainer.alpha_2_y_threebody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][pos_y] 
                                                                - myBiopolymerClassContainer.alpha_2_y_threebody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][pos_y-1])*1E+2/p_min_y));                                                             

      }
      
      if(p_min_z != 0.0) {       
      
                                                            
      abbassin_to_z = -((f_bias[2]*(myBiopolymerClassContainer.alpha_2_z_threebody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][pos_z] 
                                                                - myBiopolymerClassContainer.alpha_2_z_threebody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][pos_z-1])*1E+2/p_min_z));                                                             

      }
      
      d_tt = sqrt(pow(bodyForces[body2.getMobilizedBodyIndex()][1][0] + bodyForces[body2.getMobilizedBodyIndex()][1][1] + bodyForces[body2.getMobilizedBodyIndex()][1][2],2)); 
      
      abbassin_x = sqrt(pow(abbassin_to_x,2) + pow(abbassin_to_y,2) + pow(abbassin_to_z,2));
     
    if( (abbassin_x) > 0.0 && d_tt > 0.0) {  
      
      cosgamma = (abbassin_to_x*bodyForces[body2.getMobilizedBodyIndex()][1][0] + abbassin_to_y*bodyForces[body2.getMobilizedBodyIndex()][1][1] +
                  abbassin_to_z*bodyForces[body2.getMobilizedBodyIndex()][1][2])/(abbassin_x*d_tt);
           
      
    }   
    
     myBiopolymerClassContainer.l_max_tt[0] += cosgamma;
      myBiopolymerClassContainer.counter_l_max_tt[0] += 1.0;    
    
      bodyForces[body2.getMobilizedBodyIndex()][1][0] += abbassin_to_x*cosgamma;
      bodyForces[body2.getMobilizedBodyIndex()][1][1] += abbassin_to_y*cosgamma;       
      bodyForces[body2.getMobilizedBodyIndex()][1][2] += abbassin_to_z*cosgamma;       
    
  } else {
    
      abbassin_to_x = 0.0;
      abbassin_to_y = 0.0;
      abbassin_to_z = 0.0;
      cosgamma      = 0.0;
      
      if(p_min_x != 0.0) {     
      

      abbassin_to_x = -f_bias[0]*(myBiopolymerClassContainer.alpha_2_x_threebody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][pos_x] 
                                                                - myBiopolymerClassContainer.alpha_2_x_threebody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][pos_x-1])*1E+2/p_min_x;
                                                            
      }
                  
      if(p_min_y != 0.0) { 
                  
                                                            
       abbassin_to_y = -f_bias[1]*(myBiopolymerClassContainer.alpha_2_y_threebody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][pos_y] 
                                                                - myBiopolymerClassContainer.alpha_2_y_threebody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][pos_y-1])*1E+2/p_min_y;

      }
      
      if(p_min_z != 0.0) {       
      
                                                            
       abbassin_to_z = -f_bias[2]*(myBiopolymerClassContainer.alpha_2_z_threebody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][pos_z] 
                                                                - myBiopolymerClassContainer.alpha_2_z_threebody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][pos_z-1])*1E+2/p_min_z
                                                            ;
      
     };    
     
      d_tt = sqrt(pow(bodyForces[body3.getMobilizedBodyIndex()][1][0] + bodyForces[body3.getMobilizedBodyIndex()][1][1] + bodyForces[body3.getMobilizedBodyIndex()][1][2],2));      
      
      abbassin_x = sqrt(pow(abbassin_to_x,2) + pow(abbassin_to_y,2) + pow(abbassin_to_z,2));
     
    if( (abbassin_x) > 0.0 && d_tt > 0.0) {  
      
      cosgamma = (abbassin_to_x*bodyForces[body3.getMobilizedBodyIndex()][1][0] + abbassin_to_y*bodyForces[body3.getMobilizedBodyIndex()][1][1] +
                  abbassin_to_z*bodyForces[body3.getMobilizedBodyIndex()][1][2])/(abbassin_x*d_tt);
             
      
      }; 
      
     myBiopolymerClassContainer.l_max_tt[0] += cosgamma;
      myBiopolymerClassContainer.counter_l_max_tt[0] += 1.0;      
      
      bodyForces[body3.getMobilizedBodyIndex()][1][0] += abbassin_to_x*cosgamma;
      bodyForces[body3.getMobilizedBodyIndex()][1][1] += abbassin_to_y*cosgamma;       
      bodyForces[body3.getMobilizedBodyIndex()][1][2] += abbassin_to_z*cosgamma;         
      
      
     };
    };
   };
   
  if(isnan(myBiopolymerClassContainer.bias_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]) == 0 &&
     isnan(myBiopolymerClassContainer.bias_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]) == 0 && 
     isnan(myBiopolymerClassContainer.bias_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]) == 0) {
      
   myBiopolymerClassContainer.beta_2_x_threebody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 
   myBiopolymerClassContainer.bias_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]*
   (1.0 - myBiopolymerClassContainer.h_div_w_3body_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]); /*- 
   myBiopolymerClassContainer.beta_2_x_threebody_last[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];*/
  
   myBiopolymerClassContainer.beta_2_y_threebody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 
   myBiopolymerClassContainer.bias_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]*
   (1.0 - myBiopolymerClassContainer.h_div_w_3body_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]); /*- 
   myBiopolymerClassContainer.beta_2_y_threebody_last[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];*/
  
   myBiopolymerClassContainer.beta_2_z_threebody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 
   myBiopolymerClassContainer.bias_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]*
   (1.0 - myBiopolymerClassContainer.h_div_w_3body_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]); /*- 
   myBiopolymerClassContainer.beta_2_z_threebody_last[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];*/ 
   
   bodyForces[body1.getMobilizedBodyIndex()][1][0] *= (1.0 + myBiopolymerClassContainer.beta_2_x_threebody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count])*exp(-myBiopolymerClassContainer.beta_2_x_threebody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]);
   
   bodyForces[body1.getMobilizedBodyIndex()][1][1] *= (1.0 + myBiopolymerClassContainer.beta_2_y_threebody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count])*exp(-myBiopolymerClassContainer.beta_2_y_threebody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]);   
   
   bodyForces[body1.getMobilizedBodyIndex()][1][2] *= (1.0 + myBiopolymerClassContainer.beta_2_z_threebody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count])*exp(-myBiopolymerClassContainer.beta_2_z_threebody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]);   
   
   bodyForces[body2.getMobilizedBodyIndex()][1][0] *= (1.0 +
   myBiopolymerClassContainer.beta_2_x_threebody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count])*exp(-myBiopolymerClassContainer.beta_2_x_threebody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]);
   
   bodyForces[body2.getMobilizedBodyIndex()][1][1] *= (1.0 + myBiopolymerClassContainer.beta_2_y_threebody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count])*exp(-myBiopolymerClassContainer.beta_2_y_threebody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]);   
   
   bodyForces[body2.getMobilizedBodyIndex()][1][2] *= (1.0 + myBiopolymerClassContainer.beta_2_z_threebody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count])*exp(-myBiopolymerClassContainer.beta_2_z_threebody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]);  
       
   bodyForces[body3.getMobilizedBodyIndex()][1][0] *= (1.0 + myBiopolymerClassContainer.beta_2_x_threebody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count])*exp(-myBiopolymerClassContainer.beta_2_x_threebody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]);
  
   bodyForces[body3.getMobilizedBodyIndex()][1][1] *= (1.0 + myBiopolymerClassContainer.beta_2_y_threebody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count])*exp(-myBiopolymerClassContainer.beta_2_y_threebody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]);   
   
   bodyForces[body3.getMobilizedBodyIndex()][1][2] *= (1.0 + myBiopolymerClassContainer.beta_2_z_threebody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count])*exp(-myBiopolymerClassContainer.beta_2_z_threebody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]);
   
   };
  };
//   std::cout << " Three body " << f_bias[0] << " " << f_bias[1] << " " << f_bias[2] << "\n";        
 };
};
    
};
   

for (int r=0;r<myParameterReader.ntc_class_container.numNTC_Torsions();r++) {    
    
// if((myParameterReader.ntc_class_container.getNTC_Class(r)).count != (myParameterReader.ntc_class_container.getNTC_Class(r+1)).count) {     
    
  if((myParameterReader.ntc_class_container.getNTC_Class(r)).twobody == 1) {

   for(n22 = 0; n22 <= 0; n22 ++) {   
      
      Vec3 d12;
      Vec3 f_bias;
    
      String chainId1=(myParameterReader.ntc_class_container.getNTC_Class(r)).NtC_FirstBPChain; 
      
      ResidueID residueNumber1=(myParameterReader.ntc_class_container.getNTC_Class(r)).FirstBPResidue;      
      ResidueID residueNumber2=(myParameterReader.ntc_class_container.getNTC_Class(r)).sec_res_twobody;
      
      NTC_PAR_BondRow myNTC_PAR_BondRow = myNTC_PAR_Class.myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[(myParameterReader.ntc_class_container.getNTC_Class(r)).NTC_PAR_BondRowIndex];
      
      all_refs_true = 0;
      
    if( myBiopolymerClassContainer.updBiopolymerClass(chainId1).hasAtom(residueNumber1,"C2") == true) {    
      
      body1 = myBiopolymerClassContainer.updAtomMobilizedBody(matter,chainId1,residueNumber1,"C2");
      Vec3 stationA = myBiopolymerClassContainer.getAtomLocationInMobilizedBodyFrame(chainId1,residueNumber1,"C2");
      state_1 = myBiopolymerClassContainer.calcAtomLocationInGroundFrame(state,chainId1,residueNumber1,"C2");       
      all_refs_true += 1;
    }
      
      String chainId2=(myParameterReader.ntc_class_container.getNTC_Class(r)).sec_res_chain_twobody;
     
    if( myBiopolymerClassContainer.updBiopolymerClass(chainId2).hasAtom(residueNumber2,"C2") == true) {      
      
      body2 = myBiopolymerClassContainer.updAtomMobilizedBody(matter,chainId2,residueNumber2,"C2");
      Vec3 stationB = myBiopolymerClassContainer.getAtomLocationInMobilizedBodyFrame(chainId2,residueNumber2,"C2");
      state_2 = myBiopolymerClassContainer.calcAtomLocationInGroundFrame(state,chainId2,residueNumber2,"C2");     
      all_refs_true += 1;
            
    }
      
 // 12
     
  myBiopolymerClassContainer.h_div_w_2body_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 0.0;
  myBiopolymerClassContainer.h_div_w_2body_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 0.0;  
  myBiopolymerClassContainer.h_div_w_2body_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 0.0;
  
  if(all_refs_true == 2) {   
     
     myBiopolymerClassContainer.bias_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]  = myBiopolymerClassContainer.dL12b_x1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     /myBiopolymerClassContainer.counter_twobody_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]*
     myBiopolymerClassContainer.bias_2body_x_last[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];
     myBiopolymerClassContainer.bias_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]  = myBiopolymerClassContainer.dL12b_y1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     /myBiopolymerClassContainer.counter_twobody_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]*
     myBiopolymerClassContainer.bias_2body_y_last[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];
     myBiopolymerClassContainer.bias_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]  = myBiopolymerClassContainer.dL12b_z1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     /myBiopolymerClassContainer.counter_twobody_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]*
     myBiopolymerClassContainer.bias_2body_z_last[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];  
 
 if(myBiopolymerClassContainer.meas_12_2body_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] > 1 &&
    myBiopolymerClassContainer.meas_12_2body_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] <= 100 &&
    myBiopolymerClassContainer.meas_12_2body_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] > 1 &&
    myBiopolymerClassContainer.meas_12_2body_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] <= 100 && 
    myBiopolymerClassContainer.meas_12_2body_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] > 1 &&
    myBiopolymerClassContainer.meas_12_2body_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] <= 100) {  

     
   f_bias[0] = 0.0;
   f_bias[1] = 0.0;
   f_bias[2] = 0.0;
      
   if( myBiopolymerClassContainer.hist_global_12_2body_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][myBiopolymerClassContainer.meas_12_2body_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]] != 0.0 &&
       myBiopolymerClassContainer.weight_12_2body_x[myBiopolymerClassContainer.meas_12_2body_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]] != 0.0 &&
       myBiopolymerClassContainer.delta_x12b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] != 0.0 && p_2body_x != 0.0 &&
       myBiopolymerClassContainer.hist_global_12_2body_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][myBiopolymerClassContainer.meas_12_2body_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]-1] != 0.0 &&
       myBiopolymerClassContainer.weight_12_2body_x[myBiopolymerClassContainer.meas_12_2body_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]-1] != 0.0)
       
   {
       
  myBiopolymerClassContainer.l_max_t[0]+= myBiopolymerClassContainer.hist_global_12_2body_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][myBiopolymerClassContainer.meas_12_2body_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]]/myBiopolymerClassContainer.weight_12_2body_x[myBiopolymerClassContainer.meas_12_2body_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]];       
     
   myBiopolymerClassContainer.counter_l_max_t[0]+= 1.0;   
   
       d_x = (-log(((myBiopolymerClassContainer.hist_global_12_2body_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][myBiopolymerClassContainer.meas_12_2body_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]]/
           myBiopolymerClassContainer.weight_12_2body_x[myBiopolymerClassContainer.meas_12_2body_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]]
                    ))/p_2body_x) - myBiopolymerClassContainer.hist_global_12_2body_x_last[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]);
       
       myBiopolymerClassContainer.h_div_w_2body_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = myBiopolymerClassContainer.hist_global_12_2body_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][myBiopolymerClassContainer.meas_12_2body_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]]/
           myBiopolymerClassContainer.weight_12_2body_x[myBiopolymerClassContainer.meas_12_2body_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]];
       
       if(d_x != 0.0) {
        
           myBiopolymerClassContainer.bias_2body_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] += 
           myBiopolymerClassContainer.delta_x12b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]/d_x;
           
           myBiopolymerClassContainer.f_bias_2body_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = myBiopolymerClassContainer.delta_x12b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]/d_x;     
    
       }
           
       myBiopolymerClassContainer.hist_global_12_2body_x_last[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] =
       -log(((myBiopolymerClassContainer.hist_global_12_2body_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][myBiopolymerClassContainer.meas_12_2body_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]]/
           myBiopolymerClassContainer.weight_12_2body_x[myBiopolymerClassContainer.meas_12_2body_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]]
                    ))/p_2body_x);
       
   };
   
   if( myBiopolymerClassContainer.hist_global_12_2body_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][myBiopolymerClassContainer.meas_12_2body_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]] != 0.0 &&
       myBiopolymerClassContainer.weight_12_2body_y[myBiopolymerClassContainer.meas_12_2body_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]] != 0.0 &&
       myBiopolymerClassContainer.delta_y12b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] != 0.0 && p_2body_y != 0.0 &&
       myBiopolymerClassContainer.hist_global_12_2body_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][myBiopolymerClassContainer.meas_12_2body_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]-1] != 0.0 &&
       myBiopolymerClassContainer.weight_12_2body_y[myBiopolymerClassContainer.meas_12_2body_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]-1] != 0.0)
       
   {
     
       d_y = (-log(((myBiopolymerClassContainer.hist_global_12_2body_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][myBiopolymerClassContainer.meas_12_2body_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]]/
           myBiopolymerClassContainer.weight_12_2body_y[myBiopolymerClassContainer.meas_12_2body_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]]
                    ))/p_2body_y) - myBiopolymerClassContainer.hist_global_12_2body_y_last[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]);       
       
       myBiopolymerClassContainer.h_div_w_2body_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = myBiopolymerClassContainer.hist_global_12_2body_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][myBiopolymerClassContainer.meas_12_2body_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]]/
           myBiopolymerClassContainer.weight_12_2body_y[myBiopolymerClassContainer.meas_12_2body_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]];       
       
  myBiopolymerClassContainer.l_max_t[0]+= myBiopolymerClassContainer.hist_global_12_2body_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][myBiopolymerClassContainer.meas_12_2body_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]]/myBiopolymerClassContainer.weight_12_2body_y[myBiopolymerClassContainer.meas_12_2body_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]];        
       
   myBiopolymerClassContainer.counter_l_max_t[0]+= 1.0;   
   
       if(d_y != 0.0) {
        
           myBiopolymerClassContainer.bias_2body_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] += 
           myBiopolymerClassContainer.delta_y12b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]/d_y;
           
           myBiopolymerClassContainer.f_bias_2body_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = myBiopolymerClassContainer.delta_y12b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]/d_y;     
    
       }
           
       myBiopolymerClassContainer.hist_global_12_2body_y_last[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] =
       -log(((myBiopolymerClassContainer.hist_global_12_2body_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][myBiopolymerClassContainer.meas_12_2body_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]]/
           myBiopolymerClassContainer.weight_12_2body_y[myBiopolymerClassContainer.meas_12_2body_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]]
                    ))/p_2body_y);
       
   };   
                

   if( myBiopolymerClassContainer.hist_global_12_2body_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][myBiopolymerClassContainer.meas_12_2body_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]] != 0.0 &&
       myBiopolymerClassContainer.weight_12_2body_z[myBiopolymerClassContainer.meas_12_2body_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]] != 0.0 &&
       myBiopolymerClassContainer.delta_z12b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] != 0.0 && p_2body_z != 0.0 &&
       myBiopolymerClassContainer.hist_global_12_2body_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][myBiopolymerClassContainer.meas_12_2body_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]-1] != 0.0 &&
       myBiopolymerClassContainer.weight_12_2body_z[myBiopolymerClassContainer.meas_12_2body_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]-1] != 0.0)
       
   {
       
       d_z = (-log(((myBiopolymerClassContainer.hist_global_12_2body_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][myBiopolymerClassContainer.meas_12_2body_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]]/
           myBiopolymerClassContainer.weight_12_2body_z[myBiopolymerClassContainer.meas_12_2body_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]]
                    ))/p_2body_z) - myBiopolymerClassContainer.hist_global_12_2body_z_last[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]);
       
       myBiopolymerClassContainer.h_div_w_2body_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = myBiopolymerClassContainer.hist_global_12_2body_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][myBiopolymerClassContainer.meas_12_2body_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]]/
           myBiopolymerClassContainer.weight_12_2body_z[myBiopolymerClassContainer.meas_12_2body_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]];       
       
  myBiopolymerClassContainer.l_max_t[0]+= myBiopolymerClassContainer.hist_global_12_2body_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][myBiopolymerClassContainer.meas_12_2body_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]]/myBiopolymerClassContainer.weight_12_2body_z[myBiopolymerClassContainer.meas_12_2body_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]];        
       
   myBiopolymerClassContainer.counter_l_max_t[0]+= 1.0;
   
       if(d_z != 0.0) {
        
           myBiopolymerClassContainer.bias_2body_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] += 
           myBiopolymerClassContainer.delta_z12b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]/d_z;           
           
           myBiopolymerClassContainer.f_bias_2body_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = myBiopolymerClassContainer.delta_z12b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]/d_z;     
    
       }
           
       myBiopolymerClassContainer.hist_global_12_2body_y_last[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] =
       -log(((myBiopolymerClassContainer.hist_global_12_2body_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][myBiopolymerClassContainer.meas_12_2body_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]]/
           myBiopolymerClassContainer.weight_12_2body_y[myBiopolymerClassContainer.meas_12_2body_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]]
                    ))/p_2body_y);
       
   };       
     
   myBiopolymerClassContainer.bias_2body_norm[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 0.0;
   
   d_tt =
   sqrt(pow(myBiopolymerClassContainer.bias_2body_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count],2) +
        pow(myBiopolymerClassContainer.bias_2body_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count],2) +
        pow(myBiopolymerClassContainer.bias_2body_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count],2));     
     
   if(myBiopolymerClassContainer.bias_2body_norm[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] < d_tt) {
      
      myBiopolymerClassContainer.bias_2body_norm[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = d_tt; 
      
     };
    };
   }
  }
 }
}

for (int r=0;r<myParameterReader.ntc_class_container.numNTC_Torsions();r++) {    
    
// if((myParameterReader.ntc_class_container.getNTC_Class(r)).count != (myParameterReader.ntc_class_container.getNTC_Class(r+1)).count) {     
    
  if((myParameterReader.ntc_class_container.getNTC_Class(r)).twobody == 1) {

   for(n22 = 0; n22 <= 0; n22 ++) {   
      
      Vec3 d12;
      Vec3 f_bias;
          
    
      String chainId1=(myParameterReader.ntc_class_container.getNTC_Class(r)).NtC_FirstBPChain; 
      
      ResidueID residueNumber1=(myParameterReader.ntc_class_container.getNTC_Class(r)).FirstBPResidue;      
      ResidueID residueNumber2=(myParameterReader.ntc_class_container.getNTC_Class(r)).sec_res_twobody;
      
      NTC_PAR_BondRow myNTC_PAR_BondRow = myNTC_PAR_Class.myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[(myParameterReader.ntc_class_container.getNTC_Class(r)).NTC_PAR_BondRowIndex];
      
      all_refs_true = 0;
      
    if( myBiopolymerClassContainer.updBiopolymerClass(chainId1).hasAtom(residueNumber1,"C2") == true) {      
      
      body1 = myBiopolymerClassContainer.updAtomMobilizedBody(matter,chainId1,residueNumber1,"C2");
      Vec3 stationA = myBiopolymerClassContainer.getAtomLocationInMobilizedBodyFrame(chainId1,residueNumber1,"C2");
      state_1 = myBiopolymerClassContainer.calcAtomLocationInGroundFrame(state,chainId1,residueNumber1,"C2");       
            all_refs_true += 1;
    }
      
      String chainId2=(myParameterReader.ntc_class_container.getNTC_Class(r)).sec_res_chain_twobody;
     
    if( myBiopolymerClassContainer.updBiopolymerClass(chainId2).hasAtom(residueNumber2,"C2") == true) {     
      
      body2 = myBiopolymerClassContainer.updAtomMobilizedBody(matter,chainId2,residueNumber2,"C2");
      Vec3 stationB = myBiopolymerClassContainer.getAtomLocationInMobilizedBodyFrame(chainId2,residueNumber2,"C2");
      state_2 = myBiopolymerClassContainer.calcAtomLocationInGroundFrame(state,chainId2,residueNumber2,"C2");     
            all_refs_true += 1;
    }
      
 // 12
     
   if(all_refs_true == 2) {  
     
     myBiopolymerClassContainer.bias_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]  = myBiopolymerClassContainer.dL12b_x1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     /myBiopolymerClassContainer.counter_twobody_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]*
     myBiopolymerClassContainer.bias_2body_x_last[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];
     myBiopolymerClassContainer.bias_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]  = myBiopolymerClassContainer.dL12b_y1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     /myBiopolymerClassContainer.counter_twobody_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]*
     myBiopolymerClassContainer.bias_2body_y_last[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];
     myBiopolymerClassContainer.bias_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]  = myBiopolymerClassContainer.dL12b_z1b[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]
     /myBiopolymerClassContainer.counter_twobody_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]*
     myBiopolymerClassContainer.bias_2body_z_last[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];  
 
 if(myBiopolymerClassContainer.meas_12_2body_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] > 1 &&
    myBiopolymerClassContainer.meas_12_2body_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] <= 100 &&
    myBiopolymerClassContainer.meas_12_2body_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] > 1 &&
    myBiopolymerClassContainer.meas_12_2body_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] <= 100 && 
    myBiopolymerClassContainer.meas_12_2body_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] > 1 &&
    myBiopolymerClassContainer.meas_12_2body_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] <= 100) {
 
 if(myBiopolymerClassContainer.bias_2body_norm[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] > 0.0) { 
  
     sigma_x = myBiopolymerClassContainer.bias_2body_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]/
     myBiopolymerClassContainer.bias_2body_norm[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];
     sigma_y = myBiopolymerClassContainer.bias_2body_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]/
     myBiopolymerClassContainer.bias_2body_norm[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];
     sigma_z = myBiopolymerClassContainer.bias_2body_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]/
     myBiopolymerClassContainer.bias_2body_norm[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];
     
  if( (int)myBiopolymerClassContainer.counter_gen[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]%(int)myBiopolymerClassContainer.tau_2_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] == 0) {     
     
      range_x = -1E+0;
      range_y = -1E+0;
      range_z = -1E+0;
      
      for(k=1;k<=100;k++) {
        
      if(sigma_x != 0.0) {      
          
          myBiopolymerClassContainer.alpha_2_x_twobody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][k] -=
          (myParameterReader.ntc_class_container.getNTC_Class(r)).alpha*exp(-pow(
          (sigma_x - range_x)/0.1,2))*exp(myBiopolymerClassContainer.alpha_2_x_twobody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][k]/1000.0);
                
          range_x += 0.02;
          
       };
        
          if(sigma_y != 0.0) { 
          
          myBiopolymerClassContainer.alpha_2_y_twobody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][k] -=
          (myParameterReader.ntc_class_container.getNTC_Class(r)).alpha*exp(-pow(
          (sigma_y - range_y)/0.1,2))*exp(myBiopolymerClassContainer.alpha_2_y_twobody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][k]/1000.0);
                
          range_y += 0.02;
          
        };
        
          if(sigma_z != 0.0) {      
          
          myBiopolymerClassContainer.alpha_2_z_twobody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][k] -=
          (myParameterReader.ntc_class_container.getNTC_Class(r)).alpha*exp(-pow(
          (sigma_z - range_z)/0.1,2))*exp(myBiopolymerClassContainer.alpha_2_z_twobody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][k]/1000.0);
                
          range_z += 0.02;
          
       };
      };      
     };     
   };
     
    pos_x = (int)(sigma_x+ 1E+0)*50.0;
    pos_y = (int)(sigma_y+ 1E+0)*50.0;      
    pos_z = (int)(sigma_z+ 1E+0)*50.0; 
    
    f_bias[0] = myBiopolymerClassContainer.f_bias_2body_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];
    f_bias[1] = myBiopolymerClassContainer.f_bias_2body_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];
    f_bias[2] = myBiopolymerClassContainer.f_bias_2body_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];
    
    d_tt = sqrt(pow(f_bias[0],2) + pow(f_bias[1],2) + pow(f_bias[2],2));
    
  if(d_tt > 0.0) {
    
     f_bias[0] /= d_tt;
     f_bias[1] /= d_tt;
     f_bias[2] /= d_tt;
      
  };  
    
   double abbassin_x;
   double abbassin_to_x,abbassin_to_y,abbassin_to_z;  
  
  if(pos_x > 1 && pos_x <= 100 && pos_y > 1 && pos_y <= 100 && pos_z > 1 && pos_z <= 100) {    
      
      double p_min_x,p_min_y,p_min_z;
      
      p_min_x = 2E+0*myBiopolymerClassContainer.bias_2body_norm[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];

      abbassin_to_x = 0.0;
      abbassin_to_y = 0.0;
      abbassin_to_z = 0.0;    
      cosgamma      = 0.0;
      
      if(p_min_x != 0.0) {      
      
                                                            
       abbassin_to_x = bodyForces[body1.getMobilizedBodyIndex()][1][0] -((f_bias[0]*(myBiopolymerClassContainer.alpha_2_x_twobody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][pos_x] 
                                                                - myBiopolymerClassContainer.alpha_2_x_twobody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][pos_x-1])*1E+2/p_min_x));
                                                            
      }

      p_min_y = 2E+0*myBiopolymerClassContainer.bias_2body_norm[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];
      
      if(p_min_y != 0.0) {       
      
                                                            
       abbassin_to_y = bodyForces[body1.getMobilizedBodyIndex()][1][1] -((f_bias[1]*(myBiopolymerClassContainer.alpha_2_y_twobody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][pos_y] 
                                                                - myBiopolymerClassContainer.alpha_2_y_twobody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][pos_y-1])*1E+2/p_min_y));                                                            

      }
      
      p_min_z = 2E+0*myBiopolymerClassContainer.bias_2body_norm[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];
      
      if(p_min_z != 0.0) {        
    
                                                            
       abbassin_to_z = ((f_bias[2]*(myBiopolymerClassContainer.alpha_2_z_twobody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][pos_z] 
                                                                - myBiopolymerClassContainer.alpha_2_z_twobody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][pos_z-1])*1E+2/p_min_z));                                                     

      }
      
      d_tt = sqrt(pow(bodyForces[body1.getMobilizedBodyIndex()][1][0] + bodyForces[body1.getMobilizedBodyIndex()][1][1] + bodyForces[body1.getMobilizedBodyIndex()][1][2],2));     
      
      abbassin_x = sqrt(pow(abbassin_to_x,2) + pow(abbassin_to_y,2) + pow(abbassin_to_z,2));
      
    if( (abbassin_x) > 0.0 && d_tt > 0.0) {  
      
      cosgamma = (abbassin_to_x*bodyForces[body1.getMobilizedBodyIndex()][1][0] + abbassin_to_y*bodyForces[body1.getMobilizedBodyIndex()][1][1] +
                  abbassin_to_z*bodyForces[body1.getMobilizedBodyIndex()][1][2])/(abbassin_x*d_tt);         
           
      
    }
    
     myBiopolymerClassContainer.l_max_tt[0] += cosgamma;
      myBiopolymerClassContainer.counter_l_max_tt[0] += 1.0;    
    
      bodyForces[body1.getMobilizedBodyIndex()][1][0] += abbassin_to_x*cosgamma;
      bodyForces[body1.getMobilizedBodyIndex()][1][1] += abbassin_to_y*cosgamma;       
      bodyForces[body1.getMobilizedBodyIndex()][1][2] += abbassin_to_z*cosgamma;      
      
      abbassin_to_x = 0.0;
      abbassin_to_y = 0.0;
      abbassin_to_z = 0.0;   
      cosgamma      = 0.0;
      
      if(p_min_x != 0.0) {
      
       abbassin_to_x = -((f_bias[0]*(myBiopolymerClassContainer.alpha_2_x_twobody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][pos_x] 
                                                                - myBiopolymerClassContainer.alpha_2_x_twobody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][pos_x-1])*1E+2/p_min_x));                                                            
                                                            
      }
      
      if(p_min_y != 0.0) {      

                                                            
       abbassin_to_y = -((f_bias[1]*(myBiopolymerClassContainer.alpha_2_y_twobody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][pos_y] 
                                                                - myBiopolymerClassContainer.alpha_2_y_twobody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][pos_y-1])*1E+2/p_min_y));                                                             
                                                            
      }

      if(p_min_z != 0.0) {       
      
                                                            
       abbassin_to_z = -((f_bias[2]*(myBiopolymerClassContainer.alpha_2_z_twobody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][pos_z] 
                                                                - myBiopolymerClassContainer.alpha_2_z_twobody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count][pos_z-1])*1E+2/p_min_z));                                                               

      };
      
      d_tt = sqrt(pow(bodyForces[body2.getMobilizedBodyIndex()][1][0] + bodyForces[body2.getMobilizedBodyIndex()][1][1] + bodyForces[body2.getMobilizedBodyIndex()][1][2],2));       
     
      abbassin_x = sqrt(pow(abbassin_to_x,2) + pow(abbassin_to_y,2) + pow(abbassin_to_z,2));      
      
    if( (abbassin_x) > 0.0 && d_tt > 0.0) {  
      
      cosgamma = (abbassin_to_x*bodyForces[body2.getMobilizedBodyIndex()][1][0] + abbassin_to_y*bodyForces[body2.getMobilizedBodyIndex()][1][1] +
                  abbassin_to_z*bodyForces[body2.getMobilizedBodyIndex()][1][2])/(abbassin_x*d_tt);                         
      
    }
    
     myBiopolymerClassContainer.l_max_tt[0] += cosgamma;
      myBiopolymerClassContainer.counter_l_max_tt[0] += 1.0;    
    
      bodyForces[body2.getMobilizedBodyIndex()][1][0] += abbassin_to_x*cosgamma;
      bodyForces[body2.getMobilizedBodyIndex()][1][1] += abbassin_to_y*cosgamma;       
      bodyForces[body2.getMobilizedBodyIndex()][1][2] += abbassin_to_z*cosgamma;    
      
     };
    };
    
  if(isnan(myBiopolymerClassContainer.bias_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]) == 0 &&
     isnan(myBiopolymerClassContainer.bias_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]) == 0 && 
     isnan(myBiopolymerClassContainer.bias_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]) == 0) {
           
      
   myBiopolymerClassContainer.beta_2_x_twobody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 
    myBiopolymerClassContainer.bias_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]*
   (1.0 - myBiopolymerClassContainer.h_div_w_2body_x[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]); /*- 
           myBiopolymerClassContainer.beta_2_x_twobody_last[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];*/
   
   myBiopolymerClassContainer.beta_2_y_twobody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 
      myBiopolymerClassContainer.bias_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]*
   (1.0 - myBiopolymerClassContainer.h_div_w_2body_y[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]); /*- 
           myBiopolymerClassContainer.beta_2_y_twobody_last[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]; */
      
   myBiopolymerClassContainer.beta_2_z_twobody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count] = 
      myBiopolymerClassContainer.bias_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]*
   (1.0 - myBiopolymerClassContainer.h_div_w_2body_z[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]); /*- 
           myBiopolymerClassContainer.beta_2_z_twobody_last[(myParameterReader.ntc_class_container.getNTC_Class(r)).count];  */     
   
               
   bodyForces[body1.getMobilizedBodyIndex()][1][0] *= (1.0 + myBiopolymerClassContainer.beta_2_x_twobody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count])*exp(-myBiopolymerClassContainer.beta_2_x_twobody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]);
   
   bodyForces[body1.getMobilizedBodyIndex()][1][1] *= (1.0 + myBiopolymerClassContainer.beta_2_y_twobody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count])*exp(-myBiopolymerClassContainer.beta_2_y_twobody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]);   
   
   bodyForces[body1.getMobilizedBodyIndex()][1][2] *= (1.0 + myBiopolymerClassContainer.beta_2_z_twobody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count])*exp(-myBiopolymerClassContainer.beta_2_z_twobody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]);    
  
   
   bodyForces[body2.getMobilizedBodyIndex()][1][0] *= (1.0 + myBiopolymerClassContainer.beta_2_x_twobody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count])*exp(-myBiopolymerClassContainer.beta_2_x_twobody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]);
   
   bodyForces[body2.getMobilizedBodyIndex()][1][1] *= (1.0 + myBiopolymerClassContainer.beta_2_y_twobody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count])*exp(-myBiopolymerClassContainer.beta_2_y_twobody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]);   
   
   bodyForces[body2.getMobilizedBodyIndex()][1][2] *= (1.0 + myBiopolymerClassContainer.beta_2_z_twobody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count])*exp(-myBiopolymerClassContainer.beta_2_z_twobody[(myParameterReader.ntc_class_container.getNTC_Class(r)).count]);    
   
   };
   };
  };
  };
//  if(r == (myParameterReader.ntc_class_container.numNTC_Torsions()-1)) std::cout << " D_T_12 : " << diff_vec12 << " D_T_23 : " << diff_vec23 << " D_T_13 : " << diff_vec13 << " DIFF_ANGLE : " << diff_angle_123 << "\n";
 };
 
 if(myBiopolymerClassContainer.timesteps[0]%1000 == 0 || myBiopolymerClassContainer.timesteps[0] == 1) {
 
myBiopolymerClassContainer.l_max_t[0] /= myBiopolymerClassContainer.counter_l_max_t[0];
myBiopolymerClassContainer.l_max_tt[0] /= myBiopolymerClassContainer.counter_l_max_tt[0];
 
 myBiopolymerClassContainer.flextimestep[0] = 0.001*(1.0 + ((myBiopolymerClassContainer.l_max_tt[0])));
 
 if(rand()%100/100.0 < (1.0 -myBiopolymerClassContainer.l_max_t[0])) myBiopolymerClassContainer.flextimestep[0] = -myBiopolymerClassContainer.flextimestep[0];
 
myBiopolymerClassContainer.l_max_t[0]= 0;
myBiopolymerClassContainer.l_max_tt[0] = 0;

myBiopolymerClassContainer.counter_l_max_t[0]= 0;
myBiopolymerClassContainer.counter_l_max_tt[0] = 0;

 };
 
 myBiopolymerClassContainer.timesteps[0]++;
 
};
 
 
Real NTC_Torque::calcPotentialEnergy(const State& state) const { 

        double energy = 0.0;  
        double rms    = 0.0;
        MobilizedBody body1;
        MobilizedBody body2;
        MobilizedBody body3;
        MobilizedBody body4;
        Transform transform1;
        Transform transform2;
        double forceConstant;
        double torqueConstant;
        double dutyCycle; //must be between 0 and 1.  at 1, force is applied all the time.  at 0, basically never applied.    
        double scrubberPeriod;
        double cutoffRadius;
        double PI = 5.14159265359;
        double pot_angle;
        double angle;
        double x_d1,x_d2,x_d3,x_d4;
        double y_d1,y_d2,y_d3,y_d4;
        double z_d1,z_d2,z_d3,z_d4;                   
        double d_d1_x,d_d1_y,d_d1_z;
        double d_d2_x,d_d2_y,d_d2_z;
        double d_d3_x,d_d3_y,d_d3_z;
        double cross_1_x,cross_1_y,cross_1_z;
        double cross_2_x,cross_2_y,cross_2_z;
        double d_t,d_t2;
        double cross_3_x,cross_3_y,cross_3_z;
        double direction_x,direction_y,direction_z;
        double rgsq,fg,hg,gaa,gbb,fga,hgb;
        double dfgx,dfgy,dfgz,dthx,dthy,dthz;
        double dtfx,dtfy,dtfz;
        double s_x2,s_y2,s_z2;
        double force_a1_x,force_a1_y,force_a1_z;
        double force_a2_x,force_a2_y,force_a2_z;        
        double force_a3_x,force_a3_y,force_a3_z;
        double force_a4_x,force_a4_y,force_a4_z;
        
        for (int r=0;r<myParameterReader.ntc_class_container.numNTC_Torsions();r++) 
        { 
	    if (myParameterReader.verbose) cout<<__FILE__<<":"<<__LINE__<<"  doing base pair #"<<r<<endl;	
        
            String chainId1=(myParameterReader.ntc_class_container.getNTC_Class(r)).NtC_FirstBPChain;            
            NTC_PAR_BondRow myNTC_PAR_BondRow = myNTC_PAR_Class.myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[(myParameterReader.ntc_class_container.getNTC_Class(r)).NTC_PAR_BondRowIndex];
            String basePairIsTwoTransformForce="ntcstep";
            ResidueID residueNumber1=(myParameterReader.ntc_class_container.getNTC_Class(r)).FirstBPResidue;
            ResidueID residueNumber2=(myParameterReader.ntc_class_container.getNTC_Class(r)).SecondBPResidue;
            ResidueID myResidueNumber,myResidueNumber2;

	    if ( myNTC_PAR_BondRow.bondLength[0] == 0.0 ) {
//		cout << "torsion " << r << " is real torsion" << endl;
            
             Vec3  state_1,state_2,state_3,state_4;
        
            if(stoi(myNTC_PAR_BondRow.atom_shift[0]) == 0) myResidueNumber = residueNumber1;
            if(stoi(myNTC_PAR_BondRow.atom_shift[0]) == 1) myResidueNumber = residueNumber2;
            
            body1 = myBiopolymerClassContainer.updAtomMobilizedBody(matter,chainId1,myResidueNumber,myNTC_PAR_BondRow.residue1Atom[0]);
            Vec3 stationA = myBiopolymerClassContainer.getAtomLocationInMobilizedBodyFrame(chainId1,myResidueNumber,myNTC_PAR_BondRow.residue1Atom[0]);
            state_1 = myBiopolymerClassContainer.calcAtomLocationInGroundFrame(state,chainId1,myResidueNumber,myNTC_PAR_BondRow.residue1Atom[0]); 
            
            if(stoi(myNTC_PAR_BondRow.atom_shift[1]) == 0) myResidueNumber = residueNumber1;
            if(stoi(myNTC_PAR_BondRow.atom_shift[1]) == 1) myResidueNumber = residueNumber2;
            
            body2 = myBiopolymerClassContainer.updAtomMobilizedBody(matter,chainId1,myResidueNumber,myNTC_PAR_BondRow.residue1Atom[1]);
            Vec3 stationB = myBiopolymerClassContainer.getAtomLocationInMobilizedBodyFrame(chainId1,myResidueNumber,myNTC_PAR_BondRow.residue1Atom[1]);
            state_2 = myBiopolymerClassContainer.calcAtomLocationInGroundFrame(state,chainId1,myResidueNumber,myNTC_PAR_BondRow.residue1Atom[1]); 
            
            if(stoi(myNTC_PAR_BondRow.atom_shift[2]) == 0) myResidueNumber = residueNumber1;
            if(stoi(myNTC_PAR_BondRow.atom_shift[2]) == 1) myResidueNumber = residueNumber2;
            
            body3 = myBiopolymerClassContainer.updAtomMobilizedBody(matter,chainId1,myResidueNumber,myNTC_PAR_BondRow.residue1Atom[2]);
            Vec3 stationC = myBiopolymerClassContainer.getAtomLocationInMobilizedBodyFrame(chainId1,myResidueNumber,myNTC_PAR_BondRow.residue1Atom[2]);
            state_3 = myBiopolymerClassContainer.calcAtomLocationInGroundFrame(state,chainId1,myResidueNumber,myNTC_PAR_BondRow.residue1Atom[2]); 

            if(stoi(myNTC_PAR_BondRow.atom_shift[3]) == 0) myResidueNumber = residueNumber1;
            if(stoi(myNTC_PAR_BondRow.atom_shift[3]) == 1) myResidueNumber = residueNumber2;

            body4 = myBiopolymerClassContainer.updAtomMobilizedBody(matter,chainId1,myResidueNumber,myNTC_PAR_BondRow.residue1Atom[3]);
            Vec3 stationD = myBiopolymerClassContainer.getAtomLocationInMobilizedBodyFrame(chainId1,myResidueNumber,myNTC_PAR_BondRow.residue1Atom[3]);
            state_4 = myBiopolymerClassContainer.calcAtomLocationInGroundFrame(state,chainId1,myResidueNumber,myNTC_PAR_BondRow.residue1Atom[3]);               
            
            torqueConstant = myNTC_PAR_BondRow.torqueConstant;            

            Vec3 d_d1,d_d2,d_d3;
            
            d_d1 = state_2 - state_1;
            d_d2 = state_3 - state_2;
            d_d3 = state_4 - state_3;
            
            Vec3 cross_1, cross_2;
            
            cross_1 = d_d1 % d_d2;
            cross_2 = d_d2 % d_d3;
            
            cross_1 = cross_1 / cross_1.norm();
            cross_2 = cross_2 / cross_2.norm();
        
            Vec3 cross_3;
            
            cross_3 = cross_1 % cross_2;
        
     angle = return_angle(cross_1,cross_2,cross_3,d_d2);
     
     double dist_ang = return_dist_ang(angle,myNTC_PAR_BondRow.rotationAngle);  

    energy += torqueConstant*cos((dist_ang + 180.0)/57.295779513)*(360.0/57.295779513+1.0)/(1.0+myNTC_PAR_BondRow.CONFALVALUE)/(360.0/57.295779513);//torqueConstant*pow(dist_ang/57.295779513,2);//-torqueConstant*(-cos(dist_ang/57.295779513)+(exp(-(pow(dist_ang,2)/(2.0*pow(l_param,2))))));
        
 //   cout << " NTC sampling - CHAIN ID = " << chainId1 << ", residuenumber " << myResidueNumber.ResidueNumber  << " difference-angle = "<< dist_ang << " , CONFALVALUE = " << myNTC_PAR_BondRow.CONFALVALUE << " , " << angle*57.295779513 << " = angle at time t for atoms  = " << myNTC_PAR_BondRow.residue1Atom[0] << " , " << myNTC_PAR_BondRow.residue1Atom[1] << " , " << myNTC_PAR_BondRow.residue1Atom[2] << " , " << myNTC_PAR_BondRow.residue1Atom[3] << " , "<< myNTC_PAR_BondRow.rotationAngle*57.295779513 << " = angle_0 from  input , " << "energy = " << energy << endl;
        
    rms   += sqrt(pow(dist_ang,2));    
            
// end real torsions
	}
// bonds
	else {
	
	    Vec3  state_1;
	    Vec3  state_2;
                    
            body1 = myBiopolymerClassContainer.updAtomMobilizedBody(matter,chainId1,residueNumber1,myNTC_PAR_BondRow.residue1Atom[0]);
            Vec3 stationA = myBiopolymerClassContainer.getAtomLocationInMobilizedBodyFrame(chainId1,residueNumber1,myNTC_PAR_BondRow.residue1Atom[0]);
            state_1 = myBiopolymerClassContainer.calcAtomLocationInGroundFrame(state,chainId1,residueNumber1,myNTC_PAR_BondRow.residue1Atom[0]); 
            
            
            body2 = myBiopolymerClassContainer.updAtomMobilizedBody(matter,chainId1,residueNumber2,myNTC_PAR_BondRow.residue1Atom[1]);
            Vec3 stationB = myBiopolymerClassContainer.getAtomLocationInMobilizedBodyFrame(chainId1,residueNumber2,myNTC_PAR_BondRow.residue1Atom[1]);
            state_2 = myBiopolymerClassContainer.calcAtomLocationInGroundFrame(state,chainId1,residueNumber2,myNTC_PAR_BondRow.residue1Atom[1]); 
            
            Vec3 diff = state_2 - state_1;
            
            double d = diff.norm();

	  //    cout << "NN|CC difference: " << (d - myNTC_PAR_BondRow.bondLength[0]) << ", current value: " << d << ", equilibrium value: " << myNTC_PAR_BondRow.bondLength[0] << endl;

          energy += myNTC_PAR_BondRow.springConstant[0]*pow(1.0-exp(-(2.0*myNTC_PAR_BondRow.CONFALVALUE)*(d - myNTC_PAR_BondRow.bondLength[0])),2);
          
	     };
       };
       
       cout << rms << " = RMSD Angle sum " << endl;
       
       return energy;   
          

         
    };
    bool NTC_Torque::dependsOnlyOnPositions() const  { 
        return true; 
    };  
    
    
Real NTC_Torque::return_dist_ang(double angle,double rotationAngle) const
{
       
     double ang_diff = (angle - rotationAngle)*57.295779513; // Deg
     double dist_ang = 180.0 - abs(180.0 - abs(ang_diff));
     int angle_1 = int(round(angle*57.295779513));
     int angle_2 = int(round(rotationAngle*57.295779513));

     int interval_begin = angle_2;
     int interval_end   = (interval_begin + 180) % 360;

     if (interval_end > interval_begin) {
         if (angle_1 < interval_begin || angle_1 > interval_end ) {
             dist_ang = -dist_ang;
        }
     } else {
         if (angle_1 < interval_begin && angle_1 > interval_end ) {
             dist_ang = -dist_ang;
        };
    };    
    
  return dist_ang;  
} 

Real NTC_Torque::return_angle(Vec3 cross_1,Vec3 cross_2,Vec3 cross_3,Vec3 d_d2) const 
{
   
   double angle; 
   double PI = 5.14159265359;
           
   Vec3 direction;
            
   direction[0] = cross_3[0]*d_d2[0];
   direction[1] = cross_3[1]*d_d2[1];  
   direction[2] = cross_3[2]*d_d2[2];
            
   double scalar_product = dot(cross_1,cross_2);
        
   if(scalar_product >= 1.0) scalar_product = 1.0;
   if(scalar_product < -1.0) scalar_product = -1.0;       
      
   angle = acos(scalar_product)*180.0/PI;
        
   if(direction[0] < 0.0 && direction[1] < 0.0 && direction[2] < 0.0) {
            
      angle = -angle;
            
   };          
            
   if(angle < 0.0) angle = angle + 360.0;
        
   angle /= 57.295779513;    
    
  return angle;  
}
