typedef struct ParameterReader_wrapper{
	int mmbID;
	// ParameterReader(const ParameterReader &);
	// ParameterReader & operator = (const ParameterReader &);
	//// ParameterReader & operator = (const ParameterReader &);
	//ErrorManager & _errorManager;
	// ParameterReader();
	//// ParameterReader();
	//vector<CovalentBondClass> additionalCovalentBondVector;
	//vector<IncludeIntraChainInterface> includeIntraChainInterfaceVector;
	//BasePairContainer basePairContainer;
	//map<const ChainResidueIndex, BasePairPartner,twoIndexCmp> basePairPartners;
	//// variables previously declared and initialized in Repel.h:;
	bool addAllAtomSterics;
	bool addAllHeavyAtomSterics;
	bool addBackboneOxygenForces;
	bool addProteinBackboneSterics;
	bool addRNABackboneSterics;
	bool addSelectedAtoms;
	bool addTestSpring;
	bool applyC1pSprings;
	int calcBaseBodyFramesAtEveryTimeStep;
	bool calcEnergy;
	double totalEnergy;
	double potentialEnergy;
	double kineticEnergy;
	bool checkSatisfied;
	////bool constrainRigidSegments;
	double constraintTolerance;
	bool guessCoordinates;
	double cutoffRadius;
	double cutoffAngle;
	double densityAtomFraction;
	const char * densityFileName;
	const char * electroDensityFileName;
	double densityForceConstant;
	double electroDensityForceConstant;
	////bool densityMapActivate;
	double excludedVolumeStiffness;
	////String firstResidueMobilizerType;
	int firstStage;
	double fitDefaultTolerance;
	double globalAmberImproperTorsionScaleFactor;
	double globalBondBendScaleFactor;
	double globalBondStretchScaleFactor;
	double globalBondTorsionScaleFactor;
	double globalCoulombScaleFactor;
	double globalGbsaScaleFactor;
	double globalVdwScaleFactor;
	double hardSphereStiffnessMultiplier;
	const char * inQVectorFileName;
	double initialSeparation;
	double integratorAccuracy;
	double integratorStepSize;
	const char * integratorType;
	double kbBackboneTorsionGlobalScaleFactor;
	int lastStage;
	const char * leontisWesthofInFileName;
	bool loadTinkerParameterFile;
	const char * outQVectorFileName;
	const char * magnesiumIonChainId;
	double magnesiumIonRadius;
	int matchDefaultSkipTopLevelTransform;
	bool matchHydrogenAtomLocations;
	bool matchPurineN1AtomLocations;
	bool matchProteinCarboxylOxygenLocations;
	bool matchExact;
	bool matchIdealized;
	bool matchOptimize;
	double matchingMinimizerTolerance;
	int numReportingIntervals;
	int minimize;
	int monteCarloRun;
	double monteCarloTemperature;
	double monteCarloTemperatureIncrement;
	double nastGlobalBondTorsionScaleFactor;
	double noseHooverTime;
	int numMagnesiumIons;
	const char * outMonteCarloFileName;
	const char * outTrajectoryFileName;
	////bool physicsWhereYouWantIt;
	float physicsRadius;
	bool piecewiseRigidify;
	double planarityThreshold;
	const char * potentialType;
	bool prioritize;
	bool proteinCapping;
	double excludedVolumeRadius;
	int readInQVector;
	bool readPreviousFrameFile;
	int readMagnesiumPositionsFromFile;
	bool removeRigidBodyMomentum;
	double removeMomentumPeriod;
	double reportingInterval;
	double restrainingForceConstant;
	double restrainingTorqueConstant;
	bool rigidifyFormedHelices;
	int rigidifyTermini;
	int satisfiedBasePairs;
	int unSatisfiedBasePairs;
	double scrubberPeriod;
	bool safeParameters;
	////int setChiBondAnti;
	int setChiBondMobility;
	////int setDefaultMDParameters;
	int setDefaultStructurePredictionParameters;
	int setDefaultThreadingParameters;
	bool setForceAndStericScrubber;
	bool setForceScrubber;
	bool setHelicalStacking;
	bool setInitialVelocities;
	bool setLoopBondMobility;
	bool setOverallBondMobility;
	bool setRemoveBasePairsInRigidStretch;
	bool setRepulsiveForce;
	bool setTemperature;
	double smallGroupInertiaMultiplier;
	bool stackAllHelicalResidues;
	const char * thermostatType;
	const char * tinkerParameterFileName;
	double twoTransformForceMultiplier;
	bool useFixedStepSize;
	bool useMultithreadedComputation;
	bool useOpenMMAcceleration;
	double vanderWallSphereRadius;
	double velocityRescalingInterval;
	bool verbose;
	int vmdOutput;
	bool waterDropletMake;
	double waterDropletRadius;
	double waterDropletX;
	double waterDropletY;
	double waterDropletZ;
	double waterInertiaMultiplier;
	bool weldToGround;
	double wkdpGlobalBondTorsionScaleFactor;
	bool writeCoordinates;
	bool writeDoublePrecisionTrajectories;
	bool writeFrameFile;
	bool writeLastFrameFile;
	const char * workingDirectory;
	bool detectConvergence;
	bool converged;
	int convergenceTimeout;
	double convergenceEpsilon;
	//BondMobility::Mobility helixBondMobility;
	//BondMobility::Mobility loopBondMobility;
	//BondMobility::Mobility overallBondMobility;
	//BondMobility::Mobility chiBondMobility;
	//Vector qVector;
	const char * lastFrameFileName;
	const char * previousFrameFileName;
	//LeontisWesthofClass myLeontisWesthofClass;
	int enforceParallelness;
	//// end of variables improted from Repel.h;
	const char * sequence;
	const char * proteinSequence;
	const char * coarseNucleicAcidSequence;
	////int numChains;
	int numFirstResidues;
	int numResetBases;
	int numProteinFirstResidues;
	int numProteinChains;
	int numTemperatures;
	int numGlobalCoulombScaleFactors;
	int numGlobalVdwScaleFactors;
	////int numDutyCycles;
	double temperature;
	double dutyCycle;
	int periodicallyUpdateParameters;
	int currentStage;
	int priority;
	////vector<Biopolymer> biopolymerVector;
	////vector<String> chainId;
	//vector<int> residueNumber;
	//map<const ChainResidueIndex, int,twoIndexCmp> residueNumberTwo;
	////vector<double> temperatureArray;
	////vector<int> temperaturePriority;
	////vector<double> dutyCycleArray;
	////vector<int> dutyCyclePriority;
	///*vector<double> globalCoulombScaleFactorArray;
	//vector<int> globalCoulombScaleFactorPriority;
	//vector<double> globalVdwScaleFactorArray;
	//vector<int> globalVdwScaleFactorPriority;
	//LeontisWesthofClass _leontisWesthofClass;
	//mutable map<const String,double> userVariables;
	//DensityMap myDensityMap;
	//DensityMap myElectroDensityMap;
	//MobilizerContainer mobilizerContainer;
	//PhysicsContainer physicsContainer;
	//ConstraintToGroundContainer constraintToGroundContainer;
	//DisplacementContainer displacementContainer;
	//AtomSpringContainer atomSpringContainer;
	//BiopolymerClassContainer myBiopolymerClassContainer;
	//MoleculeClassContainer moleculeClassContainer;
	//WaterDropletContainer waterDropletContainer;
	//map<const String,String> proteinSequences;
	//map<const String,String> coarseNucleicAcidSequences;
	//map<const String, int> numRigidSegments;
	//map<const String,int>::iterator firstResidueNumbersIterator;
	// //void addRingClosingBond(const String chainID, ResidueID residueID1, String atomName1,String bondCenterName1,  ResidueID residueID2, String atomName2,String bondCenterName2);
	//// //void addRingClosingBond(const String chainID, ResidueID residueID1, String atomName1,String bondCenterName1, ResidueID residueID2, String atomName2,String bondCenterName2);
	// void addC1pSprings (LeontisWesthofClass myLeontisWesthofClass);
	//// void addC1pSprings (LeontisWesthofClass myLeontisWesthofClass);
	// void applyAtomSprings (SimbodyMatterSubsystem & matter, GeneralForceSubsystem & forces);
	//// void applyAtomSprings (SimbodyMatterSubsystem & matter, GeneralForceSubsystem & forces);
	// void configureDumm( DuMMForceFieldSubsystem & dumm);
	//// void configureDumm( DuMMForceFieldSubsystem & dumm);
	// static double myAtoF(map<const String,double> myUserVariables,const char* value );
	//// //bool chainIsBiopolymer(String myChainId );
	// static bool aToBool( const String& name, const char* value );
	// static bool compareUpper( const String& param, const char* symbol );
	//vector<BasePair> baseOperationVector;
	//ContactContainer contactContainer;
	//DensityContainer densityContainer;
	//DensityContainer electroDensityContainer;
	//vector<SingleBondMobility> singleBondMobilityVector;
	//vector<BasePairPartner> basePairPartnerVector;
	////vector<IncludeAllNonBondAtomsInResidue> includeAllNonBondAtomsInResidueVector;
	//vector<AllResiduesWithin> includeAllResiduesWithinVector;
	//vector<IncludeNonBondAtomInBiopolymerStruct> includeNonBondAtomInBiopolymerVector;
	//vector <WaterDropletAboutResidueStruct> waterDropletAboutResidueVector;
	//vector<MobilizerDomainsInterface> mobilizerDomainsInterfaceVector;
	// void removeBasePairsInRigidStretch ();
	//// // void initializeDefaults ();
	// void printAllSettings (   ostream  & myOstream = std::cout, String remarkString = "") ;
	//// void printAllSettings ( ostream & myOstream = std::cout, String remarkString = "");
	// void removeNonPriorityBasePairs (int priorityLevel);
	//// void removeNonPriorityBasePairs (int priorityLevel);
	// //int getFirstResidueNumbers(const String myChainId) const ;
	//// // int getNumBasePairs() const;
	// // int getProteinFirstResidueNumbers(const String myProteinChainId) const ;
	// //int getBasePriority(int baseResidueNumber,String baseChain, String basePairingEdge) const ;
	// // int getNumBasePairs() const;
	// void updateBasePair(int index,
	// String ch1, int res1, String edge1,
	// String ch2, int res2, String edge2,
	// String orient);
	// void updateMobilizerStretch(int index,
	// String chainId,
	// int startRes,
	// int endRes,
	// String bondMobility);
	// void addAllResiduesWithin(String chainID, int resID, double radius);
	//// void updateAllResiduesWithin(int index, String chainID, int resID, double radius);
	// void updateAllResiduesWithin(int index, String chainID, int resID, double radius);
	// void deleteAllResiduesWithin(int index);
	//// void deleteIncludeAllNonBondAtomsInResidue(int index);
	// void updateIncludeAllNonBondAtomsInResidue(int index, String chainID, int resID);
	//// void updateIncludeAllNonBondAtomsInResidue(int index, String chainID, int resID);
	// void deleteIncludeAllNonBondAtomsInResidue(int index);
	// //int calcHighestPriority();
	//// //int calcHighestPriority();
	// //int calcLowestBondingResidue(const String myChainId) ;
	//// //bool chainIsMonoAtoms(String myChainId);
	// //int calcHighestBondingResidue(const String myChainId);
	// void setLeontisWesthofBondRowIndex();
	//// void setLeontisWesthofBondRowIndex();
	// void parameterStringInterpreter(const String & paramstr);
	//// void parameterStringInterpreter(const String & paramstr);
	// void parameterStringInterpreter(const ParameterStringClass & parameterStringClass,
	// const int readStage = 0,
	// const bool readAtOneStageOnly = false,
	// const bool readOnlyUntilStage = false,
	// const bool readExcept = false);
	// void initializeFromFileOnly(const char * parameterFileName = "./commands.dat" ) ;
	// void setFirstAndLastStage(const char * parameterFileName = "./commands.dat" ) ;
	// void loadSequencesFromPdb(const char * pdbFileName);
	//// void loadSequencesFromPdb(const char * pdbFileName);
	// //void printRigidSegments();
	//// //void printRigidSegments();
	// // void printBasePairs();
	//// // void printBasePairs();
	// // void printBaseAssignments();
	//// // void printBaseAssignments();
	// void postInitialize();
	//// void postInitialize();
	// void clearContainers();
	//// void clearContainers();
	// void clearBiopolymers();
	//// void clearBiopolymers();
	// void clearForces();
	//// void clearForces();
	// void clearConstraints();
	//// void clearConstraints();
	// // void initializeDefaults ();
	// void initializeDefaults(const char * leontisWesthofInFileName = "./parameters.csv");
	//// void initializeDefaults(const char * leontisWesthofInFileName = "./parameters.csv");
	// void initialize(const char * parameterFileName = "./commands.dat" );
	// //bool chainIsBiopolymer(String myChainId );
	// //bool chainIsMonoAtoms(String myChainId);
	// //int getChainIndex(String myChainId , vector<Biopolymer> & tempChain);
	//// //int getChainIndex(String myChainId , vector<Biopolymer> & tempChain);
	//MonoAtomsContainer myMonoAtomsContainer;
	// //variables for internal use only:
	// int r;
	// int ti;
	// int gcsfi;
	// int gvsfi;
	// int d;
	// char * s;
	// //int numChains ;
	// //int numProteinChains ;
	// //int prioritize ;
	// //temperature = 300;
	// //outQVectorFileName;
	// //firstStage = 1;
	// //lastStage = 0;// calcHighestPriority();
	////dutyCycle = 1;
	////priority = 0;
}ParameterReader_wrapper;

void updateParameterReader_wrapper(ParameterReader & _struct_, ParameterReader_wrapper * _wrap_){
	// ParameterReader(const ParameterReader &);
	// ParameterReader & operator = (const ParameterReader &);
//	_wrap_->&) = _struct_.&); ////// ParameterReader & operator = (const ParameterReader
//	_wrap_->_errorManager = _struct_._errorManager; ////ErrorManager &
	// ParameterReader();
//	_wrap_->ParameterReader() = _struct_.ParameterReader(); //////
//	_wrap_->additionalCovalentBondVector = _struct_.additionalCovalentBondVector; ////vector<CovalentBondClass>
//	_wrap_->includeIntraChainInterfaceVector = _struct_.includeIntraChainInterfaceVector; ////vector<IncludeIntraChainInterface>
//	_wrap_->basePairContainer = _struct_.basePairContainer; ////BasePairContainer
//	_wrap_->basePairPartners = _struct_.basePairPartners; ////map<const ChainResidueIndex, BasePairPartner,twoIndexCmp>
//	_wrap_->Repel.h: = _struct_.Repel.h:; ////// variables previously declared and initialized in
	_wrap_->addAllAtomSterics = _struct_.addAllAtomSterics; //bool
	_wrap_->addAllHeavyAtomSterics = _struct_.addAllHeavyAtomSterics; //bool
	_wrap_->addBackboneOxygenForces = _struct_.addBackboneOxygenForces; //bool
	_wrap_->addProteinBackboneSterics = _struct_.addProteinBackboneSterics; //bool
	_wrap_->addRNABackboneSterics = _struct_.addRNABackboneSterics; //bool
	_wrap_->addSelectedAtoms = _struct_.addSelectedAtoms; //bool
	_wrap_->addTestSpring = _struct_.addTestSpring; //bool
	_wrap_->applyC1pSprings = _struct_.applyC1pSprings; //bool
	_wrap_->calcBaseBodyFramesAtEveryTimeStep = _struct_.calcBaseBodyFramesAtEveryTimeStep; //int
	_wrap_->calcEnergy = _struct_.calcEnergy; //bool
	_wrap_->totalEnergy = _struct_.totalEnergy; //double
	_wrap_->potentialEnergy = _struct_.potentialEnergy; //double
	_wrap_->kineticEnergy = _struct_.kineticEnergy; //double
	_wrap_->checkSatisfied = _struct_.checkSatisfied; //bool
//	_wrap_->constrainRigidSegments = _struct_.constrainRigidSegments; //////bool
	_wrap_->constraintTolerance = _struct_.constraintTolerance; //double
	_wrap_->guessCoordinates = _struct_.guessCoordinates; //bool
	_wrap_->cutoffRadius = _struct_.cutoffRadius; //double
	_wrap_->cutoffAngle = _struct_.cutoffAngle; //double
	_wrap_->densityAtomFraction = _struct_.densityAtomFraction; //double
	_wrap_->densityFileName = strdup(_struct_.densityFileName.c_str()); //String
	_wrap_->electroDensityFileName = strdup(_struct_.electroDensityFileName.c_str()); //String
	_wrap_->densityForceConstant = _struct_.densityForceConstant; //double
	_wrap_->electroDensityForceConstant = _struct_.electroDensityForceConstant; //double
//	_wrap_->densityMapActivate = _struct_.densityMapActivate; //////bool
	_wrap_->excludedVolumeStiffness = _struct_.excludedVolumeStiffness; //double
//	_wrap_->firstResidueMobilizerType = _struct_.firstResidueMobilizerType; //////String
	_wrap_->firstStage = _struct_.firstStage; //int
	_wrap_->fitDefaultTolerance = _struct_.fitDefaultTolerance; //double
	_wrap_->globalAmberImproperTorsionScaleFactor = _struct_.globalAmberImproperTorsionScaleFactor; //double
	_wrap_->globalBondBendScaleFactor = _struct_.globalBondBendScaleFactor; //double
	_wrap_->globalBondStretchScaleFactor = _struct_.globalBondStretchScaleFactor; //double
	_wrap_->globalBondTorsionScaleFactor = _struct_.globalBondTorsionScaleFactor; //double
	_wrap_->globalCoulombScaleFactor = _struct_.globalCoulombScaleFactor; //double
	_wrap_->globalGbsaScaleFactor = _struct_.globalGbsaScaleFactor; //double
	_wrap_->globalVdwScaleFactor = _struct_.globalVdwScaleFactor; //double
	_wrap_->hardSphereStiffnessMultiplier = _struct_.hardSphereStiffnessMultiplier; //double
	_wrap_->inQVectorFileName = strdup(_struct_.inQVectorFileName.c_str()); //String
	_wrap_->initialSeparation = _struct_.initialSeparation; //double
	_wrap_->integratorAccuracy = _struct_.integratorAccuracy; //double
	_wrap_->integratorStepSize = _struct_.integratorStepSize; //double
	_wrap_->integratorType = strdup(_struct_.integratorType.c_str()); //String
	_wrap_->kbBackboneTorsionGlobalScaleFactor = _struct_.kbBackboneTorsionGlobalScaleFactor; //double
	_wrap_->lastStage = _struct_.lastStage; //int
	_wrap_->leontisWesthofInFileName = strdup(_struct_.leontisWesthofInFileName.c_str()); //String
	_wrap_->loadTinkerParameterFile = _struct_.loadTinkerParameterFile; //bool
	_wrap_->outQVectorFileName = strdup(_struct_.outQVectorFileName.c_str()); //String
	_wrap_->magnesiumIonChainId = strdup(_struct_.magnesiumIonChainId.c_str()); //String
	_wrap_->magnesiumIonRadius = _struct_.magnesiumIonRadius; //double
	_wrap_->matchDefaultSkipTopLevelTransform = _struct_.matchDefaultSkipTopLevelTransform; //int
	_wrap_->matchHydrogenAtomLocations = _struct_.matchHydrogenAtomLocations; //bool
	_wrap_->matchPurineN1AtomLocations = _struct_.matchPurineN1AtomLocations; //bool
	_wrap_->matchProteinCarboxylOxygenLocations = _struct_.matchProteinCarboxylOxygenLocations; //bool
	_wrap_->matchExact = _struct_.matchExact; //bool
	_wrap_->matchIdealized = _struct_.matchIdealized; //bool
	_wrap_->matchOptimize = _struct_.matchOptimize; //bool
	_wrap_->matchingMinimizerTolerance = _struct_.matchingMinimizerTolerance; //double
	_wrap_->numReportingIntervals = _struct_.numReportingIntervals; //int
	_wrap_->minimize = _struct_.minimize; //int
	_wrap_->monteCarloRun = _struct_.monteCarloRun; //int
	_wrap_->monteCarloTemperature = _struct_.monteCarloTemperature; //double
	_wrap_->monteCarloTemperatureIncrement = _struct_.monteCarloTemperatureIncrement; //double
	_wrap_->nastGlobalBondTorsionScaleFactor = _struct_.nastGlobalBondTorsionScaleFactor; //double
	_wrap_->noseHooverTime = _struct_.noseHooverTime; //double
	_wrap_->numMagnesiumIons = _struct_.numMagnesiumIons; //int
	_wrap_->outMonteCarloFileName = strdup(_struct_.outMonteCarloFileName.c_str()); //String
	_wrap_->outTrajectoryFileName = strdup(_struct_.outTrajectoryFileName.c_str()); //String
//	_wrap_->physicsWhereYouWantIt = _struct_.physicsWhereYouWantIt; //////bool
	_wrap_->physicsRadius = _struct_.physicsRadius; //float
	_wrap_->piecewiseRigidify = _struct_.piecewiseRigidify; //bool
	_wrap_->planarityThreshold = _struct_.planarityThreshold; //double
	_wrap_->potentialType = strdup(_struct_.potentialType.c_str()); //String
	_wrap_->prioritize = _struct_.prioritize; //bool
	_wrap_->proteinCapping = _struct_.proteinCapping; //bool
	_wrap_->excludedVolumeRadius = _struct_.excludedVolumeRadius; //double
	_wrap_->readInQVector = _struct_.readInQVector; //int
	_wrap_->readPreviousFrameFile = _struct_.readPreviousFrameFile; //bool
	_wrap_->readMagnesiumPositionsFromFile = _struct_.readMagnesiumPositionsFromFile; //int
	_wrap_->removeRigidBodyMomentum = _struct_.removeRigidBodyMomentum; //bool
	_wrap_->removeMomentumPeriod = _struct_.removeMomentumPeriod; //double
	_wrap_->reportingInterval = _struct_.reportingInterval; //double
	_wrap_->restrainingForceConstant = _struct_.restrainingForceConstant; //double
	_wrap_->restrainingTorqueConstant = _struct_.restrainingTorqueConstant; //double
	_wrap_->rigidifyFormedHelices = _struct_.rigidifyFormedHelices; //bool
	_wrap_->rigidifyTermini = _struct_.rigidifyTermini; //int
	_wrap_->satisfiedBasePairs = _struct_.satisfiedBasePairs; //int
	_wrap_->unSatisfiedBasePairs = _struct_.unSatisfiedBasePairs; //int
	_wrap_->scrubberPeriod = _struct_.scrubberPeriod; //double
	_wrap_->safeParameters = _struct_.safeParameters; //bool
//	_wrap_->setChiBondAnti = _struct_.setChiBondAnti; //////int
	_wrap_->setChiBondMobility = _struct_.setChiBondMobility; //int
//	_wrap_->setDefaultMDParameters = _struct_.setDefaultMDParameters; //////int
	_wrap_->setDefaultStructurePredictionParameters = _struct_.setDefaultStructurePredictionParameters; //int
	_wrap_->setDefaultThreadingParameters = _struct_.setDefaultThreadingParameters; //int
	_wrap_->setForceAndStericScrubber = _struct_.setForceAndStericScrubber; //bool
	_wrap_->setForceScrubber = _struct_.setForceScrubber; //bool
	_wrap_->setHelicalStacking = _struct_.setHelicalStacking; //bool
	_wrap_->setInitialVelocities = _struct_.setInitialVelocities; //bool
	_wrap_->setLoopBondMobility = _struct_.setLoopBondMobility; //bool
	_wrap_->setOverallBondMobility = _struct_.setOverallBondMobility; //bool
	_wrap_->setRemoveBasePairsInRigidStretch = _struct_.setRemoveBasePairsInRigidStretch; //bool
	_wrap_->setRepulsiveForce = _struct_.setRepulsiveForce; //bool
	_wrap_->setTemperature = _struct_.setTemperature; //bool
	_wrap_->smallGroupInertiaMultiplier = _struct_.smallGroupInertiaMultiplier; //double
	_wrap_->stackAllHelicalResidues = _struct_.stackAllHelicalResidues; //bool
	_wrap_->thermostatType = strdup(_struct_.thermostatType.c_str()); //String
	_wrap_->tinkerParameterFileName = strdup(_struct_.tinkerParameterFileName.c_str()); //String
	_wrap_->twoTransformForceMultiplier = _struct_.twoTransformForceMultiplier; //double
	_wrap_->useFixedStepSize = _struct_.useFixedStepSize; //bool
	_wrap_->useMultithreadedComputation = _struct_.useMultithreadedComputation; //bool
	_wrap_->useOpenMMAcceleration = _struct_.useOpenMMAcceleration; //bool
	_wrap_->vanderWallSphereRadius = _struct_.vanderWallSphereRadius; //double
	_wrap_->velocityRescalingInterval = _struct_.velocityRescalingInterval; //double
	_wrap_->verbose = _struct_.verbose; //bool
	_wrap_->vmdOutput = _struct_.vmdOutput; //int
	_wrap_->waterDropletMake = _struct_.waterDropletMake; //bool
	_wrap_->waterDropletRadius = _struct_.waterDropletRadius; //double
	_wrap_->waterDropletX = _struct_.waterDropletX; //double
	_wrap_->waterDropletY = _struct_.waterDropletY; //double
	_wrap_->waterDropletZ = _struct_.waterDropletZ; //double
	_wrap_->waterInertiaMultiplier = _struct_.waterInertiaMultiplier; //double
	_wrap_->weldToGround = _struct_.weldToGround; //bool
	_wrap_->wkdpGlobalBondTorsionScaleFactor = _struct_.wkdpGlobalBondTorsionScaleFactor; //double
	_wrap_->writeCoordinates = _struct_.writeCoordinates; //bool
	_wrap_->writeDoublePrecisionTrajectories = _struct_.writeDoublePrecisionTrajectories; //bool
	_wrap_->writeFrameFile = _struct_.writeFrameFile; //bool
	_wrap_->writeLastFrameFile = _struct_.writeLastFrameFile; //bool
	_wrap_->workingDirectory = strdup(_struct_.workingDirectory.c_str()); //String
	_wrap_->detectConvergence = _struct_.detectConvergence; //bool
	_wrap_->converged = _struct_.converged; //bool
	_wrap_->convergenceTimeout = _struct_.convergenceTimeout; //int
	_wrap_->convergenceEpsilon = _struct_.convergenceEpsilon; //double
//	_wrap_->helixBondMobility = _struct_.helixBondMobility; ////BondMobility::Mobility
//	_wrap_->loopBondMobility = _struct_.loopBondMobility; ////BondMobility::Mobility
//	_wrap_->overallBondMobility = _struct_.overallBondMobility; ////BondMobility::Mobility
//	_wrap_->chiBondMobility = _struct_.chiBondMobility; ////BondMobility::Mobility
//	_wrap_->qVector = _struct_.qVector; ////Vector
	_wrap_->lastFrameFileName = strdup(_struct_.lastFrameFileName.c_str()); //String
	_wrap_->previousFrameFileName = strdup(_struct_.previousFrameFileName.c_str()); //String
//	_wrap_->myLeontisWesthofClass = _struct_.myLeontisWesthofClass; ////LeontisWesthofClass
	_wrap_->enforceParallelness = _struct_.enforceParallelness; //int
//	_wrap_->Repel.h = _struct_.Repel.h; ////// end of variables improted from
	_wrap_->sequence = strdup(_struct_.sequence.c_str()); //String
	_wrap_->proteinSequence = strdup(_struct_.proteinSequence.c_str()); //String
	_wrap_->coarseNucleicAcidSequence = strdup(_struct_.coarseNucleicAcidSequence.c_str()); //String
//	_wrap_->numChains = _struct_.numChains; //////int
	_wrap_->numFirstResidues = _struct_.numFirstResidues; //int
	_wrap_->numResetBases = _struct_.numResetBases; //int
	_wrap_->numProteinFirstResidues = _struct_.numProteinFirstResidues; //int
	_wrap_->numProteinChains = _struct_.numProteinChains; //int
	_wrap_->numTemperatures = _struct_.numTemperatures; //int
	_wrap_->numGlobalCoulombScaleFactors = _struct_.numGlobalCoulombScaleFactors; //int
	_wrap_->numGlobalVdwScaleFactors = _struct_.numGlobalVdwScaleFactors; //int
//	_wrap_->numDutyCycles = _struct_.numDutyCycles; //////int
	_wrap_->temperature = _struct_.temperature; //double
	_wrap_->dutyCycle = _struct_.dutyCycle; //double
	_wrap_->periodicallyUpdateParameters = _struct_.periodicallyUpdateParameters; //int
	_wrap_->currentStage = _struct_.currentStage; //int
	_wrap_->priority = _struct_.priority; //int
//	_wrap_->biopolymerVector = _struct_.biopolymerVector; //////vector<Biopolymer>
//	_wrap_->chainId = _struct_.chainId; //////vector<String>
//	_wrap_->residueNumber = _struct_.residueNumber; ////vector<int>
//	_wrap_->residueNumberTwo = _struct_.residueNumberTwo; ////map<const ChainResidueIndex, int,twoIndexCmp>
//	_wrap_->temperatureArray = _struct_.temperatureArray; //////vector<double>
//	_wrap_->temperaturePriority = _struct_.temperaturePriority; //////vector<int>
//	_wrap_->dutyCycleArray = _struct_.dutyCycleArray; //////vector<double>
//	_wrap_->dutyCyclePriority = _struct_.dutyCyclePriority; //////vector<int>
//	_wrap_->globalCoulombScaleFactorArray = _struct_.globalCoulombScaleFactorArray; /////*vector<double>
//	_wrap_->globalCoulombScaleFactorPriority = _struct_.globalCoulombScaleFactorPriority; ////vector<int>
//	_wrap_->globalVdwScaleFactorArray = _struct_.globalVdwScaleFactorArray; ////vector<double>
//	_wrap_->globalVdwScaleFactorPriority = _struct_.globalVdwScaleFactorPriority; ////vector<int>
//	_wrap_->_leontisWesthofClass = _struct_._leontisWesthofClass; ////LeontisWesthofClass
//	_wrap_->userVariables = _struct_.userVariables; ////mutable map<const String,double>
//	_wrap_->myDensityMap = _struct_.myDensityMap; ////DensityMap
//	_wrap_->myElectroDensityMap = _struct_.myElectroDensityMap; ////DensityMap
//	_wrap_->mobilizerContainer = _struct_.mobilizerContainer; ////MobilizerContainer
//	_wrap_->physicsContainer = _struct_.physicsContainer; ////PhysicsContainer
//	_wrap_->constraintToGroundContainer = _struct_.constraintToGroundContainer; ////ConstraintToGroundContainer
//	_wrap_->displacementContainer = _struct_.displacementContainer; ////DisplacementContainer
//	_wrap_->atomSpringContainer = _struct_.atomSpringContainer; ////AtomSpringContainer
//	_wrap_->myBiopolymerClassContainer = _struct_.myBiopolymerClassContainer; ////BiopolymerClassContainer
//	_wrap_->moleculeClassContainer = _struct_.moleculeClassContainer; ////MoleculeClassContainer
//	_wrap_->waterDropletContainer = _struct_.waterDropletContainer; ////WaterDropletContainer
//	_wrap_->proteinSequences = _struct_.proteinSequences; ////map<const String,String>
//	_wrap_->coarseNucleicAcidSequences = _struct_.coarseNucleicAcidSequences; ////map<const String,String>
//	_wrap_->numRigidSegments = _struct_.numRigidSegments; ////map<const String, int>
//	_wrap_->firstResidueNumbersIterator = _struct_.firstResidueNumbersIterator; ////map<const String,int>::iterator
	// //void addRingClosingBond(const String chainID, ResidueID residueID1, String atomName1,String bondCenterName1,  ResidueID residueID2, String atomName2,String bondCenterName2);
//	_wrap_->bondCenterName2) = _struct_.bondCenterName2); ////// //void addRingClosingBond(const String chainID, ResidueID residueID1, String atomName1,String bondCenterName1, ResidueID residueID2, String atomName2,String
	// void addC1pSprings (LeontisWesthofClass myLeontisWesthofClass);
//	_wrap_->myLeontisWesthofClass) = _struct_.myLeontisWesthofClass); ////// void addC1pSprings (LeontisWesthofClass
	// void applyAtomSprings (SimbodyMatterSubsystem & matter, GeneralForceSubsystem & forces);
//	_wrap_->forces) = _struct_.forces); ////// void applyAtomSprings (SimbodyMatterSubsystem & matter, GeneralForceSubsystem &
	// void configureDumm( DuMMForceFieldSubsystem & dumm);
//	_wrap_->dumm) = _struct_.dumm); ////// void configureDumm( DuMMForceFieldSubsystem &
	// static double myAtoF(map<const String,double> myUserVariables,const char* value );
//	_wrap_->) = _struct_.); ////// //bool chainIsBiopolymer(String myChainId
	// static bool aToBool( const String& name, const char* value );
	// static bool compareUpper( const String& param, const char* symbol );
//	_wrap_->baseOperationVector = _struct_.baseOperationVector; ////vector<BasePair>
//	_wrap_->contactContainer = _struct_.contactContainer; ////ContactContainer
//	_wrap_->densityContainer = _struct_.densityContainer; ////DensityContainer
//	_wrap_->electroDensityContainer = _struct_.electroDensityContainer; ////DensityContainer
//	_wrap_->singleBondMobilityVector = _struct_.singleBondMobilityVector; ////vector<SingleBondMobility>
//	_wrap_->basePairPartnerVector = _struct_.basePairPartnerVector; ////vector<BasePairPartner>
//	_wrap_->includeAllNonBondAtomsInResidueVector = _struct_.includeAllNonBondAtomsInResidueVector; //////vector<IncludeAllNonBondAtomsInResidue>
//	_wrap_->includeAllResiduesWithinVector = _struct_.includeAllResiduesWithinVector; ////vector<AllResiduesWithin>
//	_wrap_->includeNonBondAtomInBiopolymerVector = _struct_.includeNonBondAtomInBiopolymerVector; ////vector<IncludeNonBondAtomInBiopolymerStruct>
//	_wrap_->waterDropletAboutResidueVector = _struct_.waterDropletAboutResidueVector; ////vector <WaterDropletAboutResidueStruct>
//	_wrap_->mobilizerDomainsInterfaceVector = _struct_.mobilizerDomainsInterfaceVector; ////vector<MobilizerDomainsInterface>
	// void removeBasePairsInRigidStretch ();
//	_wrap_->() = _struct_.(); ////// // void initializeDefaults
	// void printAllSettings (   ostream  & myOstream = std::cout, String remarkString = "") ;
//	_wrap_->"") = _struct_.""); ////// void printAllSettings ( ostream & myOstream = std::cout, String remarkString =
	// void removeNonPriorityBasePairs (int priorityLevel);
//	_wrap_->priorityLevel) = _struct_.priorityLevel); ////// void removeNonPriorityBasePairs (int
	// //int getFirstResidueNumbers(const String myChainId) const ;
//	_wrap_->const = _struct_.const; ////// // int getNumBasePairs()
	// // int getProteinFirstResidueNumbers(const String myProteinChainId) const ;
	// //int getBasePriority(int baseResidueNumber,String baseChain, String basePairingEdge) const ;
	// // int getNumBasePairs() const;
	// void updateBasePair(int index,
	// String ch1, int res1, String edge1,
	// String ch2, int res2, String edge2,
	// String orient);
	// void updateMobilizerStretch(int index,
	// String chainId,
	// int startRes,
	// int endRes,
	// String bondMobility);
	// void addAllResiduesWithin(String chainID, int resID, double radius);
//	_wrap_->radius) = _struct_.radius); ////// void updateAllResiduesWithin(int index, String chainID, int resID, double
	// void updateAllResiduesWithin(int index, String chainID, int resID, double radius);
	// void deleteAllResiduesWithin(int index);
//	_wrap_->index) = _struct_.index); ////// void deleteIncludeAllNonBondAtomsInResidue(int
	// void updateIncludeAllNonBondAtomsInResidue(int index, String chainID, int resID);
//	_wrap_->resID) = _struct_.resID); ////// void updateIncludeAllNonBondAtomsInResidue(int index, String chainID, int
	// void deleteIncludeAllNonBondAtomsInResidue(int index);
	// //int calcHighestPriority();
//	_wrap_->calcHighestPriority() = _struct_.calcHighestPriority(); ////// //int
	// //int calcLowestBondingResidue(const String myChainId) ;
//	_wrap_->myChainId) = _struct_.myChainId); ////// //bool chainIsMonoAtoms(String
	// //int calcHighestBondingResidue(const String myChainId);
	// void setLeontisWesthofBondRowIndex();
//	_wrap_->setLeontisWesthofBondRowIndex() = _struct_.setLeontisWesthofBondRowIndex(); ////// void
	// void parameterStringInterpreter(const String & paramstr);
//	_wrap_->paramstr) = _struct_.paramstr); ////// void parameterStringInterpreter(const String &
	// void parameterStringInterpreter(const ParameterStringClass & parameterStringClass,
	// const int readStage = 0,
	// const bool readAtOneStageOnly = false,
	// const bool readOnlyUntilStage = false,
	// const bool readExcept = false);
	// void initializeFromFileOnly(const char * parameterFileName = "./commands.dat" ) ;
	// void setFirstAndLastStage(const char * parameterFileName = "./commands.dat" ) ;
	// void loadSequencesFromPdb(const char * pdbFileName);
//	_wrap_->pdbFileName) = _struct_.pdbFileName); ////// void loadSequencesFromPdb(const char *
	// //void printRigidSegments();
//	_wrap_->printRigidSegments() = _struct_.printRigidSegments(); ////// //void
	// // void printBasePairs();
//	_wrap_->printBasePairs() = _struct_.printBasePairs(); ////// // void
	// // void printBaseAssignments();
//	_wrap_->printBaseAssignments() = _struct_.printBaseAssignments(); ////// // void
	// void postInitialize();
//	_wrap_->postInitialize() = _struct_.postInitialize(); ////// void
	// void clearContainers();
//	_wrap_->clearContainers() = _struct_.clearContainers(); ////// void
	// void clearBiopolymers();
//	_wrap_->clearBiopolymers() = _struct_.clearBiopolymers(); ////// void
	// void clearForces();
//	_wrap_->clearForces() = _struct_.clearForces(); ////// void
	// void clearConstraints();
//	_wrap_->clearConstraints() = _struct_.clearConstraints(); ////// void
	// // void initializeDefaults ();
	// void initializeDefaults(const char * leontisWesthofInFileName = "./parameters.csv");
//	_wrap_->"./parameters.csv") = _struct_."./parameters.csv"); ////// void initializeDefaults(const char * leontisWesthofInFileName =
	// void initialize(const char * parameterFileName = "./commands.dat" );
	// //bool chainIsBiopolymer(String myChainId );
	// //bool chainIsMonoAtoms(String myChainId);
	// //int getChainIndex(String myChainId , vector<Biopolymer> & tempChain);
//	_wrap_->tempChain) = _struct_.tempChain); ////// //int getChainIndex(String myChainId , vector<Biopolymer> &
//	_wrap_->myMonoAtomsContainer = _struct_.myMonoAtomsContainer; ////MonoAtomsContainer
	// //variables for internal use only:
	// int r;
	// int ti;
	// int gcsfi;
	// int gvsfi;
	// int d;
	// char * s;
	// //int numChains ;
	// //int numProteinChains ;
	// //int prioritize ;
	// //temperature = 300;
	// //outQVectorFileName;
	// //firstStage = 1;
	// //lastStage = 0;// calcHighestPriority();
//	_wrap_->1 = _struct_.1; //////dutyCycle =
//	_wrap_->0 = _struct_.0; //////priority =

}
void updateParameterReader(ParameterReader_wrapper * _wrap_, ParameterReader & _struct_){
	// ParameterReader(const ParameterReader &);
	// ParameterReader & operator = (const ParameterReader &);
//	_struct_.&) = _wrap_->&); ////// ParameterReader & operator = (const ParameterReader
//	_struct_._errorManager = _wrap_->_errorManager; ////ErrorManager &
	// ParameterReader();
//	_struct_.ParameterReader() = _wrap_->ParameterReader(); //////
//	_struct_.additionalCovalentBondVector = _wrap_->additionalCovalentBondVector; ////vector<CovalentBondClass>
//	_struct_.includeIntraChainInterfaceVector = _wrap_->includeIntraChainInterfaceVector; ////vector<IncludeIntraChainInterface>
//	_struct_.basePairContainer = _wrap_->basePairContainer; ////BasePairContainer
//	_struct_.basePairPartners = _wrap_->basePairPartners; ////map<const ChainResidueIndex, BasePairPartner,twoIndexCmp>
//	_struct_.Repel.h: = _wrap_->Repel.h:; ////// variables previously declared and initialized in
	_struct_.addAllAtomSterics = _wrap_->addAllAtomSterics; //bool
	_struct_.addAllHeavyAtomSterics = _wrap_->addAllHeavyAtomSterics; //bool
	_struct_.addBackboneOxygenForces = _wrap_->addBackboneOxygenForces; //bool
	_struct_.addProteinBackboneSterics = _wrap_->addProteinBackboneSterics; //bool
	_struct_.addRNABackboneSterics = _wrap_->addRNABackboneSterics; //bool
	_struct_.addSelectedAtoms = _wrap_->addSelectedAtoms; //bool
	_struct_.addTestSpring = _wrap_->addTestSpring; //bool
	_struct_.applyC1pSprings = _wrap_->applyC1pSprings; //bool
	_struct_.calcBaseBodyFramesAtEveryTimeStep = _wrap_->calcBaseBodyFramesAtEveryTimeStep; //int
	_struct_.calcEnergy = _wrap_->calcEnergy; //bool
	_struct_.totalEnergy = _wrap_->totalEnergy; //double
	_struct_.potentialEnergy = _wrap_->potentialEnergy; //double
	_struct_.kineticEnergy = _wrap_->kineticEnergy; //double
	_struct_.checkSatisfied = _wrap_->checkSatisfied; //bool
//	_struct_.constrainRigidSegments = _wrap_->constrainRigidSegments; //////bool
	_struct_.constraintTolerance = _wrap_->constraintTolerance; //double
	_struct_.guessCoordinates = _wrap_->guessCoordinates; //bool
	_struct_.cutoffRadius = _wrap_->cutoffRadius; //double
	_struct_.cutoffAngle = _wrap_->cutoffAngle; //double
	_struct_.densityAtomFraction = _wrap_->densityAtomFraction; //double
	_struct_.densityFileName = String(_wrap_->densityFileName); //String
	_struct_.electroDensityFileName = String(_wrap_->electroDensityFileName); //String
	_struct_.densityForceConstant = _wrap_->densityForceConstant; //double
	_struct_.electroDensityForceConstant = _wrap_->electroDensityForceConstant; //double
//	_struct_.densityMapActivate = _wrap_->densityMapActivate; //////bool
	_struct_.excludedVolumeStiffness = _wrap_->excludedVolumeStiffness; //double
//	_struct_.firstResidueMobilizerType = _wrap_->firstResidueMobilizerType; //////String
	_struct_.firstStage = _wrap_->firstStage; //int
	_struct_.fitDefaultTolerance = _wrap_->fitDefaultTolerance; //double
	_struct_.globalAmberImproperTorsionScaleFactor = _wrap_->globalAmberImproperTorsionScaleFactor; //double
	_struct_.globalBondBendScaleFactor = _wrap_->globalBondBendScaleFactor; //double
	_struct_.globalBondStretchScaleFactor = _wrap_->globalBondStretchScaleFactor; //double
	_struct_.globalBondTorsionScaleFactor = _wrap_->globalBondTorsionScaleFactor; //double
	_struct_.globalCoulombScaleFactor = _wrap_->globalCoulombScaleFactor; //double
	_struct_.globalGbsaScaleFactor = _wrap_->globalGbsaScaleFactor; //double
	_struct_.globalVdwScaleFactor = _wrap_->globalVdwScaleFactor; //double
	_struct_.hardSphereStiffnessMultiplier = _wrap_->hardSphereStiffnessMultiplier; //double
	_struct_.inQVectorFileName = String(_wrap_->inQVectorFileName); //String
	_struct_.initialSeparation = _wrap_->initialSeparation; //double
	_struct_.integratorAccuracy = _wrap_->integratorAccuracy; //double
	_struct_.integratorStepSize = _wrap_->integratorStepSize; //double
	_struct_.integratorType = String(_wrap_->integratorType); //String
	_struct_.kbBackboneTorsionGlobalScaleFactor = _wrap_->kbBackboneTorsionGlobalScaleFactor; //double
	_struct_.lastStage = _wrap_->lastStage; //int
	_struct_.leontisWesthofInFileName = String(_wrap_->leontisWesthofInFileName); //String
	_struct_.loadTinkerParameterFile = _wrap_->loadTinkerParameterFile; //bool
	_struct_.outQVectorFileName = String(_wrap_->outQVectorFileName); //String
	_struct_.magnesiumIonChainId = String(_wrap_->magnesiumIonChainId); //String
	_struct_.magnesiumIonRadius = _wrap_->magnesiumIonRadius; //double
	_struct_.matchDefaultSkipTopLevelTransform = _wrap_->matchDefaultSkipTopLevelTransform; //int
	_struct_.matchHydrogenAtomLocations = _wrap_->matchHydrogenAtomLocations; //bool
	_struct_.matchPurineN1AtomLocations = _wrap_->matchPurineN1AtomLocations; //bool
	_struct_.matchProteinCarboxylOxygenLocations = _wrap_->matchProteinCarboxylOxygenLocations; //bool
	_struct_.matchExact = _wrap_->matchExact; //bool
	_struct_.matchIdealized = _wrap_->matchIdealized; //bool
	_struct_.matchOptimize = _wrap_->matchOptimize; //bool
	_struct_.matchingMinimizerTolerance = _wrap_->matchingMinimizerTolerance; //double
	_struct_.numReportingIntervals = _wrap_->numReportingIntervals; //int
	_struct_.minimize = _wrap_->minimize; //int
	_struct_.monteCarloRun = _wrap_->monteCarloRun; //int
	_struct_.monteCarloTemperature = _wrap_->monteCarloTemperature; //double
	_struct_.monteCarloTemperatureIncrement = _wrap_->monteCarloTemperatureIncrement; //double
	_struct_.nastGlobalBondTorsionScaleFactor = _wrap_->nastGlobalBondTorsionScaleFactor; //double
	_struct_.noseHooverTime = _wrap_->noseHooverTime; //double
	_struct_.numMagnesiumIons = _wrap_->numMagnesiumIons; //int
	_struct_.outMonteCarloFileName = String(_wrap_->outMonteCarloFileName); //String
	_struct_.outTrajectoryFileName = String(_wrap_->outTrajectoryFileName); //String
//	_struct_.physicsWhereYouWantIt = _wrap_->physicsWhereYouWantIt; //////bool
	_struct_.physicsRadius = _wrap_->physicsRadius; //float
	_struct_.piecewiseRigidify = _wrap_->piecewiseRigidify; //bool
	_struct_.planarityThreshold = _wrap_->planarityThreshold; //double
	_struct_.potentialType = String(_wrap_->potentialType); //String
	_struct_.prioritize = _wrap_->prioritize; //bool
	_struct_.proteinCapping = _wrap_->proteinCapping; //bool
	_struct_.excludedVolumeRadius = _wrap_->excludedVolumeRadius; //double
	_struct_.readInQVector = _wrap_->readInQVector; //int
	_struct_.readPreviousFrameFile = _wrap_->readPreviousFrameFile; //bool
	_struct_.readMagnesiumPositionsFromFile = _wrap_->readMagnesiumPositionsFromFile; //int
	_struct_.removeRigidBodyMomentum = _wrap_->removeRigidBodyMomentum; //bool
	_struct_.removeMomentumPeriod = _wrap_->removeMomentumPeriod; //double
	_struct_.reportingInterval = _wrap_->reportingInterval; //double
	_struct_.restrainingForceConstant = _wrap_->restrainingForceConstant; //double
	_struct_.restrainingTorqueConstant = _wrap_->restrainingTorqueConstant; //double
	_struct_.rigidifyFormedHelices = _wrap_->rigidifyFormedHelices; //bool
	_struct_.rigidifyTermini = _wrap_->rigidifyTermini; //int
	_struct_.satisfiedBasePairs = _wrap_->satisfiedBasePairs; //int
	_struct_.unSatisfiedBasePairs = _wrap_->unSatisfiedBasePairs; //int
	_struct_.scrubberPeriod = _wrap_->scrubberPeriod; //double
	_struct_.safeParameters = _wrap_->safeParameters; //bool
//	_struct_.setChiBondAnti = _wrap_->setChiBondAnti; //////int
	_struct_.setChiBondMobility = _wrap_->setChiBondMobility; //int
//	_struct_.setDefaultMDParameters = _wrap_->setDefaultMDParameters; //////int
	_struct_.setDefaultStructurePredictionParameters = _wrap_->setDefaultStructurePredictionParameters; //int
	_struct_.setDefaultThreadingParameters = _wrap_->setDefaultThreadingParameters; //int
	_struct_.setForceAndStericScrubber = _wrap_->setForceAndStericScrubber; //bool
	_struct_.setForceScrubber = _wrap_->setForceScrubber; //bool
	_struct_.setHelicalStacking = _wrap_->setHelicalStacking; //bool
	_struct_.setInitialVelocities = _wrap_->setInitialVelocities; //bool
	_struct_.setLoopBondMobility = _wrap_->setLoopBondMobility; //bool
	_struct_.setOverallBondMobility = _wrap_->setOverallBondMobility; //bool
	_struct_.setRemoveBasePairsInRigidStretch = _wrap_->setRemoveBasePairsInRigidStretch; //bool
	_struct_.setRepulsiveForce = _wrap_->setRepulsiveForce; //bool
	_struct_.setTemperature = _wrap_->setTemperature; //bool
	_struct_.smallGroupInertiaMultiplier = _wrap_->smallGroupInertiaMultiplier; //double
	_struct_.stackAllHelicalResidues = _wrap_->stackAllHelicalResidues; //bool
	_struct_.thermostatType = String(_wrap_->thermostatType); //String
	_struct_.tinkerParameterFileName = String(_wrap_->tinkerParameterFileName); //String
	_struct_.twoTransformForceMultiplier = _wrap_->twoTransformForceMultiplier; //double
	_struct_.useFixedStepSize = _wrap_->useFixedStepSize; //bool
	_struct_.useMultithreadedComputation = _wrap_->useMultithreadedComputation; //bool
	_struct_.useOpenMMAcceleration = _wrap_->useOpenMMAcceleration; //bool
	_struct_.vanderWallSphereRadius = _wrap_->vanderWallSphereRadius; //double
	_struct_.velocityRescalingInterval = _wrap_->velocityRescalingInterval; //double
	_struct_.verbose = _wrap_->verbose; //bool
	_struct_.vmdOutput = _wrap_->vmdOutput; //int
	_struct_.waterDropletMake = _wrap_->waterDropletMake; //bool
	_struct_.waterDropletRadius = _wrap_->waterDropletRadius; //double
	_struct_.waterDropletX = _wrap_->waterDropletX; //double
	_struct_.waterDropletY = _wrap_->waterDropletY; //double
	_struct_.waterDropletZ = _wrap_->waterDropletZ; //double
	_struct_.waterInertiaMultiplier = _wrap_->waterInertiaMultiplier; //double
	_struct_.weldToGround = _wrap_->weldToGround; //bool
	_struct_.wkdpGlobalBondTorsionScaleFactor = _wrap_->wkdpGlobalBondTorsionScaleFactor; //double
	_struct_.writeCoordinates = _wrap_->writeCoordinates; //bool
	_struct_.writeDoublePrecisionTrajectories = _wrap_->writeDoublePrecisionTrajectories; //bool
	_struct_.writeFrameFile = _wrap_->writeFrameFile; //bool
	_struct_.writeLastFrameFile = _wrap_->writeLastFrameFile; //bool
	_struct_.workingDirectory = String(_wrap_->workingDirectory); //String
	_struct_.detectConvergence = _wrap_->detectConvergence; //bool
	_struct_.converged = _wrap_->converged; //bool
	_struct_.convergenceTimeout = _wrap_->convergenceTimeout; //int
	_struct_.convergenceEpsilon = _wrap_->convergenceEpsilon; //double
//	_struct_.helixBondMobility = _wrap_->helixBondMobility; ////BondMobility::Mobility
//	_struct_.loopBondMobility = _wrap_->loopBondMobility; ////BondMobility::Mobility
//	_struct_.overallBondMobility = _wrap_->overallBondMobility; ////BondMobility::Mobility
//	_struct_.chiBondMobility = _wrap_->chiBondMobility; ////BondMobility::Mobility
//	_struct_.qVector = _wrap_->qVector; ////Vector
	_struct_.lastFrameFileName = String(_wrap_->lastFrameFileName); //String
	_struct_.previousFrameFileName = String(_wrap_->previousFrameFileName); //String
//	_struct_.myLeontisWesthofClass = _wrap_->myLeontisWesthofClass; ////LeontisWesthofClass
	_struct_.enforceParallelness = _wrap_->enforceParallelness; //int
//	_struct_.Repel.h = _wrap_->Repel.h; ////// end of variables improted from
	_struct_.sequence = String(_wrap_->sequence); //String
	_struct_.proteinSequence = String(_wrap_->proteinSequence); //String
	_struct_.coarseNucleicAcidSequence = String(_wrap_->coarseNucleicAcidSequence); //String
//	_struct_.numChains = _wrap_->numChains; //////int
	_struct_.numFirstResidues = _wrap_->numFirstResidues; //int
	_struct_.numResetBases = _wrap_->numResetBases; //int
	_struct_.numProteinFirstResidues = _wrap_->numProteinFirstResidues; //int
	_struct_.numProteinChains = _wrap_->numProteinChains; //int
	_struct_.numTemperatures = _wrap_->numTemperatures; //int
	_struct_.numGlobalCoulombScaleFactors = _wrap_->numGlobalCoulombScaleFactors; //int
	_struct_.numGlobalVdwScaleFactors = _wrap_->numGlobalVdwScaleFactors; //int
//	_struct_.numDutyCycles = _wrap_->numDutyCycles; //////int
	_struct_.temperature = _wrap_->temperature; //double
	_struct_.dutyCycle = _wrap_->dutyCycle; //double
	_struct_.periodicallyUpdateParameters = _wrap_->periodicallyUpdateParameters; //int
	_struct_.currentStage = _wrap_->currentStage; //int
	_struct_.priority = _wrap_->priority; //int
//	_struct_.biopolymerVector = _wrap_->biopolymerVector; //////vector<Biopolymer>
//	_struct_.chainId = _wrap_->chainId; //////vector<String>
//	_struct_.residueNumber = _wrap_->residueNumber; ////vector<int>
//	_struct_.residueNumberTwo = _wrap_->residueNumberTwo; ////map<const ChainResidueIndex, int,twoIndexCmp>
//	_struct_.temperatureArray = _wrap_->temperatureArray; //////vector<double>
//	_struct_.temperaturePriority = _wrap_->temperaturePriority; //////vector<int>
//	_struct_.dutyCycleArray = _wrap_->dutyCycleArray; //////vector<double>
//	_struct_.dutyCyclePriority = _wrap_->dutyCyclePriority; //////vector<int>
//	_struct_.globalCoulombScaleFactorArray = _wrap_->globalCoulombScaleFactorArray; /////*vector<double>
//	_struct_.globalCoulombScaleFactorPriority = _wrap_->globalCoulombScaleFactorPriority; ////vector<int>
//	_struct_.globalVdwScaleFactorArray = _wrap_->globalVdwScaleFactorArray; ////vector<double>
//	_struct_.globalVdwScaleFactorPriority = _wrap_->globalVdwScaleFactorPriority; ////vector<int>
//	_struct_._leontisWesthofClass = _wrap_->_leontisWesthofClass; ////LeontisWesthofClass
//	_struct_.userVariables = _wrap_->userVariables; ////mutable map<const String,double>
//	_struct_.myDensityMap = _wrap_->myDensityMap; ////DensityMap
//	_struct_.myElectroDensityMap = _wrap_->myElectroDensityMap; ////DensityMap
//	_struct_.mobilizerContainer = _wrap_->mobilizerContainer; ////MobilizerContainer
//	_struct_.physicsContainer = _wrap_->physicsContainer; ////PhysicsContainer
//	_struct_.constraintToGroundContainer = _wrap_->constraintToGroundContainer; ////ConstraintToGroundContainer
//	_struct_.displacementContainer = _wrap_->displacementContainer; ////DisplacementContainer
//	_struct_.atomSpringContainer = _wrap_->atomSpringContainer; ////AtomSpringContainer
//	_struct_.myBiopolymerClassContainer = _wrap_->myBiopolymerClassContainer; ////BiopolymerClassContainer
//	_struct_.moleculeClassContainer = _wrap_->moleculeClassContainer; ////MoleculeClassContainer
//	_struct_.waterDropletContainer = _wrap_->waterDropletContainer; ////WaterDropletContainer
//	_struct_.proteinSequences = _wrap_->proteinSequences; ////map<const String,String>
//	_struct_.coarseNucleicAcidSequences = _wrap_->coarseNucleicAcidSequences; ////map<const String,String>
//	_struct_.numRigidSegments = _wrap_->numRigidSegments; ////map<const String, int>
//	_struct_.firstResidueNumbersIterator = _wrap_->firstResidueNumbersIterator; ////map<const String,int>::iterator
	// //void addRingClosingBond(const String chainID, ResidueID residueID1, String atomName1,String bondCenterName1,  ResidueID residueID2, String atomName2,String bondCenterName2);
//	_struct_.bondCenterName2) = _wrap_->bondCenterName2); ////// //void addRingClosingBond(const String chainID, ResidueID residueID1, String atomName1,String bondCenterName1, ResidueID residueID2, String atomName2,String
	// void addC1pSprings (LeontisWesthofClass myLeontisWesthofClass);
//	_struct_.myLeontisWesthofClass) = _wrap_->myLeontisWesthofClass); ////// void addC1pSprings (LeontisWesthofClass
	// void applyAtomSprings (SimbodyMatterSubsystem & matter, GeneralForceSubsystem & forces);
//	_struct_.forces) = _wrap_->forces); ////// void applyAtomSprings (SimbodyMatterSubsystem & matter, GeneralForceSubsystem &
	// void configureDumm( DuMMForceFieldSubsystem & dumm);
//	_struct_.dumm) = _wrap_->dumm); ////// void configureDumm( DuMMForceFieldSubsystem &
	// static double myAtoF(map<const String,double> myUserVariables,const char* value );
//	_struct_.) = _wrap_->); ////// //bool chainIsBiopolymer(String myChainId
	// static bool aToBool( const String& name, const char* value );
	// static bool compareUpper( const String& param, const char* symbol );
//	_struct_.baseOperationVector = _wrap_->baseOperationVector; ////vector<BasePair>
//	_struct_.contactContainer = _wrap_->contactContainer; ////ContactContainer
//	_struct_.densityContainer = _wrap_->densityContainer; ////DensityContainer
//	_struct_.electroDensityContainer = _wrap_->electroDensityContainer; ////DensityContainer
//	_struct_.singleBondMobilityVector = _wrap_->singleBondMobilityVector; ////vector<SingleBondMobility>
//	_struct_.basePairPartnerVector = _wrap_->basePairPartnerVector; ////vector<BasePairPartner>
//	_struct_.includeAllNonBondAtomsInResidueVector = _wrap_->includeAllNonBondAtomsInResidueVector; //////vector<IncludeAllNonBondAtomsInResidue>
//	_struct_.includeAllResiduesWithinVector = _wrap_->includeAllResiduesWithinVector; ////vector<AllResiduesWithin>
//	_struct_.includeNonBondAtomInBiopolymerVector = _wrap_->includeNonBondAtomInBiopolymerVector; ////vector<IncludeNonBondAtomInBiopolymerStruct>
//	_struct_.waterDropletAboutResidueVector = _wrap_->waterDropletAboutResidueVector; ////vector <WaterDropletAboutResidueStruct>
//	_struct_.mobilizerDomainsInterfaceVector = _wrap_->mobilizerDomainsInterfaceVector; ////vector<MobilizerDomainsInterface>
	// void removeBasePairsInRigidStretch ();
//	_struct_.() = _wrap_->(); ////// // void initializeDefaults
	// void printAllSettings (   ostream  & myOstream = std::cout, String remarkString = "") ;
//	_struct_."") = _wrap_->""); ////// void printAllSettings ( ostream & myOstream = std::cout, String remarkString =
	// void removeNonPriorityBasePairs (int priorityLevel);
//	_struct_.priorityLevel) = _wrap_->priorityLevel); ////// void removeNonPriorityBasePairs (int
	// //int getFirstResidueNumbers(const String myChainId) const ;
//	_struct_.const = _wrap_->const; ////// // int getNumBasePairs()
	// // int getProteinFirstResidueNumbers(const String myProteinChainId) const ;
	// //int getBasePriority(int baseResidueNumber,String baseChain, String basePairingEdge) const ;
	// // int getNumBasePairs() const;
	// void updateBasePair(int index,
	// String ch1, int res1, String edge1,
	// String ch2, int res2, String edge2,
	// String orient);
	// void updateMobilizerStretch(int index,
	// String chainId,
	// int startRes,
	// int endRes,
	// String bondMobility);
	// void addAllResiduesWithin(String chainID, int resID, double radius);
//	_struct_.radius) = _wrap_->radius); ////// void updateAllResiduesWithin(int index, String chainID, int resID, double
	// void updateAllResiduesWithin(int index, String chainID, int resID, double radius);
	// void deleteAllResiduesWithin(int index);
//	_struct_.index) = _wrap_->index); ////// void deleteIncludeAllNonBondAtomsInResidue(int
	// void updateIncludeAllNonBondAtomsInResidue(int index, String chainID, int resID);
//	_struct_.resID) = _wrap_->resID); ////// void updateIncludeAllNonBondAtomsInResidue(int index, String chainID, int
	// void deleteIncludeAllNonBondAtomsInResidue(int index);
	// //int calcHighestPriority();
//	_struct_.calcHighestPriority() = _wrap_->calcHighestPriority(); ////// //int
	// //int calcLowestBondingResidue(const String myChainId) ;
//	_struct_.myChainId) = _wrap_->myChainId); ////// //bool chainIsMonoAtoms(String
	// //int calcHighestBondingResidue(const String myChainId);
	// void setLeontisWesthofBondRowIndex();
//	_struct_.setLeontisWesthofBondRowIndex() = _wrap_->setLeontisWesthofBondRowIndex(); ////// void
	// void parameterStringInterpreter(const String & paramstr);
//	_struct_.paramstr) = _wrap_->paramstr); ////// void parameterStringInterpreter(const String &
	// void parameterStringInterpreter(const ParameterStringClass & parameterStringClass,
	// const int readStage = 0,
	// const bool readAtOneStageOnly = false,
	// const bool readOnlyUntilStage = false,
	// const bool readExcept = false);
	// void initializeFromFileOnly(const char * parameterFileName = "./commands.dat" ) ;
	// void setFirstAndLastStage(const char * parameterFileName = "./commands.dat" ) ;
	// void loadSequencesFromPdb(const char * pdbFileName);
//	_struct_.pdbFileName) = _wrap_->pdbFileName); ////// void loadSequencesFromPdb(const char *
	// //void printRigidSegments();
//	_struct_.printRigidSegments() = _wrap_->printRigidSegments(); ////// //void
	// // void printBasePairs();
//	_struct_.printBasePairs() = _wrap_->printBasePairs(); ////// // void
	// // void printBaseAssignments();
//	_struct_.printBaseAssignments() = _wrap_->printBaseAssignments(); ////// // void
	// void postInitialize();
//	_struct_.postInitialize() = _wrap_->postInitialize(); ////// void
	// void clearContainers();
//	_struct_.clearContainers() = _wrap_->clearContainers(); ////// void
	// void clearBiopolymers();
//	_struct_.clearBiopolymers() = _wrap_->clearBiopolymers(); ////// void
	// void clearForces();
//	_struct_.clearForces() = _wrap_->clearForces(); ////// void
	// void clearConstraints();
//	_struct_.clearConstraints() = _wrap_->clearConstraints(); ////// void
	// // void initializeDefaults ();
	// void initializeDefaults(const char * leontisWesthofInFileName = "./parameters.csv");
//	_struct_."./parameters.csv") = _wrap_->"./parameters.csv"); ////// void initializeDefaults(const char * leontisWesthofInFileName =
	// void initialize(const char * parameterFileName = "./commands.dat" );
	// //bool chainIsBiopolymer(String myChainId );
	// //bool chainIsMonoAtoms(String myChainId);
	// //int getChainIndex(String myChainId , vector<Biopolymer> & tempChain);
//	_struct_.tempChain) = _wrap_->tempChain); ////// //int getChainIndex(String myChainId , vector<Biopolymer> &
//	_struct_.myMonoAtomsContainer = _wrap_->myMonoAtomsContainer; ////MonoAtomsContainer
	// //variables for internal use only:
	// int r;
	// int ti;
	// int gcsfi;
	// int gvsfi;
	// int d;
	// char * s;
	// //int numChains ;
	// //int numProteinChains ;
	// //int prioritize ;
	// //temperature = 300;
	// //outQVectorFileName;
	// //firstStage = 1;
	// //lastStage = 0;// calcHighestPriority();
//	_struct_.1 = _wrap_->1; //////dutyCycle =
//	_struct_.0 = _wrap_->0; //////priority =

}
