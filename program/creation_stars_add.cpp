 #include <iostream>
#include <vector>
#include <bitset>


#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/feature/FeatureAttributes.h>
#include <LeMonADE/feature/FeatureBondsetUnsaveCheck.h>
#include <LeMonADE/utility/RandomNumberGenerators.h>
#include <LeMonADE/analyzer/AnalyzerWriteBfmFile.h>
#include <LeMonADE/updater/UpdaterSimpleSimulator.h>
#include <LeMonADE/updater/UpdaterReadBfmFile.h>
#include <LeMonADE/utility/TaskManager.h>
#include <LeMonADE/feature/FeatureMoleculesIOUnsaveCheck.h>
#include <LeMonADE/feature/FeatureExcludedVolumeScIdOnLattice.h>
#include <LeMonADE/utility/LatticePredicates.h>

#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/feature/FeatureAttributes.h>
#include <LeMonADE/feature/FeatureBondset.h>
#include <LeMonADE/feature/FeatureExcludedVolumeSc.h>
#include <LeMonADE/feature/FeatureMoleculesIO.h>
#include <LeMonADE/utility/RandomNumberGenerators.h>
#include <LeMonADE/analyzer/AnalyzerWriteBfmFile.h>
#include <LeMonADE/utility/TaskManager.h>
#include <LeMonADE/updater/UpdaterReadBfmFile.h>

#include "../../updater/UpdaterAddstars.h"
#include "../../feature/FeatureNetworkInformation.h"



/******************************************************************************
 *
 ******************************************************************************/
int main(int argc, char* argv[])
{
  //default arguments and simulation parameters	
  
 
  std::string outfile;
  
  std::string infile; 
  int N_arms;
  int length_arm;
  uint32_t N_stars;
  uint32_t box_length;

  if(!((argc=6)))
{
  std::string errormessage;
		errormessage="usage: ./creation_stars N_stars N_arms box_length length_arm outfile\n";
		throw std::runtime_error(errormessage);
}
else
{
N_stars=atoi(argv[1]);
N_arms=atoi(argv[2]); 
box_length=atoi(argv[3]);
length_arm=atoi(argv[4]); // in our case length of the arm or the degree of polymersiation of the arm is counted by including the crosslink
outfile=argv[5];
}
    RandomNumberGenerators rng;
    rng.seedAll();

        uint32_t xBox(box_length);
	uint32_t yBox(box_length);
	uint32_t zBox(box_length);
	
	bool periodic(true);
	bool readIn(false);
	bool BoxSetX(false);
	bool BoxSetY(false);
	bool BoxSetZ(false);
	
	//define the fetures to be used
	
	typedef FeatureExcludedVolumeScIdOnLattice<FeatureLatticePowerOfTwo<uint32_t>, AttributeANDMonomerID > MyExcludedVolumeFeature;
    typedef LOKI_TYPELIST_4(FeatureMoleculesIOUnsaveCheck, 
			    FeatureMoleculesInformation, 
			    FeatureAttributes, 
			    MyExcludedVolumeFeature
			   ) Features;
    typedef ConfigureSystem<VectorInt3,Features, 7> Config;
    typedef Ingredients<Config> Ing;
	
	Ing myIngredients;

	 myIngredients.setBoxX(xBox);
		    myIngredients.setBoxY(yBox);
		    myIngredients.setBoxZ(zBox);
		    myIngredients.setPeriodicX(periodic);
		    myIngredients.setPeriodicY(periodic);
		    myIngredients.setPeriodicZ(periodic);
		    myIngredients.modifyBondset().addBFMclassicBondset();
		    myIngredients.synchronize(myIngredients);
		    


    //define the fetures to be used
        //set up the tasks
    TaskManager taskmanager;
    //Read in config 
    
    taskmanager.addUpdater(new UpdaterAddstars<Ing>(myIngredients,N_stars,N_arms,length_arm),1);
    taskmanager.addAnalyzer(new AnalyzerWriteBfmFile<Ing>(outfile.c_str(), myIngredients
			, AnalyzerWriteBfmFile<Ing>::APPEND),1);
    taskmanager.initialize();
    taskmanager.run(1);
    taskmanager.cleanup();
	   
	  
	  
  
  
  
  return 0;
	
}

