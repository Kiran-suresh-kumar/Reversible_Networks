  /*--------------------------------------------------------------------------------
    ooo      L   attice-based  |
  o\.|./o    e   xtensible     | LeMonADE: An Open Source Implementation of the
 o\.\|/./o   Mon te-Carlo      |           Bond-Fluctuation-Model for Polymers
oo---0---oo  A   lgorithm and  |
 o/./|\.\o   D   evelopment    | Copyright (C) 2013-2015 by
  o/.|.\o    E   nvironment    | LeMonADE Principal Developers
    ooo                        |
----------------------------------------------------------------------------------

This file is part of LeMonADE.

LeMonADE is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

LeMonADE is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with LeMonADE.  If not, see <http://www.gnu.org/licenses/>.

--------------------------------------------------------------------------------*/

#ifndef LEMONADE_UPDATER_SETUP_stars_
#define LEMONADE_UPDATER_SETUP_stars
/**
 * @file
 *
 * @class UpdaterAddstars
 *
 * @brief Updater to create a solution of four-armed stars.
 *
 * @details This is a simple implementation of a system setup starting from an empty ingredients
 * or a system with some monomers inside. This updater requires FeatureAttributes.
 * Three tags are added to the monomers, one for the central monomer(attribute tag- 1), one for the end monomer /reactive monomer (attribute tag- 3) and another for the middle monomers (attribute tag- 2) .
 *
 * @tparam IngredientsType
 *
 * @param ingredients_ The system, holding eigther an empty simulation box for system setup
 * @param N_stars_ number of stars that should be added to the box
 * @param N_arms_ number of arms for each star
 * @param length_arm_ length of each star arm
 * @param type1_ attribute tag of "middle" monomers
 * @param type2_ attribute tag of "central" monomer
 * @param type_3_  attribute tag of "reactive(end)" monomers
 **/

#include <LeMonADE/updater/UpdaterAbstractCreate.h>
#include <LeMonADE/utility/Vector3D.h>
#include <cmath>

template<class IngredientsType>
class UpdaterAddstars: public UpdaterAbstractCreate<IngredientsType>
{
  typedef UpdaterAbstractCreate<IngredientsType> BaseClass;

public:
  UpdaterAddstars(IngredientsType& ingredients_, uint32_t N_stars_, uint32_t N_arms_, uint32_t length_arm_, int32_t type1_=2, int32_t type2_=3,int32_t type3_=1, bool IsSolvent=false);

  virtual void initialize();
  virtual bool execute();
  virtual void cleanup();

  //! getter function for write compressed solvent bool
  const bool getIsSolvent() const {return IsSolvent;}
  
  //! getter function for calculated density
  const double getDensity() const {return density;}

private:
  // provide access to functions of UpdaterAbstractCreate used in this updater
  using BaseClass::ingredients;
  using BaseClass::addMonomerToParent;
  using BaseClass::addSingleMonomer;
  using BaseClass::linearizeSystem;

  //! number of stars in the box
  uint32_t N_stars;

  //! number of arms of the stars in the box
  uint32_t N_arms;

  //! length of the arm of the stars
  uint32_t length_arm;
  
  //! lattice occupation density
  double density;

  //! bool for execution
  bool wasExecuted;

  //! attribute tag of inner monomers
  int32_t type1;

  //! getAttributeTag of central monomers
  int32_t type2;

  //! getAttributeTag of reactive Monomers, the reactive one
  int32_t type3;
  
  //! bool to check if chains of size 1 should be compressed to solvent
  bool IsSolvent;
};

/**
* @brief Constructor handling the new systems paramters
*
* @param ingredients_ a reference to the IngredientsType - mainly the system
* @param N_arms number of arms of star
* @param N_stars number of stars in the box
* @param length_arm degree of polymerisation of the star arm
*/
template < class IngredientsType >
UpdaterAddstars<IngredientsType>::UpdaterAddstars(IngredientsType& ingredients_, uint32_t N_stars_, uint32_t N_arms_, uint32_t length_arm_, int32_t type1_, int32_t type2_,int32_t type3_, bool IsSolvent_):
BaseClass(ingredients_), N_stars(N_stars_), N_arms(N_arms_),length_arm(length_arm_), type3(type3_), density(0.0), wasExecuted(false),
type1(type1_), type2(type2_), IsSolvent(IsSolvent_)
{}

/**
* @brief initialise function, calculate the target density to compare with at the end.
*
* @tparam IngredientsType Features used in the system. See Ingredients.
*/
template < class IngredientsType >
void UpdaterAddstars<IngredientsType>::initialize()
{
  std::cout<< "initialize UpdaterAddStars" << std::endl;

  // get the target density from the sum of existing monomers and the newly added stars
  density=(double)( ingredients.getMolecules().size() + (((N_stars)*(N_arms)*(length_arm))+N_stars) ) * 8  /(double)( ingredients.getBoxX()*ingredients.getBoxY()*ingredients.getBoxZ() );

  std::cout << "add "<<(((N_stars)*(N_arms)*(length_arm))+N_stars)<<" monomers to the box"<<std::endl;

  execute();
}

/**
* @brief Execution of the system creation
*
* @tparam IngredientsType Features used in the system. See Ingredients.
*/
template < class IngredientsType >
bool UpdaterAddstars<IngredientsType>::execute()
{
  if(wasExecuted)
    return true;

  std::cout << "execute UpdaterAddstars" << std::endl;

  //loop over No. of stars and arms of stars and length of the arm and build it up


  for(uint32_t i_star=0; i_star<N_stars; i_star++)  // N_stars is the no. of stars
  {
    addSingleMonomer(type2);// type 2 is for the center monomer
  for(uint32_t j_=0; j_<N_arms;j_++) // N_arms is the number of arms
  {
    addMonomerToParent((((length_arm)*N_arms)+1)*i_star,type1);   // for adding the new arm to the center monomer each time
  for(uint32_t j_1=2; j_1<=(length_arm);j_1++) 
  {
    if(j_1==(length_arm))
    {
      addMonomerToParent(ingredients.getMolecules().size()-1, type3); // reactive end monomers are labelled with type 3
    }
    else
    { 
     addMonomerToParent(ingredients.getMolecules().size()-1, type1);  // type 1 for the inner monomer of the arm
    }
  }
  
    
  }
    
  }
  
  
  
  
  
  
  ingredients.synchronize();
  double lattice_volume(ingredients.getBoxX()*ingredients.getBoxY()*ingredients.getBoxZ());
  
  if(std::abs(density - ( (double)(ingredients.getMolecules().size()*8) / lattice_volume )) > 0.0000000001 )
  {
    std::cout << density << " " <<( (ingredients.getMolecules().size()*8) / lattice_volume)<<std::endl;
    throw std::runtime_error("some issue with the density, plese check that!");
  }
  else
  {
    std::cout << "real lattice occupation density =" << (8*ingredients.getMolecules().size()) / lattice_volume<<std::endl;
    wasExecuted=true;
    return true;
  }
}

/**
* @brief Standard clean up.
*
* @tparam IngredientsType Features used in the system. See Ingredients.
*/
template < class IngredientsType >
void UpdaterAddstars<IngredientsType>::cleanup()
{

}


#endif /* LEMONADE_UPDATER_SETUP_ADD_STARS */