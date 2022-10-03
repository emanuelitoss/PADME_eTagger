//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file DetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4LogicalVolume.hh"
#include "globals.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4UserLimits;

/// Detector construction class to define materials and geometry.

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:

    DetectorConstruction();
    virtual ~DetectorConstruction();

    virtual G4VPhysicalVolume* Construct();

    // getters for physical volumes of the apparatus
    const G4VPhysicalVolume* GetScintillator() const;
    const G4VPhysicalVolume* GetPlasticScintillator_1() const;
    const G4VPhysicalVolume* GetPlasticScintillator_2() const;
    const G4VPhysicalVolume* GetCerenkovVolume() const;
    const G4VPhysicalVolume* GetScintillatorVolume() const;

    // return the scoring volume
    G4LogicalVolume* GetScoringVolume() const { return fScoringVolume; }

  private:

    // physical volumes of detector
    G4VPhysicalVolume* fScintillator;
    G4VPhysicalVolume* fPlasticScintillator_1;
    G4VPhysicalVolume* fPlasticScintillator_2;
    G4VPhysicalVolume* fCerenkovPMT;
    G4VPhysicalVolume* fScintillatorPMT;
    
    G4double minimal_radius;
    
    G4UserLimits* fStepLimit;            // pointer to user step limits

  protected:

    G4LogicalVolume* fScoringVolume;

    // materials of the experiment
    G4Material* CreatePlasticMaterial() const;
    G4Material* CreatePyrex() const;

    // set the optical surface between the scintillator and SiPMs
    void OpticalSurfacePlastic_SiPM(G4VPhysicalVolume*, G4VPhysicalVolume*) const;

};

// inline functions

inline const G4VPhysicalVolume* DetectorConstruction::GetScintillator() const { 
  return fScintillator; 
}

inline const G4VPhysicalVolume* DetectorConstruction::GetPlasticScintillator_1() const { 
  return fPlasticScintillator_1; 
}

inline const G4VPhysicalVolume* DetectorConstruction::GetPlasticScintillator_2() const { 
  return fPlasticScintillator_2; 
}

inline const G4VPhysicalVolume* DetectorConstruction::GetCerenkovVolume() const { 
  return fCerenkovPMT; 
}

inline const G4VPhysicalVolume* DetectorConstruction::GetScintillatorVolume() const { 
  return fScintillatorPMT; 
}

#endif