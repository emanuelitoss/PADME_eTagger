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
#include <vector>

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4UserLimits;

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:

    DetectorConstruction();
    virtual ~DetectorConstruction();

    virtual G4VPhysicalVolume* Construct();

    // getters for physical volumes
    const G4VPhysicalVolume* GetScintillator() const;
    const std::vector <G4VPhysicalVolume*> GetSiPMs() const;

    // scoring volume
    G4LogicalVolume* GetScoringVolume() const { return fScoringVolume; }

  private:

    // physical volumes
    G4VPhysicalVolume* fTagger;
    std::vector <G4VPhysicalVolume*> fSiPMs
      = {nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr};

    G4UserLimits* fStepLimit;

  protected:

    G4LogicalVolume* fScoringVolume;

    G4Material* CreatePlasticMaterial() const;
    G4Material* CreatePyrex() const;
    G4Material* CreateLowVacuumAir() const;

    void OpticalSurfaceTagger_SiPM(G4VPhysicalVolume*, G4VPhysicalVolume*) const;
    void OpticalSurfaceTagger_Vacuum(G4VPhysicalVolume*, G4VPhysicalVolume*) const;

};

// inline functions

inline const G4VPhysicalVolume* DetectorConstruction::GetScintillator() const { 
  return fTagger;
}

inline const std::vector <G4VPhysicalVolume*> DetectorConstruction::GetSiPMs() const {
  return fSiPMs;
}

#endif