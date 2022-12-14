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
/// \file SteppingAction.cc
/// \brief Implementation of the SteppingAction class

#include "SteppingAction.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"
#include "RunData.hh"
#include "OutputColors.hh"

#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"

#include "G4Cerenkov.hh"
#include "G4EventManager.hh"
#include "G4Scintillation.hh"
#include "G4OpticalPhoton.hh"
#include "G4ProcessManager.hh"
#include "G4SteppingManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4UnitsTable.hh"

#include <vector>

#include "Randomize.hh"
#define NOISE_STD_DEV 0.6 //[ns]

SteppingAction::SteppingAction(EventAction* eventAction, const DetectorConstruction* detConstruction)
: G4UserSteppingAction(),
  fEventAction(eventAction),
  fnoise_generator(new CLHEP::MTwistEngine(), 0, NOISE_STD_DEV),
  fDetConstruction(detConstruction)
{}

SteppingAction::~SteppingAction(){}

void SteppingAction::UserSteppingAction(const G4Step* step){

  auto runData = static_cast<RunData*> (G4RunManager::GetRunManager()->GetNonConstCurrentRun());

  G4VPhysicalVolume* PostStepPV = step->GetPostStepPoint()->GetTouchableHandle()->GetVolume();
  G4VPhysicalVolume* PreStepPV = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();

  // Store of arrival time for optical photons
  G4bool IsOpticalPhoton = ( step->GetTrack()->GetParticleDefinition() == G4OpticalPhoton::OpticalPhotonDefinition() );
  G4bool IsPhotDetectedInSiPM;

  G4double Global_arrival_time = 0., rndm_noise;
  G4double phot_energy = 0.*eV;

  if(IsOpticalPhoton){
  
    //G4StepPoint* endPoint   = step->GetPostStepPoint();
    //G4StepPoint* startPoint = step->GetPreStepPoint();
    //const G4VProcess* pds = endPoint->GetProcessDefinedStep();
    //if(pds->GetProcessName() && (pds->GetProcessName()!="Transportation") )
    //{
    //  G4cout << "Endpoint = " << endPoint->GetPosition() << "\t\t" << pds->GetProcessName() << endl;
    //}

    // get the photon energy
    phot_energy = step->GetTrack()->GetKineticEnergy();

    // loop over SiPMs
    for(int id = 0; id < 8; ++id){

      IsPhotDetectedInSiPM = (PreStepPV == fDetConstruction->GetScintillator()) && (PostStepPV == fDetConstruction->GetSiPMs()[id]);
      if(!IsPhotDetectedInSiPM) continue;

      IsPhotDetectedInSiPM = IsPhotDetectedInSiPM && ( ApplyDetectionEfficiency(phot_energy) );

      if(IsPhotDetectedInSiPM)
      {
        // get global time
        Global_arrival_time = step->GetTrack()->GetGlobalTime();

        // add noise
        rndm_noise =  (fnoise_generator.fire())*ns;

        Global_arrival_time += rndm_noise;

        runData->FillTimePerPhoton(id, Global_arrival_time);
        fEventAction->SetMinTimeIfLess(id, Global_arrival_time);
        fEventAction->SetDetectedPhoton(id);
      }

      // kill the photon in order to avoid double counting
      step->GetTrack()->SetTrackStatus(fStopAndKill);

    }
  }
}

G4bool SteppingAction::ApplyDetectionEfficiency(G4double photon_energy){

  G4double conversionFactor = 1239.8*nm*eV;
  G4double wavelength = conversionFactor/photon_energy;
  
  // Here photodetection efficiency of the SiPMs is defined.
  // See datasheet of Hamamatsu 13360-3050PE

  std::vector <G4double> bin_borders = {300.*nm, 325.*nm, 350.*nm, 375.*nm, 400.*nm, 425.*nm, 450.*nm, 475.*nm,
                                        500.*nm, 525.*nm, 550.*nm, 575.*nm, 600.*nm, 625.*nm, 650.*nm, 675.*nm,
                                        700.*nm, 725.*nm, 750.*nm, 775.*nm, 800.*nm, 825.*nm, 850.*nm, 875.*nm,
                                        900.*nm}; // 25 entries

  std::vector <G4double> efficiencies = {0.0125, 0.115, 0.250, 0.330, 0.367, 0.387, 0.405, 0.390, 0.375, 0.345,
                                         0.3125, 0.278, 0.250, 0.216, 0.185, 0.162, 0.140, 0.121, 0.105, 0.0875,
                                         0.0761, 0.060, 0.047, 0.0375}; // 25-1 = 24 entries.

  G4bool detection = 0;
  G4bool check_right_bin = 0;
  G4int idx;

  if(wavelength < 300.*nm || wavelength >= 900.*nm) detection = 0;
  else{

    idx = (int)((wavelength-300*nm)/(25*nm));

    if(wavelength >= bin_borders[idx] && wavelength < bin_borders[idx+1]) check_right_bin = 1;

    if(G4UniformRand() < efficiencies[idx]) detection = 1;
    else detection = 0;

  }

  return detection && check_right_bin;

}