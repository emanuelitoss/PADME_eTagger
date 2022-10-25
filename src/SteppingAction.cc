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
  fScoringVolume(nullptr),
  fDetConstruction(detConstruction)
{}

SteppingAction::~SteppingAction(){}

void SteppingAction::UserSteppingAction(const G4Step* step){

  if (!fScoringVolume) {
    const DetectorConstruction* detectorConstruction
      = static_cast<const DetectorConstruction*> (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    fScoringVolume = detectorConstruction->GetScoringVolume();   
  }

  auto runData = static_cast<RunData*> (G4RunManager::GetRunManager()->GetNonConstCurrentRun());

  G4VPhysicalVolume* PostStepPV = step->GetPostStepPoint()->GetTouchableHandle()->GetVolume();
  G4VPhysicalVolume* PreStepPV = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
  G4LogicalVolume* volume = PreStepPV->GetLogicalVolume();

  // Store of arrival time for optical photons
  G4bool IsOpticalPhoton = ( step->GetTrack()->GetParticleDefinition() == G4OpticalPhoton::OpticalPhotonDefinition() );
  G4bool IsPhotDetectedInSiPM;

  // I add a random noise extrapolated from a Gaussian with 0.6 ns std deviation
  CLHEP::RandGauss noise_generator(new CLHEP::MTwistEngine(), 0*ns, NOISE_STD_DEV*ns);

  if(IsOpticalPhoton){
    for(int id = 0; id < 8; ++id){

      IsPhotDetectedInSiPM = (PreStepPV == fDetConstruction->GetScintillator()) && (PostStepPV == fDetConstruction->GetSiPMs()[id]);
      if (IsPhotDetectedInSiPM)
      {

        // get global time
        G4double G_arrival_time = step->GetTrack()->GetGlobalTime();

        // adding noise
        G_arrival_time += noise_generator.fire();

        runData->FillTimePerPhoton(id, G_arrival_time);
        step->GetTrack()->SetTrackStatus(fStopAndKill);
        fEventAction->SetMinTimeIfLess(id, G_arrival_time);

      }
    }
  }

  // collect energy deposited in this step
  G4double edepStep = step->GetTotalEnergyDeposit();

  if ( PreStepPV == fDetConstruction->GetScintillator() ) {
    runData->AddEnergy(kBGO, edepStep);
  }

  // check if we are not in scoring volume
  if (volume != fScoringVolume) return;
  // else score
  fEventAction->AddEdep(edepStep);

}
