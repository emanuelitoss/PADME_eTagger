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

#include "../include/SteppingAction.hh"
#include "../include/EventAction.hh"
#include "../include/DetectorConstruction.hh"
#include "../include/RunData.hh"
#include "../include/OutputColors.hh"

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

SteppingAction::SteppingAction(EventAction* eventAction, const DetectorConstruction* detConstruction)
: G4UserSteppingAction(),
  fEventAction(eventAction),
  fScoringVolume(0),
  fDetConstruction(detConstruction)
{}

SteppingAction::~SteppingAction(){}

void SteppingAction::UserSteppingAction(const G4Step* step){

  if (!fScoringVolume) {
    const DetectorConstruction* detectorConstruction
      = static_cast<const DetectorConstruction*>
        (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    fScoringVolume = detectorConstruction->GetScoringVolume();   
  }

  // get volume of the current step
  G4LogicalVolume* volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
  G4VPhysicalVolume* physicalVolume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();

  /////////////////////////////////////////////// TRIGGER //////////////////////////////////////////////
  
  // check that it is a e-
  G4Track* primary = step->GetTrack();
  G4bool check_electron = primary->GetParticleDefinition()->GetParticleName() == "e-";

  //////////////////////////////////////////// SCORE ENERGY ////////////////////////////////////////////

  // collect energy deposited in this step
  G4double edepStep = step->GetTotalEnergyDeposit();

  auto runData = static_cast<RunData*> (G4RunManager::GetRunManager()->GetNonConstCurrentRun());

  // for each PV:
  // modify a boolean value to check the passage of particle through PV and add energy
  if ( physicalVolume == fDetConstruction->GetScintillator() ) {
    fEventAction->PassedThroughBGO();
    runData->Add(kBGO, edepStep);
    fEventAction->AddEdepBGO(edepStep);
  }

  // Cherenkov and scintillation photons in the two PMTs
  
  if(!check_electron) // if it is not an e-
  {
    G4String creator_process_thisparticle = primary->GetCreatorProcess()->GetProcessName();

    // boolean variables
    // a photon is red if passes through the surface between BGO and PMT
    // notice: G4string::compare returns 0 if it is true and 1 if it is false.
    //          >> additive NOT ! to get coherent definitions
    G4bool check_cerenkov = !creator_process_thisparticle.compare("Cerenkov");
    G4bool check_scintillation = !creator_process_thisparticle.compare("Scintillation");

    G4bool prestep = step->GetPreStepPoint()->GetPhysicalVolume() == fDetConstruction->GetScintillator();
    G4bool poststep_cerenkovPMT = step->GetPostStepPoint()->GetPhysicalVolume() == fDetConstruction->GetCerenkovVolume();
    G4bool poststep_scintillationPMT = step->GetPostStepPoint()->GetPhysicalVolume() == fDetConstruction->GetScintillatorVolume();

    if(check_cerenkov && prestep && poststep_cerenkovPMT){
      // add a detection in PMT1 - Cherenkov light - boolean variable
      fEventAction->DetectionInPMT1();

      // add a photon and its energy
      G4double cher_photon_energy = primary->GetKineticEnergy();
      runData->Add(kBGO_Cherenkov, cher_photon_energy);
      fEventAction->AddEdepBGOCerenkov(cher_photon_energy);

      // delete the photon in order to avoid double counting
      primary->SetKineticEnergy(0.);
      primary->SetTrackStatus(fStopAndKill);
    }

    if(check_scintillation && prestep && poststep_scintillationPMT){
      // add a detection in PMT2 - scintillation - boolean variable
      fEventAction->DetectionInPMT2();

      // add a photon and its energy
      G4double scint_photon_energy = primary->GetKineticEnergy();
      runData->Add(kBGO_Scintillation, scint_photon_energy);
      fEventAction->AddEdepBGOScint(scint_photon_energy);

      // delete the photon in order to avoid double counting
      primary->SetKineticEnergy(0.);
      primary->SetTrackStatus(fStopAndKill);
    }
  } else {
    // if we are in the BGO
    if(step->GetPreStepPoint()->GetPhysicalVolume() == fDetConstruction->GetScintillator()){
    const std::vector <const G4Track*> * secondaries = step->GetSecondaryInCurrentStep();
      // loop over secondaries
      for( auto sec : *secondaries ){

        G4String creator_process = sec->GetCreatorProcess()->GetProcessName();

        if(!creator_process.compare("Cerenkov")) fEventAction->AddProducedCerenkovPhoton();
        else if(!creator_process.compare("Scintillation")) fEventAction->AddProducedScintillationPhoton();
      }
    }
  }

  // check if we are not in scoring volume
  if (volume != fScoringVolume) return;
  // else score
  fEventAction->AddEdep(edepStep);

}
