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

  // get volume of the current step
  G4LogicalVolume* volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
  G4VPhysicalVolume* PreStepPV = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
  G4VPhysicalVolume* PostStepPV = step->GetPostStepPoint()->GetTouchableHandle()->GetVolume();

  // collect energy deposited in this step
  G4double edepStep = step->GetTotalEnergyDeposit();

  G4bool IsPhotDetected = (PreStepPV == fDetConstruction->GetScintillator()) && (PostStepPV == fDetConstruction->GetSiPMs()[0]);
  if (IsPhotDetected)
  {
    G4double L_arrival_time = step->GetTrack()->GetLocalTime();
    G4double G_arrival_time = step->GetTrack()->GetGlobalTime();
    std::cout << OBOLDCYAN
      << "Local arrival time:\t" << G4BestUnit(L_arrival_time,"Time")
      << "\tGlobal arrival time:\t" << G4BestUnit(G_arrival_time,"Time") << ORESET << std::endl;
  }
  else {
  }

  if ( PreStepPV == fDetConstruction->GetScintillator() ) {
    fEventAction->PassedThroughScintillator();
    runData->Add(kBGO, edepStep);
    fEventAction->AddEdepScintillator(edepStep);
  }

  // check if we are not in scoring volume
  if (volume != fScoringVolume) return;
  // else score
  fEventAction->AddEdep(edepStep);

}
