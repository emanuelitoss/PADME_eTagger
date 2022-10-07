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
/// \file EventAction.cc
/// \brief Implementation of the EventAction class

#include "EventAction.hh"
#include "RunAction.hh"
#include "RunData.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"

#include "Randomize.hh"
#include <iomanip>

EventAction::EventAction(RunAction* runAction)
: G4UserEventAction(),
  fRunAction(runAction),
  fEdep(0.)
{}

EventAction::~EventAction(){}

void EventAction::BeginOfEventAction(const G4Event*){

  fEdep = 0.;
  fEdep_BGO = 0.;
  fEdep_BGO_Cherenkov = 0.;
  fEdep_BGO_Scintillation = 0.;
  Nphotons_Cerenkov = 0;
  Nphotons_Scint = 0;

  IsInBGO = false;

  auto runData = static_cast<RunData*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
  runData->Reset();

}

void EventAction::EndOfEventAction(const G4Event* event){
  
  // accumulate statistics in run action
  fRunAction->AddEdep(fEdep);

  auto eventID = event->GetEventID();
  if (( eventID % 10000 == 0 )) {
    G4cout << "-------> End of event: " << eventID << G4endl;
  }

  if ( IsInBGO ){
    
    auto runData = static_cast<RunData*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
    runData->FillPerEvent();

    // adding a particle
    fRunAction->addDetectedParticle();
    
  }

}

void EventAction::AddEdep(G4double edep){
  fEdep += edep;
}

void EventAction::AddEdepScintillator(G4double edep){
  fEdep_BGO += edep;
}