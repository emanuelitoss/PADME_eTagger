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
#include "G4SystemOfUnits.hh"

#include "Randomize.hh"
#include <iomanip>

EventAction::EventAction(RunAction* runAction)
: G4UserEventAction(),
  fRunAction(runAction),
  fEdep(0.)
{
  min_times = {0.,0.,0.,0.,0.,0.,0.,0.};
  signals_charges = {0.,0.,0.,0.,0.,0.,0.,0.};
}

EventAction::~EventAction(){}

void EventAction::BeginOfEventAction(const G4Event*){

  fEdep = 0.;

  auto runData = static_cast<RunData*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
  runData->Reset();

  G4double max_t = 2000.*ns;
  min_times = {max_t,max_t,max_t,max_t,max_t,max_t,max_t,max_t};

}

void EventAction::EndOfEventAction(const G4Event* event){
  
  // accumulate statistics in run action
  fRunAction->AddEdep(fEdep);

  auto eventID = event->GetEventID();
  if (( eventID % 10000 == 0 )) {
    G4cout << "-------> End of event: " << eventID << G4endl;
  }

  auto runData = static_cast<RunData*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());

  for(int i = 0; i<8; ++i){
    runData->FIllFirstTimes(i, min_times[i]);
  }

  runData->FillPerEvent(signals_charges);

}

void EventAction::AddEdep(G4double edep){
  fEdep += edep;
}

void EventAction::SetMinTimeIfLess(G4int channel, G4double time){
  if(time != 0 && time < min_times[channel])  min_times[channel] = time;
}