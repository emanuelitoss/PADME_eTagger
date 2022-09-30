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
/// \file EventAction.hh
/// \brief Definition of the EventAction class

#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class RunAction;

/// Event action class
class EventAction : public G4UserEventAction
{

  public:

  EventAction(RunAction* runAction);
  virtual ~EventAction();
  
  virtual void BeginOfEventAction(const G4Event* event);
  virtual void EndOfEventAction(const G4Event* event);

  // to check the passage in the detectors
  void PassedThroughBGO();
  void PassedThroughScint1();
  void PassedThroughScint2();
  void DetectionInPMT1();
  void DetectionInPMT2();
  G4bool BoolTrigger1() const { return IsInTrg1; }
  G4bool BoolTrigger2() const { return IsInTrg2; }
  G4bool BoolPMT1() const { return PMT1detection; }
  G4bool BoolPMT2() const { return PMT2detection; }
  G4bool Bool() const { return IsInBGO; }

  // add energy losses in materials 
  void AddEdep(G4double edep);
  void AddEdepBGO(G4double edep);
  void AddEdepScint1(G4double edep);
  void AddEdepScint2(G4double edep);
  void AddEdepBGOCerenkov(G4double edep);
  void AddEdepBGOScint(G4double edep);

  void AddProducedCerenkovPhoton() { ++Nproduced_Cerenkov; }
  void AddProducedScintillationPhoton() { ++Nproduced_Scintillation; }

  void PrintStatus();

  private:

  RunAction* fRunAction;
  
  // deposited energies in: ScoringVolume, BGO crystal and two plastic Scintillators
  G4double fEdep;
  G4double fEdep_BGO;
  G4double fEdep_Scint1;
  G4double fEdep_Scint2;
  G4double fEdep_BGO_Cherenkov;
  G4double fEdep_BGO_Scintillation;
  G4int Nphotons_Cerenkov;
  G4int Nphotons_Scint;
  G4int Nproduced_Cerenkov;
  G4int Nproduced_Scintillation;
  
  // boolean variables to check if the particle passes through physical volumes
  G4bool IsInBGO = false;
  G4bool IsInTrg1 = false;
  G4bool IsInTrg2 = false;
  G4bool PMT1detection = false;
  G4bool PMT2detection = false;
  
};

// inline functions

inline void EventAction::PassedThroughBGO(){
  // old: if(PMT1detection && PMT2detection) IsInBGO = true;
  if(PMT1detection) IsInBGO = true;
}

inline void EventAction::PassedThroughScint1(){
  IsInTrg1 = true;
}

inline void EventAction::PassedThroughScint2(){
  IsInTrg2 = true;
}

inline void EventAction::DetectionInPMT1(){
  PMT1detection = true;
}

inline void EventAction::DetectionInPMT2(){
  PMT2detection = true;
}

#endif
