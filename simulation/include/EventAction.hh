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
#include <vector>

class RunAction;

/// Event action class
class EventAction : public G4UserEventAction
{

  public:

  EventAction();
  virtual ~EventAction();
  
  virtual void BeginOfEventAction(const G4Event* event);
  virtual void EndOfEventAction(const G4Event* event);

  void SetMinTimeIfLess(G4int channel, G4double time);
  void SetDetectedPhoton(G4int id);
  G4double GetMinTime(G4int channel) { return min_times[channel]; }
  std::vector <G4double> GetSignalCharges() { return signals_charges; }
  
  private:

  std::vector <G4double> min_times;
  std::vector <G4double> signals_charges;
  G4double position_x, position_y;
  
};

inline void EventAction::SetDetectedPhoton(G4int id){
  signals_charges[id]++;
}

#endif
