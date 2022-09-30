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
/// \file ActionInitialization.cc
/// \brief Implementation of the ActionInitialization class

#include "../include/ActionInitialization.hh"
#include "../include/PrimaryGeneratorAction.hh"
#include "../include/RunAction.hh"
#include "../include/EventAction.hh"
#include "../include/SteppingAction.hh"

#include "DetectorConstruction.hh"
#include "G4MTRunManager.hh"

#include <fstream>
#include <iostream>
using namespace std;

ActionInitialization::ActionInitialization(DetectorConstruction* detConstruction)
 : G4VUserActionInitialization(),
 fDetConstruction(detConstruction)
{}

ActionInitialization::~ActionInitialization(){}

void ActionInitialization::BuildForMaster() const {
  RunAction* runAction = new RunAction;
  SetUserAction(runAction);
}

void ActionInitialization::Build() const {

  G4String outputFile = "../eTag_data/positions.txt";

  SetUserAction(new PrimaryGeneratorAction(outputFile));

  RunAction* runAction = new RunAction;
  SetUserAction(runAction);

  EventAction* eventAction = new EventAction(runAction);
  SetUserAction(eventAction);

  SetUserAction(new SteppingAction(eventAction, fDetConstruction));

}