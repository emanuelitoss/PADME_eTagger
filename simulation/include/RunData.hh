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
/// \file RunData.hh
/// \brief Definition of the B4b::RunData class
  
#ifndef B4bRunData_h
#define B4bRunData_h 1

#include "G4Run.hh"
#include "globals.hh"
#include "g4root.hh"

#include "OutputColors.hh"
#include "RunAction.hh"

#include <array>

///  Run data class
///
/// It defines data members to hold the energy deposit and track lengths
/// of charged particles in Absober and Gap layers.
///
/// The data are collected step by step in SteppingAction, and
/// the accumulated values are filled in histograms and entuple
/// event by event in EventAction.

class RunData : public G4Run
{
    public:
        RunData();
        virtual ~RunData();

        void FillPerEvent(std::vector <G4double> charges, G4double x, G4double y);
        void FillTimePerPhoton(G4int id, G4double t);
        void FIllFirstTimes(std::vector <G4double> time);

    private:
};
  
#endif