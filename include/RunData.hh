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

#include "../include/OutputColors.hh"

#include <array>

const G4int kBGO = 0;
//const G4int kScint1 = 1;
//const G4int kScint2 = 2;
const G4int kBGO_Cherenkov = 1;
const G4int kBGO_Scintillation = 2;
//const G4int kNum_Cerenkov = 5;
//const G4int kNum_Scint = 6;
const G4int kDim = 3;//7;

const G4int kDimVolumes = 5;

///  Run data class
///
/// It defines data members to hold the energy deposit and track lengths
/// of charged particles in Absober and Gap layers.
///
/// In order to reduce the number of data members a 2-dimensions array
/// is introduced for each quantity:
/// - fEdep[], fTrackLength[].
///
/// The data are collected step by step in SteppingAction, and
/// the accumulated values are filled in histograms and entuple
/// event by event in EventAction.

class RunData : public G4Run
{
    public:
        RunData();
        virtual ~RunData();

        void Add(G4int id, G4double de);
        void FillPerEvent();

        void Reset();

        // Get methods
        G4String GetVolumeName(G4int id) const;
        G4double GetEdep(G4int id) const;

    private:
        std::array<G4String, kDimVolumes> fVolumeNames = { "BGO", "Scintillator_1", "Scintillator_2", "Cherenkov PMT", "Scintillation PMT"};
        std::array<G4double, kDim> fEdep = { 0., 0., 0.}; //, 0., 0., 0., 0. };
};

// inline functions
inline void RunData::Add(G4int id, G4double de) {
    fEdep[id] += de;
}

inline G4String RunData::GetVolumeName(G4int id) const {
    return fVolumeNames[id];
}

inline G4double RunData::GetEdep(G4int id) const {
    return fEdep[id];
}

  
#endif