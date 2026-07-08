/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2020-2023 OpenCFD Ltd.
    Copyright (c) 2026 Egor Pavlov
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::kodesChemistryModel

Description
    Custom chemistry model that extends StandardChemistryModel with
    additional logging capabilities.

SourceFiles
    kodesChemistryModel.C

\*---------------------------------------------------------------------------*/

#ifndef kodesChemistryModel_H
#define kodesChemistryModel_H

#include "StandardChemistryModel.H"

#include "basic_linalg.cuh"
#include "SeulexDeviceResources.cuh"
#include "pyJacSystem.cuh"
#include "HostResources.cuh"
#include "Operator.cuh"
#include "Seulex.cuh"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cuda/cmath>
#include <cuda_runtime.h>
#include "gpu_memory.cuh"
#include "mechanism.cuh"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class kodesChemistryModel Declaration
\*---------------------------------------------------------------------------*/

template<class ReactionThermo, class ThermoType>
class kodesChemistryModel
:
    public StandardChemistryModel<ReactionThermo, ThermoType>
{
    kodes::HostResources HostResource;
    kodes::HostResources HostResourceCopy;

    mechanism_memory *HostpyJacMechnismMemory;
    mechanism_memory *DevicepyJacMechnismMemory = nullptr;

    kodes::SeulexDeviceResources   DeviceResourceHostCopy;
    kodes::SeulexDeviceResources*  DeviceResourcePrt;

    kodes::pyJacSystem* pyJacSystemPrt;

    kodes::Operator<kodes::HostResources, kodes::SeulexDeviceResources> ResourceOperator;

    kodes::Seulex<kodes::pyJacSystem>* Integrator;

    label inertIndex;

public:

    //- Runtime type information
    TypeName("kodesChemistryModel");


    // Constructors

        //- Construct from thermo
        kodesChemistryModel(ReactionThermo& thermo);


    //- Destructor
    virtual ~kodesChemistryModel();

    // Member Functions

    virtual scalar solve(const scalar deltaT) override;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "kodesChemistryModel.C"
#endif

#endif

// ************************************************************************* //


