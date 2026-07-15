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

\*---------------------------------------------------------------------------*/

#include "kodesChemistryModel.H"
#include "reactingMixture.H"
#include "UniformField.H"
#include "extrapolatedCalculatedFvPatchFields.H"

// i \in \[ 1, ..., NSP-1 \]

#define INERTINDEX this->thermo().composition().species().find("N2")
#define INDEXINGKODES2FOAM(i) ((i) <= (INERTINDEX) ? (i-1) : (i))

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ReactionThermo, class ThermoType>
Foam::kodesChemistryModel<ReactionThermo, ThermoType>::kodesChemistryModel
(
    ReactionThermo& thermo
)
:
    StandardChemistryModel<ReactionThermo, ThermoType>(thermo),
    HostResource(this->thermo().T().size(), NSP, 1),
    HostResourceCopy(HostResource.numOfSystems(), HostResource.sizeOfSystem(), HostResource.numOfParameters()),
    DeviceResourceHostCopy(HostResource.numOfSystems(), HostResource.sizeOfSystem(), HostResource.numOfParameters()),
    ResourceOperator(&HostResourceCopy, &DeviceResourceHostCopy)
{
    HostResource.vectors[0] = this->thermo().T().data();
    HostResource.parameters[0] = this->thermo().p().data();

    for (int i=1; i < HostResource.sizeOfSystem(); ++i)
    {
        HostResource.vectors[i] = this->Y_[INDEXINGKODES2FOAM(i)].data();
    }

    HostResourceCopy.vectors = (scalar**)malloc(HostResource.sizeOfSystem() * sizeof(scalar*));
    for (label i = 0; i < HostResource.sizeOfSystem(); ++i)
    {
        HostResourceCopy.vectors[i] = (scalar*)malloc(HostResource.numOfSystems() * sizeof(scalar));
    }
    HostResourceCopy.parameters = (scalar**)malloc(HostResource.sizeOfSystem() * sizeof(scalar*));
    for (label i = 0; i < HostResource.numOfParameters(); ++i)
    {
        HostResourceCopy.parameters[i] = (scalar*)malloc(HostResource.numOfSystems() * sizeof(scalar));
    }

    HostpyJacMechnismMemory = (mechanism_memory*)malloc(sizeof(mechanism_memory));
    initialize_gpu_memory(HostResourceCopy.numOfSystems(), &HostpyJacMechnismMemory, &DevicepyJacMechnismMemory);

    DeviceResourcePrt = kodes::SeulexDeviceResources::create(HostResource.numOfSystems(), HostResource.sizeOfSystem(), 1, &DeviceResourceHostCopy);

    pyJacSystemPrt = kodes::pyJacSystem::createGPU(DevicepyJacMechnismMemory);

    Integrator = new kodes::Seulex<kodes::pyJacSystem>(pyJacSystemPrt, DeviceResourcePrt, HostResource.numOfSystems());

    Info<< "\n========================================" << nl
        << "  kodesChemistryModel: CUSTOM MODEL LOADED!" << nl
        << "  Using custom chemistry model by user" << nl
        << "  Number of species = " << this->nSpecie_ << nl
        << "  Number of reactions = " << this->nReaction_ << nl
        << "  Thermophysical type: " << ThermoType::typeName() << nl
        << "  Time: " << this->mesh().time().timeName() << nl
        << "========================================\n" << endl;

    // const word inertSpecie(this->thermo().get<word>("inertSpecie"));
    inertIndex = this->thermo().composition().species().find("N2");
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ReactionThermo, class ThermoType>
Foam::kodesChemistryModel<ReactionThermo, ThermoType>::~kodesChemistryModel()
{
    Info<< "kodesChemistryModel: destructor called" << endl;

    kodes::pyJacSystem::destroyGPU(pyJacSystemPrt);
    kodes::SeulexDeviceResources::destroy(DeviceResourcePrt, &DeviceResourceHostCopy);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ReactionThermo, class ThermoType>
Foam::scalar Foam::kodesChemistryModel<ReactionThermo, ThermoType>::solve
(
    const scalar deltaT
)
{
    // return this->StandardChemistryModel<ReactionThermo, ThermoType>::solve(deltaT);

    BasicChemistryModel<ReactionThermo>::correct();

    if (!this->chemistry_)
    {
        return GREAT;
    }

    HostResourceCopy = HostResource;

    ResourceOperator.cpyHostToDevice();

    stepState step(deltaT);

    Integrator->solve(step);

    ResourceOperator.cpyDeviceToHost();

    tmp<volScalarField> trho(this->thermo().rho());
    const scalarField& rho = trho();
    const scalarField& T = this->thermo().T();

    forAll(rho, celli)
    {
        if (T[celli] > this->Treact_)
        {
            scalar Ysum = 0;
            scalar tmpYicelli;

            for (label i=1; i<HostResource.sizeOfSystem(); i++)
            {
                tmpYicelli = HostResourceCopy.vectors[i][celli];
                Ysum += tmpYicelli;
                this->RR_[INDEXINGKODES2FOAM(i)][celli] = rho[celli] * (tmpYicelli - this->Y_[INDEXINGKODES2FOAM(i)][celli])/deltaT;
            }

            this->RR_[inertIndex][celli] = rho[celli] * ((1-Ysum) - this->Y_[inertIndex][celli])/deltaT;
        } else 
        {
            for (label i = 0; i < HostResource.sizeOfSystem(); i++)
            {
                this->RR_[i][celli] = 0;
            }
        }
    }

    return deltaT;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

