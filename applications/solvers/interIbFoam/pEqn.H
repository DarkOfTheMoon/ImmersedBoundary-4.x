{
#if 0 //TODO fix me
    if (nOuterCorr != 1)
    {
        pd.storePrevIter();
    }
#endif

    volScalarField rUA = 1.0/UEqn.A();
    surfaceScalarField rUAf = fvc::interpolate(rUA);

    U = rUA*UEqn.H();
    // Immersed boundary update
    U.correctBoundaryConditions();

    surfaceScalarField phiU
    (
        "phiU",
        faceIbMask*(fvc::interpolate(U) & mesh.Sf())
    );

    // Adjust immersed boundary fluxes
    immersedBoundaryAdjustPhi(phiU, U);
    adjustPhi(phiU, U, pd);

    phi = phiU
      + faceIbMask*
        (
            fvc::interpolate(interface.sigmaK())*fvc::snGrad(alpha1)
          - ghf*fvc::snGrad(rho)
        )*rUAf*mesh.magSf();


    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pdEqn
        (
            fvm::laplacian(rUAf, pd) == fvc::div(phi)
        );

        pdEqn.setReference(pdRefCell, pdRefValue);

        pdEqn.boundaryManipulate(pd.boundaryField());

        pdEqn.solve(mesh.solver(pd.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            phi -= pdEqn.flux();
        }
    }

    pd.relax();

    U += rUA*fvc::reconstruct((phi - phiU)/rUAf);
    U.correctBoundaryConditions();
}