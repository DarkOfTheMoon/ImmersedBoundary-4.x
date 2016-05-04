{
#if 0 //TODO fix me
    if (nOuterCorr != 1)
    {
        pd.storePrevIter();
    }
#endif

    volScalarField rAU("rAU", 1.0/UEqn.A());
    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));

    U = rAU*UEqn.H();
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
        )*rAUf*mesh.magSf();


    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pdEqn
        (
            fvm::laplacian(rAUf, pd) == fvc::div(phi)
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

    U += rAU*fvc::reconstruct((phi - phiU)/rAUf);
    U.correctBoundaryConditions();
}
