{
    p.boundaryField().updateCoeffs();

    volScalarField rAU(1.0/UEqn().A());

    U = UEqn().H()*rAU;
    // Immersed boundary update
    U.correctBoundaryConditions();

    UEqn.clear();
    phi = faceIbMask*(fvc::interpolate(U) & mesh.Sf());

    // Adjust immersed boundary fluxes
    immersedBoundaryAdjustPhi(phi, U);
    adjustPhi(phi, U, p);

    // Non-orthogonal pressure corrector loop
    while (simple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(rAU, p) == fvc::div(phi)
        );

        pEqn.setReference(pRefCell, pRefValue);

        pEqn.boundaryManipulate(p.boundaryField());

        // Retain the residual from the first iteration
        if (simple.nNonOrthCorr() == 0)
        {
            eqnResidual = pEqn.solve().initialResidual();
            maxResidual = max(eqnResidual, maxResidual);
        }
        else
        {
            pEqn.solve();
        }

        if (simple.finalNonOrthogonalIter())
        {
            phi -= pEqn.flux();
        }
    }

    #include "immersedBoundaryContinuityErrs.H"

    // Explicitly relax pressure for momentum corrector
    p.relax();

    // Momentum corrector
    U -= rAU*fvc::grad(p);
    U.correctBoundaryConditions();
    fvOptions.correct(U);
}
