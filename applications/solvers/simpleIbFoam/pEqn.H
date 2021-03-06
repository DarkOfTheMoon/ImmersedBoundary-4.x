{
    p.boundaryField().updateCoeffs();

    volScalarField rAU(1.0/UEqn().A());
    volVectorField HbyA("HbyA", U);
    HbyA = rAU*UEqn().H();
    UEqn.clear();

    HbyA.correctBoundaryConditions();

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        faceIbMask*
        (fvc::interpolate(HbyA) & mesh.Sf())
    );

    // Adjust immersed boundary fluxes
    immersedBoundaryAdjustPhi(phiHbyA, U);
    adjustPhi(phiHbyA, U, p);

    // Non-orthogonal pressure corrector loop
    while (simple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(rAU, p) == fvc::div(phiHbyA)
        );

        pEqn.setReference(pRefCell, pRefValue);

        pEqn.boundaryManipulate(p.boundaryField());

        pEqn.solve();

        if (simple.finalNonOrthogonalIter())
        {
            phi = phiHbyA - pEqn.flux();
        }
    }

    #include "immersedBoundaryContinuityErrs.H"

    // Explicitly relax pressure for momentum corrector
    p.relax();

    // Momentum corrector
    U = HbyA - rAU*fvc::grad(p);
    U.correctBoundaryConditions();
    fvOptions.correct(U);
}
