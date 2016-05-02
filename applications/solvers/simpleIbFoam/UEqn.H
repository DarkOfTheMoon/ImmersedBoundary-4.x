    // Momentum predictor

    tmp<fvVectorMatrix> UEqn
    (
        fvm::div(phi, U)
      + turbulence->divDevReff(U)
      ==
        fvOptions(U)
    );

    UEqn().relax();

    fvOptions.constrain(UEqn());

    UEqn().boundaryManipulate(U.boundaryField());

    solve(UEqn() == -cellIbMask*fvc::grad(p));

    fvOptions.correct(U);
