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

    eqnResidual = solve
    (
        UEqn() == -cellIbMask*fvc::grad(p)
    ).initialResidual();

    maxResidual = max(eqnResidual, maxResidual);

    fvOptions.correct(U);
