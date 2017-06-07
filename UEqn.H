    // Momentum predictor

    tmp<fvVectorMatrix> tUEqn
    (
	fvm::ddt(rho, U)
        fvm::laplacian(mu, U)
    );

    fvVectorMatrix& UEqn = tUEqn.ref();

    UEqn.relax();

    if (simple.momentumPredictor())
    {
        solve(UEqn == fvc::grad(p));
    }
