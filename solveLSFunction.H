// solve Level-Set function as the re-initialization equation
Info<< "solve the reinitialization equation"     
      << nl << endl;

psi = psi0;

// Re-initialization of level set function
for (int corr=0; corr<int(epsilon.value()/deltaTau.value()); corr++)
{
   psi = psi + psi0/(mag(psi0)+SMALL) * (double(1)-mag(fvc::grad(psi)*dimChange)) * deltaTau;
   psi.correctBoundaryConditions();
}


// update delta function
Info<< "calculating Dirac function..."     
      << nl << endl;
forAll(psi,celli)
{
   if(mag(psi[celli]) > epsilon.value())
      delta[celli] = double(0);
   else
      delta[celli] = double(1.0)/(double(2.0)*epsilon.value())*(double(1.0)+Foam::cos(M_PI*psi[celli]/epsilon.value()));
};

// update Heaviside function
forAll(H,celli)
{
   if(psi[celli] < -epsilon.value())
      H[celli] = double(0);
   else if(epsilon.value() < psi[celli])
      H [celli] = double(1);
   else
      // normal
      // H[celli] = double(1.0)/double(2.0)*(
      //                  double(1.0) 
      //                + psi[celli] / epsilon.value() 
      //                + Foam::sin(M_PI*psi[celli] / (epsilon.value()) * M_PI)
      //             );
      // density-scaled
      H[celli] = double(1.0)/double(2.0)*(
                       double(1.0) 
                     + psi[celli] / epsilon.value() 
                     + (psi[celli] * psi[celli]) / (2.0 * epsilon.value() * epsilon.value())
                     - 1.0 / (4.0 * M_PI * M_PI) * (Foam::cos((2.0 * M_PI * psi[celli])/epsilon.value()) - 1.0)
                     + (epsilon.value() + psi[celli])/(epsilon.value()*M_PI) * Foam::sin(M_PI*psi[celli]/epsilon.value())
                  );
};

