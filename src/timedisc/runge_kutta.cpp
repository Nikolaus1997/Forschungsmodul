 #include <runge_kutta.h>

LowStorageRungeKutta3::LowStorageRungeKutta3(double dt) : 
dt_(dt) {

}

// Advance the solution using the LSRK3 scheme
void LowStorageRungeKutta3::advance(Variable u, Variable ut) {
    // Coefficients for LSRK3
    const std::vector<double> a = {0.0, -5.0 / 9.0, -153.0 / 128.0};
    const std::vector<double> b = {1.0 / 3.0, 15.0 / 16.0, 8.0 / 15.0};

    // Temporary storage for the right-hand side
    std::vector<double> rhs(u.size()[0], 0.0);

    // Start with the initial solution
    Variable uTemp = u;

    for (int stage = 0; stage < 3; ++stage) {
        // Compute the right-hand side
        rhs = rhsFunc(uTemp);

        // Update the solution in-place
        for (size_t i = 0; i < u.size()[0]; ++i) {
            uTemp(i) = u(i) + a[stage] * (uTemp(i) - u(i)) + dt_ * b[stage] * rhs[i];
        }
    }

    // Update the solution after all stages
    for (int i = 0; i < u.size()[0]; i++)
    {
        u(i) = uTemp(i);
    }

}


