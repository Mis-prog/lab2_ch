

double adaptive_integration(double delta, double start, double end, const int n, double* xi, const std::function<double(double)> f) {
    double h = (end - start) / 2.;
    double S = 0.;

    xi[0] = start;

    for (int i = 1; i < n; ++i) {
        xi[i] = xi[i - 1] + h;
        while (true) {
            double Itr = h / 2 * (f(xi[i - 1]) + f(xi[i]));
            double Itr_sost = h / 4 * (f(xi[i - 1]) + 2 * f((xi[i - 1] + xi[i]) / 2.) + f(xi[i]));
            double eps = 1. / 3. * (Itr_sost - Itr);


            if (abs(eps) > h * delta / (end - start)) {
                h /= 2.;
                xi[i] = xi[i - 1] + h;
            }
            else {
                xi[i] = xi[i - 1] + h;
                S += Itr_sost;
                break;
            }
        }
        if (xi[i] + h > end) {
            h = end - xi[i];
        }
    }

    return S;
}

void task4_main(){
    double* adaptive_grid_xi;
    double adaptive_value = 0;
    int adaptive_n = 2;
    {
        // murlib::adaptive_integration();
        for (;;adaptive_n++) {
            double* xi = new double[adaptive_n];
            adaptive_value = adaptive_integration(delta, a, b, adaptive_n, xi, individual_func);
            if (abs(adaptive_value - true_value) < delta) {
                adaptive_grid_xi = new double[adaptive_n];
                std::copy(xi, xi + adaptive_n, adaptive_grid_xi);
                break;
            }
            delete[] xi;
        }
    }
    adaptive_n--;

    std::cout << std::endl << "Adaptive n: " << adaptive_n <<" value " << adaptive_value << " with error " << abs(adaptive_value - true_value) << " (" << delta << ")" << std::endl;
    double* adaptive_grid_yi = new double[adaptive_n];
    for (int i = 0; i < adaptive_n; ++i) {
        adapt
}