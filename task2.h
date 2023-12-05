
double parable(int n0) {
    vector<double> x(2 * n0 + 1), y(2 * n0 + 1);
    double sum_rect = 0, h = nodeEqualsFill(x, y, 2 * n0);
    sum_rect = acos_foo(a) + acos_foo(b);
    for (int i = 1; i <= 2 * n0 - 1; i++) {
        if (i % 2 == 0) {
            sum_rect+= 2*acos_foo(x[i]);
        }else{
            sum_rect+= 4*acos_foo(x[i]);
        }
    }
    sum_rect = ((b-a) * sum_rect)/(6.0*n0);
    return abs(answer - sum_rect);
}

void task2_main() {
    int best_parab;
    double parab_cur;
    bool check_parab = 0;
    ofstream out2;
    out2.open("task2.txt");
    for (int i = 0; i < 1e4; i++) {
        parab_cur = parable(i);
        out2 << i << "," << parab_cur << endl;
        if (parab_cur < epsilon and !check_parab) {
            best_parab = i;
            check_parab = 1;
        }
    }
    out2.close();
    cout << "Парабола: " << best_parab << endl;
    system("python D:/5sem/numerical/Lab2/plot_2.py");
}