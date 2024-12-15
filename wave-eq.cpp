#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

const double dx = 0.001;    // 空間の刻み幅
const double dt = 0.001;   // 時間の刻み幅
const double c = 1.0;      // 波の速度
const int nx = 1000;        // 空間方向の格子点数
const int nt = 1000;       // 時間方向のステップ数

void initialize(std::vector<double>& u) {
    for (int i = 0; i < nx; ++i) {
        double x = i * dx;
        u[i] = sin(M_PI * x); // 初期条件: 正弦波
    }
}

void save_to_file(const std::vector<std::vector<double>>& data, const std::string& filename) {
    std::ofstream file(filename);
    for (const auto& row : data) {
        for (double val : row) {
            file << val << " ";
        }
        file << "\n";
    }
    file.close();
}

void plot_gnuplot(FILE* gp, const std::vector<double>& data, const std::string& label) {
    for (int i = 0; i < data.size(); ++i) {
        fprintf(gp, "%f %f\n", i*dx, data[i]);
    }
    fprintf(gp, "e\n");
}

std::vector<std::vector<double>> matrix_maker(double k){

    std::vector<std::vector<double>> matrix(nx,std::vector<double>(nx,0.0));
    for(int i=0;i<nx;i++){
        matrix[i][i]=2*k+1;

        if(i-1>=0){
            matrix[i][i-1]=-k;
        }
        if(i+1<=nx-1){
            matrix[i][i+1]=-k;
        }
    }
    return matrix;
}

int main() {
    std::vector<double> u(nx, 0.0);       // 現在の時刻の値
    std::vector<double> u_prev(nx, 0.0); // 前の時刻の値
    std::vector<double> u_next(nx, 0.0); // 次の時刻の値

    initialize(u_prev); // 初期条件を設定
    u = u_prev;

    double r = (c * dt / dx) * (c * dt / dx); // 安定条件に基づく係数
    if (r > 1.0) {
        std::cerr << "安定条件を満たしていません (r = " << r << ")\n";
        return -1;
    }

    std::vector<std::vector<double>> results;
    results.push_back(u_prev);

    for (int t = 0; t < nt; ++t) {
        for (int i = 1; i < nx - 1; ++i) {
            u_next[i] = 2 * (1 - r) * u[i] - u_prev[i] + r * (u[i + 1] + u[i - 1]);
        }
        // 境界条件 (固定端)
        u_next[0] = 0.0;
        u_next[nx - 1] = 0.0;

        // 時間を進める
        u_prev = u;
        u = u_next;

        // 保存
        results.push_back(u);
    }

    FILE *gp;

    gp=popen("gnuplot -persist","w");

    // Gnuplotの設定
    fprintf(gp, "set title '2D Data Plot'\n");
    fprintf(gp, "set xlabel 'X'\n");
    fprintf(gp, "set ylabel 'Y'\n");


    int display_step = 10;

    std::string plot_command ;
    plot_command.append("plot ");
    for(int i=0; i<nt; i=i+display_step){
        plot_command.append("'-' with lines, ");
    }
    plot_command.append("\n");

    fprintf(gp,"%s",plot_command.c_str());

    for(int i=0; i < nt; i=i+display_step){
        std::cout <<"plotting "<<std::to_string(i)<<"\n";
        plot_gnuplot(gp,results[i],std::to_string(i));
    }

    fflush(gp);
    pclose(gp);

    std::cout << "シミュレーションが終了しました。\n";

    //save_to_file(results, "wave_simulation.txt");
    //std::cout << "シミュレーションが終了しました。結果は 'wave_simulation.txt' に保存されました。\n";
    return 0;
}
