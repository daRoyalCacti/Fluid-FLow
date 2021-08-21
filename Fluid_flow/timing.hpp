//
// Created by jacob on 21/8/21.
//

#ifndef CODE_TIMING_HPP
#define CODE_TIMING_HPP

#include <chrono>

struct flow_timer{
    std::chrono::time_point<std::chrono::high_resolution_clock> start{};
    std::chrono::time_point<std::chrono::high_resolution_clock> end{};
    double s_create_time{};
    double b_create_time{};
    double p_solve_time{};
    double vx_solve_time{};
    double vy_solve_time{};
    double vz_solve_time{};
    std::ofstream output;

    flow_timer() = delete;
    explicit flow_timer(const char* output_loc) {
        output.open(output_loc);
    }

    ~flow_timer() {
        output.close();
    }

    void save_s_create_time() { s_create_time = static_cast<std::chrono::duration<double>>(end - start).count();}
    void save_b_create_time() {b_create_time = static_cast<std::chrono::duration<double>>(end - start).count();}
    void save_p_solve_time() {p_solve_time = static_cast<std::chrono::duration<double>>(end - start).count();}
    void save_vx_solve_time() {vx_solve_time = static_cast<std::chrono::duration<double>>(end - start).count();}
    void save_vy_solve_time() {vy_solve_time = static_cast<std::chrono::duration<double>>(end - start).count();}
    void save_vz_solve_time() {vz_solve_time = static_cast<std::chrono::duration<double>>(end - start).count();}

    void write_times(const double t) {
        output << t << " " << s_create_time << " " << p_solve_time << " " << b_create_time << " " << vx_solve_time
             << " " << vy_solve_time << " " << vz_solve_time << " " << total_time() << "\n";
        output.flush();
    }


    void set_start(std::chrono::time_point<std::chrono::high_resolution_clock> s) noexcept {start = s;}
    void set_end(std::chrono::time_point<std::chrono::high_resolution_clock> e) noexcept {end = e;}

private:
    double total_time() const noexcept{
        return s_create_time + p_solve_time + b_create_time + vx_solve_time + vy_solve_time + vz_solve_time;
    }
};

#endif //CODE_TIMING_HPP
