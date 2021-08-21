//
// Created by jacob on 21/8/21.
//

#ifndef CODE_TIMING_HPP
#define CODE_TIMING_HPP

#include <chrono>

struct flow_timer final {
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
    explicit flow_timer(const char* output_loc) noexcept {
        output.open(output_loc);
    }

    ~flow_timer() noexcept {
        output.close();
    }

    void save_s_create_time() noexcept { s_create_time = static_cast<std::chrono::duration<double>>(end - start).count();}
    void save_b_create_time() noexcept {b_create_time = static_cast<std::chrono::duration<double>>(end - start).count();}
    void save_p_solve_time() noexcept {p_solve_time = static_cast<std::chrono::duration<double>>(end - start).count();}
    void save_vx_solve_time() noexcept {vx_solve_time = static_cast<std::chrono::duration<double>>(end - start).count();}
    void save_vy_solve_time() noexcept {vy_solve_time = static_cast<std::chrono::duration<double>>(end - start).count();}
    void save_vz_solve_time() noexcept {vz_solve_time = static_cast<std::chrono::duration<double>>(end - start).count();}

    void write_times(const double t) noexcept {
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
