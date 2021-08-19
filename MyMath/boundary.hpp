//
// Created by jacob on 19/8/21.
//

#ifndef CODE_BOUNDARY_HPP
#define CODE_BOUNDARY_HPP

#include <vector>

template <class T>
class no_init_alloc
        : public std::allocator<T>
{
public:
    using std::allocator<T>::allocator;

    template <class U, class... Args> void construct(U*, Args&&...) {}
};

//only have to set those that are false
struct edge_point {
    bool has_left = true;
    bool has_right = true;
    bool has_front = true;
    bool has_back = true;
    bool has_up = true;
    bool has_down = true;

    edge_point() = default;
};


template <unsigned N, unsigned M, unsigned P>
struct boundary_points {
    std::vector<edge_point, no_init_alloc<edge_point>> v;

    boundary_points() {
        v.resize( (N+1)*(M+1)*(P+1) );
    }

    [[nodiscard]] edge_point& operator()(const unsigned i, const unsigned j, const unsigned k) { return v[get_index(i,j,k)]; }
    [[nodiscard]] const edge_point& operator()(const unsigned i, const unsigned j, const unsigned k) const { return v[get_index(i,j,k)]; }

    [[nodiscard]] constexpr inline bool has_left(const unsigned i, const unsigned j, const unsigned k) const {return v[get_index(i,j,k)].has_left;}
    [[nodiscard]] constexpr inline bool has_right(const unsigned i, const unsigned j, const unsigned k) const {return v[get_index(i,j,k)].has_right;}
    [[nodiscard]] constexpr inline bool has_down(const unsigned i, const unsigned j, const unsigned k) const {return v[get_index(i,j,k)].has_down;}
    [[nodiscard]] constexpr inline bool has_up(const unsigned i, const unsigned j, const unsigned k) const {return v[get_index(i,j,k)].has_up;}
    [[nodiscard]] constexpr inline bool has_front(const unsigned i, const unsigned j, const unsigned k) const {return v[get_index(i,j,k)].has_front;}
    [[nodiscard]] constexpr inline bool has_back(const unsigned i, const unsigned j, const unsigned k) const {return v[get_index(i,j,k)].has_back;}


    [[nodiscard]] constexpr inline bool has_2left(const unsigned i, const unsigned j, const unsigned k) const {return v[get_index(i,j,k)].has_left && v[get_index(i-1,j,k)].has_left;}
    [[nodiscard]] constexpr inline bool has_2right(const unsigned i, const unsigned j, const unsigned k) const {return v[get_index(i,j,k)].has_right && v[get_index(i+1,j,k)].has_right;}
    [[nodiscard]] constexpr inline bool has_2down(const unsigned i, const unsigned j, const unsigned k) const {return v[get_index(i,j,k)].has_down && v[get_index(i,j-1,k)].has_down;}
    [[nodiscard]] constexpr inline bool has_2up(const unsigned i, const unsigned j, const unsigned k) const {return v[get_index(i,j,k)].has_up && v[get_index(i,j+1,k)].has_up;}
    [[nodiscard]] constexpr inline bool has_2front(const unsigned i, const unsigned j, const unsigned k) const {return v[get_index(i,j,k)].has_front && v[get_index(i,j,k-1)].has_front;}
    [[nodiscard]] constexpr inline bool has_2back(const unsigned i, const unsigned j, const unsigned k) const {return v[get_index(i,j,k)].has_back && v[get_index(i,j,k+1)].has_back;}

private:
    [[nodiscard]] constexpr inline unsigned get_index(const unsigned i, const unsigned j, const unsigned k) const {
        return i + (N+1)*j + (N+1)*(M+1)*k;
    }

};

#endif //CODE_BOUNDARY_HPP
