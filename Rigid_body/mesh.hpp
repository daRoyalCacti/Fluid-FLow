//
// Created by jacob on 24/8/21.
//

#ifndef CODE_MESH_HPP
#define CODE_MESH_HPP

#include <utility>
#include <vector>
#include <numeric>
#include "../MyMath/vec3.hpp"
#include "triangle.hpp"
#include <numeric>

//for loading a model
#include <assimp/Importer.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>

struct bounding_box {
    vec3 min, max;
};


struct mesh final{
    template <class T>
    class no_init_alloc : public std::allocator<T> {
    public:
        using std::allocator<T>::allocator;

        template <class U, class... Args> void construct(U*, Args&&...) {}
    };

    std::vector<vec3> vertices;
    std::vector<vec3, no_init_alloc<vec3>> velocities;
    std::vector<unsigned> indices;
    std::vector<double> mass;
    std::vector<vec3> normals;
    vec3 v, w;  //velocity of center of mass and angular velocity about center of mass
    bounding_box bounds;
    vec3 pos_cm;  //the position of the center of mass

    mesh() = delete;

    mesh(const std::vector<vec3> &vertices_, std::vector<unsigned> indices_, std::vector<double> mass_, std::vector<vec3> normals_, vec3 v_, vec3 w_) noexcept :
    vertices(vertices_), indices(std::move(indices_)), mass(std::move(mass_)), velocities{}, normals(std::move(normals_)), v(v_), w(w_),
    pos_cm{ std::inner_product(mass.begin(), mass.end(), vertices.begin(), vec3{}) / std::accumulate(mass.begin(), mass.end(), 0.0)} {
        velocities.resize(vertices_.size());
        //filling the initial velocities based off the cm velocity and angular velocity
        //pos_cm =  std::inner_product(mass.begin(), mass.end(), vertices.begin(), vec3{}) / std::accumulate(mass.begin(), mass.end(), 0.0);
        for (size_t i = 0; i < velocities.size(); i++) {
            velocities[i] = v_ + cross(w_, (vertices_[i] - pos_cm) );
        }
        update_bounding_box();


#ifndef NDEBUG
        if (vertices.size() != mass.size()) {
            std::cerr << "vertices and mass must be the same size\n";
        }
        if (vertices.size() != normals.size()) {
            std::cerr << "vertices and normals must be the same size\n";
        }
        if (indices.size() < vertices.size()) {
            std::cerr << "having less indices that vertices doesn't make sense\n";
        }
#endif
    }


    static void process_mesh(const aiMesh *mesh, std::vector<vec3> &vert, std::vector<unsigned> &ind, std::vector<vec3> &norm) {

        for (unsigned i = 0; i < mesh->mNumVertices; i++) {
            //process vertex positions, normals and texture coordinates
            vert.emplace_back(mesh->mVertices[i].x, mesh->mVertices[i].y, mesh->mVertices[i].z );
            norm.emplace_back( mesh->mNormals[i].x, mesh->mNormals[i].y, mesh->mNormals[i].z );
        }

        for (unsigned i = 0; i < mesh->mNumFaces; i++) {
            //assimp defines eah mesh as having an array of faces
            // - aiProcess_Triangulate means these faces are always triangles
            //Iterate over all the faces and store the face's indices

            const aiFace face = mesh->mFaces[i];
            for (unsigned j = 0; j < face.mNumIndices; j++) {
                ind.push_back(face.mIndices[j]);
            }
        }

    }



    static void process_node(const aiNode *node, const aiScene *scene, std::vector<vec3> &vert, std::vector<unsigned> &ind, std::vector<vec3> &norm) {
        //process all the node's meshes
        for (unsigned i = 0; i < node->mNumMeshes; i++) {
            const aiMesh *mesh = scene->mMeshes[node->mMeshes[i]];
            process_mesh(mesh, vert, ind, norm);
        }

        //process the meshes of all the nodes children
        for (unsigned i = 0; i < node->mNumChildren; i++) {
            process_node(node->mChildren[i], scene, vert, ind, norm);
        }
    }

    explicit mesh(const std::string &file_name, const vec3 &v_cm, const vec3 &w_cm) : v(v_cm), w(w_cm) {
        const unsigned assimp_settings = aiProcess_Triangulate | aiProcess_GenNormals;
        //aiProcess_Triangulate tells assimp to make the model entirely out of triangles
        //aiProcess_GenNormals creates normal vectors for each vertex if one doesn't already exist

        Assimp::Importer importer;  //cannot be const
        const aiScene *scene = importer.ReadFile(file_name, assimp_settings);
        //reading the data in
        if (!scene || scene->mFlags & AI_SCENE_FLAGS_INCOMPLETE || !scene->mRootNode) {
            std::cerr << "Assimp Error:\n\t" << importer.GetErrorString() << std::endl;
        }
        process_node(scene->mRootNode, scene, vertices, indices, normals);

        //swapping x,y, and z around to get the mesh sittting right
        /*for (unsigned i = 0; i < vertices.size(); i++) {
            auto cpy = vertices[i];
            vertices[i].y() = cpy.z();
            vertices[i].z() = cpy.y();

            cpy = normals[i];
            normals[i].y() = cpy.z();
            normals[i].z() = cpy.y();
            normals[i] = unit_vector(normals[i]);
        }*/

        //scaling the data so its fits in a cube of unit length
        update_bounding_box();

        const double lx = bounds.max.x() - bounds.min.x();
        const double ly = bounds.max.y() - bounds.min.y();
        const double lz = bounds.max.z() - bounds.min.z();

        const double div = std::max({lx, ly, lz});
        const auto mid = vec3(bounds.max.x() + bounds.min.x(), bounds.max.y() + bounds.min.y(), bounds.max.z() + bounds.min.z()  )/2;

        //constexpr vec3 middle_pos = {0,0,0};
        const vec3 middle_pos = vec3(lx, ly, lz)/div;

        //dividing the vertex positions by div in order to make the longest dimension 1
        // and shifting it so the middle aligns with middle_pos
        for (auto &vert : vertices) {
            vert = (vert - mid) / div + middle_pos;
            //vert /= div;
        }
        //updating the bounds from this
        bounds.max = (bounds.max - mid) / div  + middle_pos;
        bounds.min = (bounds.min - mid) / div+  middle_pos;
        //bounds.max/=div;
        //bounds.min/=div;


        //setting mass
        mass.resize(vertices.size());
        for (auto &m : mass) {
            m = 1;
        }

        velocities.resize(vertices.size());
        //filling the initial velocities based off the cm velocity and angular velocity
        pos_cm =  std::inner_product(mass.begin(), mass.end(), vertices.begin(), vec3{}) / std::accumulate(mass.begin(), mass.end(), 0.0);
        for (size_t i = 0; i < velocities.size(); i++) {
            normals[i] = unit_vector(normals[i]);
            velocities[i] = v_cm + cross(w_cm, (vertices[i] - pos_cm) );
        }


    }

    [[nodiscard]] vec3 get_domain() const {
        const double lx = bounds.max.x() - bounds.min.x();
        const double ly = bounds.max.y() - bounds.min.y();
        const double lz = bounds.max.z() - bounds.min.z();

        return {2*lx, 2*ly, 2*lz};
    }



    mesh(const mesh &m) : mass(m.mass), vertices(m.vertices), indices(m.indices), normals(m.normals), v(m.v), w(m.w), pos_cm(m.pos_cm), bounds(m.bounds)  {
        velocities.resize(m.velocities.size() );
        for (size_t i = 0; i < m.velocities.size(); i++) {
            velocities[i] = m.velocities[i];
        }
        update_bounding_box();
    }



    void update_bounding_box() {
        double minx=std::numeric_limits<double>::max(), miny=std::numeric_limits<double>::max(), minz=std::numeric_limits<double>::max();
        double maxx=std::numeric_limits<double>::min(), maxy=std::numeric_limits<double>::min(), maxz=std::numeric_limits<double>::min();
        for (const auto &vert : vertices) {
            if (vert.x() > maxx) {maxx = vert.x();}
            if (vert.x() < minx) {minx = vert.x();}
            if (vert.y() > maxy) {maxy = vert.y();}
            if (vert.y() < miny) {miny = vert.y();}
            if (vert.z() > maxz) {maxz = vert.z();}
            if (vert.z() < minz) {minz = vert.z();}
        }
        bounds = bounding_box{ {minx, miny, minz}, {maxx, maxy, maxz} };
    }

    [[nodiscard]] const vec3& get_vertice_index( const unsigned i) const noexcept {
        return vertices[indices[i]];
    }

    [[nodiscard]] const vec3& get_normal_index( const unsigned i) const noexcept {
        return normals[indices[i]];
    }

    [[nodiscard]] const vec3& get_velocity_index( const unsigned i) const noexcept {
        return velocities[indices[i]];
    }


};

#endif //CODE_MESH_HPP
