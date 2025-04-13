/***************************************************************************************************
 * @file  synthese.cpp
 * @brief Contains the main program of the project
 **************************************************************************************************/

 #include <iostream>
 #include <stdexcept>
 #include "vec.h"
 #include <cmath>
 #include "color.h"
 #include "image.h"
 #include "image_io.h"
 #include <random>
 #include <mesh_io.h>
 #include <limits>
 #include <omp.h>
 #include <chrono>
 using namespace std::chrono;
 #include "mat.h" // pour radians()
 
 // constantes globales
 const float inf = std::numeric_limits<float>::infinity();
 const float focale = 4.0f;      // distance de mise au point
 const float aperture = 0.02f;    // ouverture de la lentille (flou de profondeur)

// structure de base pour stocker une intersection
struct Hit
{
    float t;        // distance depuis l'origine du rayon
    Point p;        // point d'intersection
    Vector n;       // normale au point
    Color color;    // couleur de l'objet
    bool mirror = false; // objet miroir ?


    Hit() : t(inf) {} // aucun hit / pas d'intersection par défaut
    Hit(float _t, const Point& _p, const Vector& _n, const Color& _color, bool _mirror)
        : t(_t), p(_p), n(_n), color(_color), mirror(_mirror) {}

    operator bool() const { return t < inf; } // permet d'utiliser : if(hit) ...
};

// structure Trianglez avec intersection rayon/triangle
struct Triangle {
    Point a, b, c;
    Color color;

    Hit intersect(const Point& o, const Vector& d) const {
        Vector ab = Vector(a, b);
        Vector ac = Vector(a, c);
        Vector n = normalize(cross(ab, ac)); // normale du triangle

        float denom = dot(n, d);
        if (std::abs(denom) < 1e-6f) return Hit(); // rayon // au triangle

        float t = dot(n, Vector(o, a)) / denom;
        if (t < 0) return Hit(); // intersection derrière la caméra

        Point p = o + t * d;

        // test de barycentrie
        Vector ap = Vector(a, p);
        float area = length(cross(ab, ac));
        float u = length(cross(Vector(p, b), Vector(p, c))) / area;
        float v = length(cross(Vector(p, c), Vector(p, a))) / area;
        float w = 1 - u - v;

        if (u < 0 || v < 0 || w < 0) return Hit(); // dehors du triangle

        return Hit(t, p, n, color, false);
    }
};

// structure Plan (damier)
struct Plan
{
    Point a;
    Vector n;
    Color color1;
    Color color2;
    bool mirror = false;

    Color get_color_at(const Point& p) const
    {
        // coordonnées x et z du point d'intersection
        int tx = int(std::floor(p.x));
        int tz = int(std::floor(p.z));
        if ((tx + tz) % 2 == 0)
            return color1;
        else
            return color2;
    }

    Hit intersect(const Point& o, const Vector& d) const
    {
        float denom = dot(n, d);
        if (std::abs(denom) < 1e-6f)
            return Hit();

        float t = dot(n, Vector(o, a)) / denom;
        if (t < 0)
            return Hit();

        Point p = o + t * d;
        return Hit(t, p, normalize(n), get_color_at(p), mirror);
    }
};

// structure Sphere classique
struct Sphere
{
    Point c;
    float r;
    Color color;
    bool mirror = false;

    Hit intersect(const Point& o, const Vector& d) const
    {
        Vector co = Vector(c, o);
        float a = dot(d, d);
        float b = 2 * dot(d, co);
        float k = dot(co, co) - r * r;
        float delta = b * b - 4 * a * k;

        if (delta < 0) return Hit();

        float t1 = (-b - sqrt(delta)) / (2 * a);
        float t2 = (-b + sqrt(delta)) / (2 * a);
        float t = std::min(t1, t2);

        if (t < 0) return Hit();

        Point p = o + t * d;
        Vector n = normalize(p - c);
        //return Hit(t, p, n, texture(p), mirror);
        return Hit(t, p, n, color, mirror);

    }
};

// boite englobante axis-alignée
struct AABB {
    Point min, max;

    AABB() : min(inf, inf, inf), max(-inf, -inf, -inf) {}

    void expand(const Point& p) {
        min = Point(std::min(min.x, p.x), std::min(min.y, p.y), std::min(min.z, p.z));
        max = Point(std::max(max.x, p.x), std::max(max.y, p.y), std::max(max.z, p.z));
    }

    bool intersect(const Point& o, const Vector& d, float& tmin, float& tmax) const {
        tmin = 0;
        tmax = inf;
        for (int i = 0; i < 3; i++) {
            float invD = 1.0f / d(i);
            float t0 = (min(i) - o(i)) * invD;
            float t1 = (max(i) - o(i)) * invD;
            if (invD < 0) std::swap(t0, t1);
            tmin = std::max(tmin, t0);
            tmax = std::min(tmax, t1);
            if (tmin > tmax) return false;
        }
        return true;
    }
};

// noeud BVH
struct BVHNode {
    AABB box;
    BVHNode* left = nullptr;
    BVHNode* right = nullptr;
    std::vector<const Triangle*> tris;

    bool is_leaf() const { return !left && !right; }
};

//construction récursive du BVH
BVHNode* build_bvh(std::vector<const Triangle*>& tris, int depth = 0) {
    if (tris.empty()) return nullptr;

    BVHNode* node = new BVHNode();

    for (const Triangle* t : tris) {
        node->box.expand(t->a);
        node->box.expand(t->b);
        node->box.expand(t->c);
    }

    if (tris.size() <= 4 || depth > 16) {
        node->tris = tris;
        return node;
    }

    int axis = depth % 3;
    std::sort(tris.begin(), tris.end(), [axis](const Triangle* a, const Triangle* b) {
        float ca = (a->a(axis) + a->b(axis) + a->c(axis)) / 3.0f;
        float cb = (b->a(axis) + b->b(axis) + b->c(axis)) / 3.0f;
        return ca < cb;
    });

    size_t mid = tris.size() / 2;
    std::vector<const Triangle*> left(tris.begin(), tris.begin() + mid);
    std::vector<const Triangle*> right(tris.begin() + mid, tris.end());

    node->left = build_bvh(left, depth + 1);
    node->right = build_bvh(right, depth + 1);

    return node;
}

// intersection BVH récursive
Hit intersect_bvh(const BVHNode* node, const Point& o, const Vector& d) {
    float tmin, tmax;
    if (!node || !node->box.intersect(o, d, tmin, tmax)) return Hit();

    if (node->is_leaf()) {
        Hit best;
        for (const Triangle* t : node->tris) {
            Hit h = t->intersect(o, d);
            if (h && h.t < best.t) best = h;
        }
        return best;
    }

    Hit left_hit = intersect_bvh(node->left, o, d);
    Hit right_hit = intersect_bvh(node->right, o, d);
    return (left_hit.t < right_hit.t) ? left_hit : right_hit;
}

// structure Scene avec tous les objets
struct Scene
{
    std::vector<Sphere> spheres;
    std::vector<Triangle> triangles;
    BVHNode* bvh = nullptr;

    Plan sol;

    Hit intersect(const Point& o, const Vector& d) const
    {
        Hit best;
        for (const Sphere& s : spheres)
        {
            Hit h = s.intersect(o, d);
            if (h && h.t < best.t) best = h;
        }

        Hit h_tri = intersect_bvh(bvh, o, d);
        if (h_tri && h_tri.t < best.t) best = h_tri;
                   

        Hit h_sol = sol.intersect(o, d);
        if (h_sol && h_sol.t < best.t) best = h_sol;

        return best;
    }
};

// calcul de réflexion
Vector reflect(const Vector& d, const Vector& n)
{
    return d - 2 * dot(d, n) * n;
}

int main()
{
    auto start = high_resolution_clock::now();

    // Crée la scène
    Scene scene;
    // sphère principale
    //scene.spheres.push_back(Sphere{ Point(0, 0, -3), 2, Color(1.0f, 0.0f, 0.75f, 1.0f), true }); // magenta néon
    // sphère 2
    //scene.spheres.push_back(Sphere{ Point(-2.5f, -0.75f, -1.0f), 0.5f, Color(0.0f, 1.0f, 1.0f, 1.0f) }); // cyan électrique

    scene.sol = Plan{
        Point(0, -1, 0),
        Vector(0, 1, 0),
        Color(1.0f, 1.0f, 1.0f, 1.0f), // blanc
        //Color(0.0f, 0.0f, 0.0f, 1.0f)  // noir
        Color(1.0f, 1.0f, 1.0f, 1.0f), // blanc
    };  
    
    // chargement de l'objet blender 'bea_skeleton_psx'
    /*
    MeshIOData mesh;
    if (!read_meshio_data("data/synthese/bea_skeleton_psx.obj", mesh)) {
        std::cerr << "Erreur chargement objet" << std::endl;
        return 1;
    }
    */

    // chargement de l'objet blender 'cube'
    
    MeshIOData mesh;
    if (!read_meshio_data("data/synthese/cube.obj", mesh)) {
        std::cerr << "Erreur chargement objet" << std::endl;
        return 1;
    }
    

    // afficher les coordonnées de la mesh
    std::cout << "Exemple de coordonnées mesh :\n";
    for (int i = 0; i < 3; ++i)
        std::cout << mesh.positions[i] << std::endl;
        
    // convertir les triangles de 'bea_skeleton_psx' en mesh
    /*
    for (size_t i = 0; i + 2 < mesh.indices.size(); i += 3) {
        Point a = mesh.positions[mesh.indices[i]] + Vector(0, 9.5f, 0);
        Point b = mesh.positions[mesh.indices[i + 1]] + Vector(0, 9.5f, 0);
        Point c = mesh.positions[mesh.indices[i + 2]] + Vector(0, 9.5f, 0);
        int mat_id = mesh.material_indices[i / 3];
        Color col = mesh.materials(mat_id).diffuse;
        scene.triangles.push_back(Triangle{a, b, c, col});
    }
    */

    // convertir les triangles de 'cube' en mesh

    for (size_t i = 0; i + 2 < mesh.indices.size(); i += 3) {
        Point a = mesh.positions[mesh.indices[i]] + Vector(0, -2.0f, 0);
        Point b = mesh.positions[mesh.indices[i + 1]] + Vector(0, -2.0f, 0);
        Point c = mesh.positions[mesh.indices[i + 2]] + Vector(0, -2.0f, 0);        
        int mat_id = mesh.material_indices[i / 3];
        Color col = mesh.materials(mat_id).diffuse;
        scene.triangles.push_back(Triangle{a, b, c, col});
    }
    

    // construire le BVH
    std::vector<const Triangle*> triangle_ptrs;
    for (const Triangle& t : scene.triangles) {
        triangle_ptrs.push_back(&t);
    }
    scene.bvh = build_bvh(triangle_ptrs);
        
     

    // Deux lumières directionnelles
    std::vector<Vector> lumières = {
        normalize(Vector(-4, 6, 1)),   // gauche
        normalize(Vector(4, 6, -1))    // droite
    };

    Color emission = Color(1.2f, 1.2f, 1.2f, 1.0f); // lumière blanche douce et forte
    // rendu
    //Image image(1024, 512);
    // dev
    Image image(512, 512);

    #pragma omp parallel for schedule(dynamic, 1)
    for (int py = 0; py < image.height(); py++)
    {
        std::random_device seed;
        std::default_random_engine rng(seed());
        std::uniform_real_distribution<float> distrib(0.f, 1.f);

        for (int px = 0; px < image.width(); px++)
        {
            // rendu
            int samples = 128;
            // dev
            //int samples = 8;
            Color final_color = Black();

            for (int s = 0; s < samples; s++)
            {
                float dx = distrib(rng);
                float dy = distrib(rng);

                // position de la lentille (œil décalé dans un disque)
                /*
                Point o = Point(0, 0, 0);
                Point e = Point(px + dx - image.width() / 2, py + dy - image.height() / 2, -100);
                Vector d = normalize(Vector(o, e));
                */

                // position de la lentille (œil décalé dans un disque)
                float angle = 2 * M_PI * distrib(rng);
                float radius = aperture * sqrt(distrib(rng));
                float ox = radius * std::cos(angle);
                float oy = radius * std::sin(angle);
                Point o = Point(ox, oy, 20.0f); // position de la lentille
                // point sur le plan focal
                Point f = Point(px + dx - image.width() / 2, py + dy - image.height() / 2, -focale);
                Vector d = normalize(Vector(o, f));
                
                Hit h = scene.intersect(o, d);
                if (!h)
                {
                    // lumière dôme directionnelle
                    // rendu
                    int nb_directions = 512;
                    Color ciel = Black();
                    Color ciel_base = Color(0.7f, 0.85f, 1.0f, 1.0f); // bleu ciel très clair
                    
                    Vector L = normalize(Vector(0, 1, 0)); // direction principale (zenith)
                    float angle_max = radians(30.f);
                    float cos_max = std::cos(angle_max);              
                    for (int k = 0; k < nb_directions; ++k)
                    {
                        float x = 2 * distrib(rng) - 1;
                        float y = distrib(rng);
                        float z = 2 * distrib(rng) - 1;
                        Vector l = normalize(Vector(x, y, z));
                    
                        float cos_theta = dot(l, L);
                        if (cos_theta > cos_max){
                            ciel += ciel_base * cos_theta;
                        }               
                    }                   

                    ciel = ciel / float(nb_directions);
                    final_color += ciel * 10.0f;
                    //final_color += Black(); // fond noir = aucune lumière du ciel
                    continue;
                }

                Color total_light = Black();

                for (const Vector& light_dir : lumières)
                {
                    Vector to_light = light_dir;
                    Point shadow_origin = h.p + 0.001f * to_light;
                    Hit shadow_hit = scene.intersect(shadow_origin, to_light);

                    bool shadow = shadow_hit && shadow_hit.t < inf;

                    if (!shadow)
                    {
                        float cos_theta = std::max(dot(h.n, to_light), 0.0f);
                        total_light += h.color * emission * cos_theta;
                    }
                    else
                    {
                        total_light += h.color * 0.1f;
                    }
                }

                //final_color += total_light;

                final_color += total_light * 1.0f;

                if (h.mirror)
                {
                    Vector r = reflect(d, h.n);
                    Point ro = h.p + 0.001f * r;
                    Hit rh = scene.intersect(ro, r);
                
                    if (rh)
                    {
                        for (const Vector& light_dir : lumières)
                        {
                            Vector to_light = light_dir;
                            Point shadow_origin = rh.p + 0.001f * to_light;
                            Hit shadow_hit = scene.intersect(shadow_origin, to_light);
                            bool shadow = shadow_hit && shadow_hit.t < inf;
                
                            if (!shadow)
                            {
                                float cos_theta = std::max(dot(rh.n, to_light), 0.0f);
                                total_light += rh.color * emission * cos_theta;
                            }
                            else
                            {
                                total_light += rh.color * 0.1f;
                            }
                        }
                
                        final_color += total_light * 0.5f; // réflexion atténuée
                    }
                }
                                
            }

            Color corrected = srgb(final_color / float(samples), 2.2f);
            corrected.a = 1.0f;
            image(px, py) = corrected;
        }
    }

    write_image(image, "image.png");

    auto end = high_resolution_clock::now();
    auto duration = duration_cast<seconds>(end - start);
    std::cout << "temps de rendu : " << duration.count() << " secondes" << std::endl;

    return 0;
}



