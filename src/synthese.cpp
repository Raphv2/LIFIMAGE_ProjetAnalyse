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

#include <limits>

#include <omp.h>

#include <chrono>
using namespace std::chrono;

#include "mat.h"  // pour avoir accès à radians()



const float inf = std::numeric_limits<float>::infinity();

struct Hit
{
    float t;        // distance au point d'intersection
    Point p;        // point d'intersection
    Vector n;       // normale au point
    Color color;    // couleur de l'objet
    bool mirror = false;


    Hit() : t(inf) {} // pas d'intersection par défaut
    Hit(float _t, const Point& _p, const Vector& _n, const Color& _color, bool _mirror)
        : t(_t), p(_p), n(_n), color(_color), mirror(_mirror) {}

    operator bool() const { return t < inf; } // permet d'écrire : if(hit) ...
};

struct Plan
{
    Point a;
    Vector n;
    Color color;
    bool mirror = false;


    Hit intersect(const Point& o, const Vector& d) const
    {
        float denom = dot(n, d);
        if (std::abs(denom) < 1e-6f)
            return Hit();

        float t = dot(n, Vector(o, a)) / denom;
        if (t < 0)
            return Hit();

        Point p = o + t * d;
        return Hit(t, p, normalize(n), color, mirror);

    }
};

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
        return Hit(t, p, n, color, mirror);

    }
};

struct Scene
{
    std::vector<Sphere> spheres;
    Plan sol;

    Hit intersect(const Point& o, const Vector& d) const
    {
        Hit best;
        for (const Sphere& s : spheres)
        {
            Hit h = s.intersect(o, d);
            if (h && h.t < best.t) best = h;
        }

        Hit h_sol = sol.intersect(o, d);
        if (h_sol && h_sol.t < best.t) best = h_sol;

        return best;
    }
};

/*
float intersect_plan(const Point& a, const Vector& n, const Point& o, const Vector& d)
{
    float denom = dot(n, d);
    if (std::abs(denom) < 1e-6f) // rayon // au plan ?
        return inf;

    float t = dot(n, Vector(o, a)) / denom;

    if (t < 0)
        return inf;

    return t;
}

float intersect_sphere(const Point& c, const float r, const Point& o, const Vector& d)
{
    Vector co = Vector(c, o); // ← direction depuis le centre vers le rayon
    float a = dot(d, d);
    float b = 2 * dot(d, co);
    float k = dot(co, co) - r * r;

    float delta = b * b - 4 * a * k;
    if (delta < 0)
        return inf;

    float t1 = (-b - sqrt(delta)) / (2 * a);
    float t2 = (-b + sqrt(delta)) / (2 * a);
    float t = std::min(t1, t2);

    if (t < 0)
        return inf;

    return t;
}

<<<<<<< Updated upstream

void test_intersections()
=======
int main( )
>>>>>>> Stashed changes
{
    std::cout << "- Test intersection plan -" << std::endl;
    Point o = Point(0, 0, 0);
    Vector d = Vector(0, 0, -1);
    Point a = Point(0, 0, -1);
    Vector n = Vector(0, 0, 1);
    float t1 = intersect_plan(a, n, o, d);
    std::cout << "Plan : t = " << t1 << " (attendu : 1)" << std::endl;

    std::cout << "- Test intersection sphere -" << std::endl;
    Point c = Point(0, 0, -3);
    float r = 2;
    float t2 = intersect_sphere(c, r, o, d);
    std::cout << "Sphere : t = " << t2 << " (attendu : 1)" << std::endl;
}
*/

Vector reflect(const Vector& d, const Vector& n)
{
    return d - 2 * dot(d, n) * n;
}

int main()
{
    auto start = high_resolution_clock::now();

    // Crée la scène
    Scene scene;
    scene.spheres.push_back(Sphere{ Point(0, 0, -3), 2, Color(1.0f, 0.2f, 0.5f, 1.0f), true }); // sphère miroir rose
    scene.spheres.push_back(Sphere{ Point(-2.5f, -0.75f, -1.0f), 0.5f, Color(0.2f, 1.0f, 0.2f, 1.0f) }); // vert fluo


    scene.sol = Plan{ Point(0, -1, 0), Vector(0, 1, 0), Color(0.5f, 0.5f, 0.5f, 1.0f) }; // gris moyen


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
            //int samples = 32;
            // dev
            int samples = 16;
            Color final_color = Black();

            for (int s = 0; s < samples; s++)
            {
                float dx = distrib(rng);
                float dy = distrib(rng);

                Point o = Point(0, 0, 0);
                Point e = Point(px + dx - image.width() / 2, py + dy - image.height() / 2, -100);
                Vector d = normalize(Vector(o, e));

                Hit h = scene.intersect(o, d);
                if (!h)
                {
                    // lumière dôme directionnelle
                    // rendu
                    //int nb_directions = 512;
                    // dev
                    int nb_directions = 128;
                    Color ciel = Black();
                    Color ciel_base = Color(0.7f, 0.85f, 1.0f, 1.0f); // bleu ciel très clair
                    
                    Vector L = normalize(Vector(0, 1, 0)); // direction principale (zenith)
                    float angle_max = radians(30.f);
                    float cos_max = std::cos(angle_max);
                    
                    /*
                    #pragma omp parallel
                    {
                        Color ciel_local = Black();
                        std::default_random_engine local_rng(seed() + omp_get_thread_num() + py * 10000 + px);
                        std::uniform_real_distribution<float> distrib(0.f, 1.f);
                    
                        #pragma omp for nowait
                        for (int k = 0; k < nb_directions; ++k)
                        {
                            float x = 2 * distrib(local_rng) - 1;
                            float y = distrib(local_rng);
                            float z = 2 * distrib(local_rng) - 1;
                            Vector l = normalize(Vector(x, y, z));
                    
                            float cos_theta = dot(l, L);
                            if (cos_theta > cos_max)
                                ciel_local += ciel_base * cos_theta;
                        }
                    
                        #pragma omp critical
                        ciel += ciel_local;
                    }
                    */                 
                   for (int k = 0; k < nb_directions; ++k)
                   {
                       float x = 2 * distrib(rng) - 1;
                       float y = distrib(rng);
                       float z = 2 * distrib(rng) - 1;
                       Vector l = normalize(Vector(x, y, z));
                   
                       float cos_theta = dot(l, L);
                       if (cos_theta > cos_max)
                           ciel += ciel_base * cos_theta;
                   }                   

                    ciel = ciel / float(nb_directions);
                    //final_color += ciel * 10.0f;
                    final_color += Black(); // fond noir = aucune lumière du ciel
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



