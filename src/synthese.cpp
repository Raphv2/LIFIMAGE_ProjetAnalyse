/***************************************************************************************************
 * @file  synthese.cpp
 * @brief projet de synthèse d'image - ray tracing avec flou,
 *        modèle d'éclairage phong, bvh, statistiques, et code refactoré.
 **************************************************************************************************/

 #include <omp.h>

 #include <chrono>
 #include <cmath>
 #include <fstream>
 #include <iostream>
 #include <limits>
 #include <opencv2/opencv.hpp>
 #include <random>
 #include <vector>
 
 #include "color.h"
 #include "image.h"
 #include "image_io.h"
 #include "mat.h"
 #include "mesh_io.h"
 #include "vec.h"
 
 using namespace std::chrono;
 
 // constantes globales
 const float inf =
     std::numeric_limits<float>::infinity();  // valeur infinie pour les calculs
 const float focale =
     10.0f;  // distance focale de la caméra (plus grand = plan plus éloigné)
 const float aperture = 0.02f;  // taille de l'ouverture pour le flou (plus grand
                                // = flou plus étendu)
 
 // statistiques pour le bvh
 int bvh_node_count = 0;  // nombre total de nœuds dans le bvh
 int bvh_leaf_count = 0;  // nombre de feuilles dans le bvh
 int bvh_max_depth = 0;   // profondeur maximale du bvh
 
 // structure représentant une intersection entre un rayon et un objet
 struct Hit {
    float t = inf;        // distance de l'intersection (par défaut infinie)
    Point p;              // point d'intersection
    Vector n;             // normale au point d'intersection
    Color color;          // couleur de l'objet intersecté
    bool mirror = false;  // indique si l'objet est réfléchissant
    operator bool() const {
       return t < inf;
    }  // permet de tester si une intersection est valide
 };
 
 // structure représentant un triangle dans la scène
 struct Triangle {
    Point a, b, c;  // sommets du triangle
    Color color;    // couleur du triangle
 
    // méthode pour calculer l'intersection entre un rayon et le triangle
    Hit intersect(const Point &o, const Vector &d) const {
       Vector ab = Vector(a, b);
       Vector ac = Vector(a, c);
       Vector n = normalize(cross(ab, ac));  // normale au triangle
 
       float denom = dot(n, d);
       if (std::abs(denom) < 1e-6f) return {};  // rayon parallèle au triangle
 
       float t = dot(n, Vector(o, a)) / denom;
       if (t < 0) return {};  // intersection derrière l'origine
 
       Point p = o + t * d;  // point d'intersection
 
       // vérification si le point est à l'intérieur du triangle
       Vector ap = Vector(a, p);
       Vector bp = Vector(b, p);
       Vector cp = Vector(c, p);
       Vector bc = Vector(b, c);
       Vector ca = Vector(c, a);
 
       if (dot(n, cross(ab, ap)) >= 0 && dot(n, cross(bc, bp)) >= 0 &&
           dot(n, cross(ca, cp)) >= 0) {
          return Hit{t, p, n, color};
       }
 
       return {};  // pas d'intersection valide
    }
 };
 
 // structure représentant un plan infini dans la scène
 struct Plan {
    Point a;               // point d'origine du plan
    Vector n;              // normale au plan
    Color color1, color2;  // couleurs pour un damier
    bool mirror = false;   // indique si le plan est réfléchissant
 
    // méthode pour obtenir la couleur au point donné (damier)
    Color get_color_at(const Point &p) const {
       int tx = int(std::floor(p.x));
       int tz = int(std::floor(p.z));
       return ((tx + tz) % 2 == 0) ? color1 : color2;
    }
 
    // méthode pour calculer l'intersection entre un rayon et le plan
    Hit intersect(const Point &o, const Vector &d) const {
       float denom = dot(n, d);
       if (std::abs(denom) < 1e-6f) return {};  // rayon parallèle au plan
       float t = dot(n, Vector(o, a)) / denom;
       if (t < 0) return {};  // intersection derrière l'origine
       Point p = o + t * d;
       return Hit{t, p, normalize(n), get_color_at(p), mirror};
    }
 };
 
 // structure représentant une sphère dans la scène
 struct Sphere {
    Point c;              // centre de la sphère
    float r;              // rayon de la sphère
    Color color;          // couleur de la sphère
    bool mirror = false;  // indique si la sphère est réfléchissante
 
    // méthode pour calculer l'intersection entre un rayon et la sphère
    Hit intersect(const Point &o, const Vector &d) const {
       Vector co = Vector(c, o);
       float a = dot(d, d);
       float b = 2 * dot(d, co);
       float k = dot(co, co) - r * r;
       float delta = b * b - 4 * a * k;
       if (delta < 0) return {};  // pas d'intersection
       float t1 = (-b - std::sqrt(delta)) / (2 * a);
       float t2 = (-b + std::sqrt(delta)) / (2 * a);
       float t = std::min(t1, t2);
       if (t < 0) return {};  // intersection derrière l'origine
       Point p = o + t * d;
       Vector n = normalize(p - c);
       return Hit{t, p, n, color, mirror};
    }
 };
 
 // structure représentant une boîte englobante alignée sur les axes
 struct AABB {
    Point min = Point(inf, inf, inf);     // coin minimal de la boîte
    Point max = Point(-inf, -inf, -inf);  // coin maximal de la boîte
 
    // méthode pour étendre la boîte pour inclure un point donné
    void expand(const Point &p) {
       min = Point(std::min(min.x, p.x), std::min(min.y, p.y),
                   std::min(min.z, p.z));
       max = Point(std::max(max.x, p.x), std::max(max.y, p.y),
                   std::max(max.z, p.z));
    }
 
    // méthode pour tester si un rayon intersecte la boîte
    bool intersect(const Point &o, const Vector &d, float &tmin,
                   float &tmax) const {
       tmin = 0;
       tmax = inf;
       for (int i = 0; i < 3; i++) {
          float invD = 1.0f / d(i);
          float t0 = (min(i) - o(i)) * invD;
          float t1 = (max(i) - o(i)) * invD;
          if (invD < 0) std::swap(t0, t1);
          tmin = std::max(tmin, t0);
          tmax = std::min(tmax, t1);
          if (tmin > tmax) return false;  // pas d'intersection
       }
       return true;  // intersection valide
    }
 };
 
 struct BVHNode {
    AABB box;                                   // boîte englobante du nœud
    BVHNode *left = nullptr, *right = nullptr;  // enfants du nœud
    std::vector<const Triangle *>
        tris;  // triangles contenus dans le nœud (si feuille)
 
    bool is_leaf() const {
       return !left && !right;
    }  // vérifie si le nœud est une feuille
 };
 
 BVHNode *build_bvh(std::vector<const Triangle *> &tris, int depth = 0) {
    if (tris.empty()) return nullptr;  // aucun triangle, retourne un nœud vide
    bvh_node_count++;                  // incrémente le compteur de nœuds
    bvh_max_depth =
        std::max(bvh_max_depth, depth);  // met à jour la profondeur maximale
    BVHNode *node = new BVHNode();
 
    // calcule la boîte englobante pour tous les triangles
    for (auto t : tris)
       node->box.expand(t->a), node->box.expand(t->b), node->box.expand(t->c);
 
    // si le nombre de triangles est faible ou la profondeur trop grande, crée
    // une feuille
    if (tris.size() <= 4 || depth > 16)
       return bvh_leaf_count++, node->tris = tris, node;
 
    // divise les triangles en deux groupes selon un axe
    int axis = depth % 3;  // alterne entre x, y, z
    std::sort(tris.begin(), tris.end(),
              [axis](const Triangle *a, const Triangle *b) {
                 float ca = (a->a(axis) + a->b(axis) + a->c(axis)) / 3.0f;
                 float cb = (b->a(axis) + b->b(axis) + b->c(axis)) / 3.0f;
                 return ca < cb;
              });
 
    size_t mid = tris.size() / 2;
    std::vector<const Triangle *> left(tris.begin(), tris.begin() + mid);
    std::vector<const Triangle *> right(tris.begin() + mid, tris.end());
 
    // construit récursivement les enfants gauche et droit
    node->left = build_bvh(left, depth + 1);
    node->right = build_bvh(right, depth + 1);
    return node;
 }
 
 void delete_bvh(BVHNode *node) {
    if (!node) return;        // si le nœud est nul, rien à supprimer
    delete_bvh(node->left);   // supprime récursivement l'enfant gauche
    delete_bvh(node->right);  // supprime récursivement l'enfant droit
    delete node;              // supprime le nœud courant
 }
 
 Hit intersect_bvh(const BVHNode *node, const Point &o, const Vector &d) {
    float tmin, tmax;
    if (!node || !node->box.intersect(o, d, tmin, tmax))
       return {};  // pas d'intersection avec la boîte
 
    if (node->is_leaf()) {
       Hit best;
       for (auto t : node->tris)  // teste tous les triangles de la feuille
          if (auto h = t->intersect(o, d); h && h.t < best.t) best = h;
       return best;
    }
 
    // teste récursivement les enfants gauche et droit
    auto h1 = intersect_bvh(node->left, o, d);
    auto h2 = intersect_bvh(node->right, o, d);
    return (h1.t < h2.t) ? h1 : h2;  // retourne l'intersection la plus proche
 }
 
 struct Scene {
    std::vector<Sphere> spheres;      // liste des sphères dans la scène
    std::vector<Triangle> triangles;  // liste des triangles dans la scène
    Plan sol;                         // plan représentant le sol
    BVHNode *bvh = nullptr;           // BVH pour accélérer les intersections
 
    // méthode pour tester l'intersection d'un rayon avec tous les objets de la
    // scène
    Hit intersect(const Point &o, const Vector &d) const {
       Hit best;
       for (const Sphere &s : spheres) {  // teste les sphères
          if (Hit h = s.intersect(o, d); h && h.t < best.t) best = h;
       }
       
       if (Hit h = intersect_bvh(bvh, o, d); h && h.t < best.t)
          best = h;  // teste le BVH
          
       if (Hit h = sol.intersect(o, d); h && h.t < best.t)
          best = h;  // teste le sol
       return best;
    }
 };
 
 Vector reflect(const Vector &d, const Vector &n) {
    return d - 2 * dot(d, n) * n;  // calcule le vecteur réfléchi
 }
 
 Color fr(const Vector &d, const Vector &n, const Vector &l, const Color &color,
          const Color &specular, float alpha) {
    Vector h = normalize(d + l);  // vecteur demi-angle
    float cos_theta_h = dot(n, h);
    float phong_term = (alpha + 8) / (8 * float(M_PI)) *
                       std::pow(cos_theta_h, alpha);  // terme de brillance phong
    return color +
           phong_term *
               specular;  // retourne la couleur avec la brillance ajoutée
 }
 
 Color compute_lighting(const Hit &h, const Scene &scene,
                        const std::vector<Vector> &lights, const Vector &d,
                        const Color &emission) {
    Color total = Black();
    for (const Vector &light : lights) {
       Point origin =
           h.p +
           0.001f *
               light;  // déplace légèrement le point pour éviter l'auto-ombre
       Hit shadow = scene.intersect(origin, light);
       if (shadow && shadow.t < inf) {  // si le point est dans l'ombre
          total += h.color * 0.1f;      // ajoute une faible lumière ambiante
          continue;
       }
 
       float diff = std::max(dot(h.n, light), 0.0f);  // intensité diffuse
       Color diffuse = h.color * emission * diff;
 
       float alpha = 32.0f;           // exposant de brillance
       Color specular = Color(1.0f);  // blanc = diélectrique
       total += fr(-d, h.n, light, diffuse, specular,
                   alpha);  // ajoute la lumière réfléchie
    }
    return total;
 }
 
 Color trace_ray(const Scene &scene, const Point &o, const Vector &d,
                 const std::vector<Vector> &lights, const Color &emission,
                 std::default_random_engine &rng,
                 std::uniform_real_distribution<float> &distrib) {
    Hit h = scene.intersect(o, d);
    if (!h) {
       // si le rayon ne touche rien, retourne la couleur du ciel
       Vector L_moon =
           normalize(Vector(0.2f, 1.0f, -0.3f));  // direction vers la lune
       float angle_max = radians(15.0f);
       float cos_max = std::cos(angle_max);
 
       float cos_theta = dot(normalize(d), L_moon);
 
       if (cos_theta > cos_max) {
          return Color(1.0f, 0.95f, 0.85f) * 5.0f;  // lumière blanche chaude
       } else {
          float t = std::max(0.f, dot(normalize(d), Vector(0, 1, 0)));
          return Color(0.2f, 0.3f, 0.5f) * t;  // bleu ciel
       }
    }
 
    Color result = compute_lighting(h, scene, lights, d, emission);
    if (h.mirror) {  // si l'objet est réfléchissant
       Vector r = reflect(d, h.n);
       Point ro = h.p + 0.001f * r;
       Hit rh = scene.intersect(ro, r);
       if (rh) result += compute_lighting(rh, scene, lights, d, emission) * 0.5f;
    }
    return result;
 }
 
 void render_scene(const Scene &scene, Image &image,
                   const std::vector<Vector> &lights, const Color &emission) {
 #pragma omp parallel for schedule(guided)
    for (int py = 0; py < image.height(); py++) {
       std::random_device seed;
       std::default_random_engine rng(seed() +
                                      py);  // générateur aléatoire par ligne
       std::uniform_real_distribution<float> distrib(0.f, 1.f);
 
       int samples = 4096;  // nombre d'échantillons par pixel
       for (int px = 0; px < image.width(); px++) {
          Color final_color = Black();
          for (int s = 0; s < samples; s++) {
             float dx = distrib(rng);
             float dy = distrib(rng);
 
             float angle = 2 * M_PI * distrib(rng);
             float radius = aperture * std::sqrt(distrib(rng));
             float ox = radius * std::cos(angle);
             float oy = radius * std::sin(angle);
             Point o(ox, oy, 2.0f);
             Point f(px + dx - image.width() / 2, py + dy - image.height() / 2,
                     -focale);
             Vector d = normalize(Vector(o, f));
 
             final_color +=
                 trace_ray(scene, o, d, lights, emission, rng, distrib);
          }
 
          Color corrected =
              srgb(final_color / float(samples), 2.2f);  // correction gamma
          corrected.a = 1.0f;
          image(px, py) = corrected;
       }
    }
 }
 
 // Rotation d’un vecteur (ou point relatif à une origine)
 Vector rotateX(const Vector &v, float angle) {
    float c = std::cos(angle);
    float s = std::sin(angle);
    return Vector(v.x, c * v.y - s * v.z, s * v.y + c * v.z);
 }
 
 void load_scene(Scene &scene) {
    // Sol damier
    scene.sol = Plan{Point(0, -1, 0), Vector(0, 1, 0), White(), Black()};
 
    // Sphères déplacées : l'une à gauche, l'autre à droite du mesh
    scene.spheres.push_back(
        Sphere{Point(-3.5f, -1.0f + 1.0f, 1.5f),  // rayon 1.0f, sol à y=-1
               1.0f, Color(1.0f, 0.0f, 0.75f, 1.0f), true});
    scene.spheres.push_back(Sphere{Point(3.5f, -1.0f + 0.7f, 1.0f), 0.7f,
                                   Color(0.0f, 1.0f, 1.0f, 1.0f), false});
 
    
    // Mesh OBJ (Blender)
    std::vector<Point> meshPoints;
    if (!read_positions("data/synthese/monkey_head.obj", meshPoints)) {
       std::cerr << "Erreur chargement .obj" << std::endl;
       exit(1);
    }
 
    // Calculer la boîte englobante du mesh
    Point pmin = meshPoints[0];
    Point pmax = meshPoints[0];
    for (const Point &p : meshPoints) {
       pmin = min(pmin, p);
       pmax = max(pmax, p);
    }
    Vector center = Vector((pmin + pmax) * 0.5f);  // centre de la boîte du mesh
 
 
    // Calculer la diagonale pour normaliser
    // float scale = (2.0f / length(pmax - pmin)) * 2.0f;
    float scale = 1.0f / length(pmax - pmin);  // réduction plus douce
    float rotation_angle = radians(-90.0f);
    Vector offset(0.0f, 0.0f, 1.2f); // vers le centre visuel
 
    for (size_t i = 0; i + 2 < meshPoints.size(); i += 3) {
      Point a = Point(rotateX(Vector(meshPoints[i]     - center) * scale,
  rotation_angle)) + offset; Point b = Point(rotateX(Vector(meshPoints[i + 1] -
  center) * scale, rotation_angle)) + offset; Point c =
  Point(rotateX(Vector(meshPoints[i + 2] - center) * scale, rotation_angle)) +
  offset;
 
 
      Vector n = normalize(cross(Vector(a, b), Vector(a, c)));
      if (dot(n, Vector(0, 1, 0)) < 0) std::swap(b, c);
 
      scene.triangles.push_back(Triangle{a, b, c, Color(0.2f, 0.6f, 0.4f)});
  }
 
 
    // BVH
    std::vector<const Triangle *> ptrs;
    for (const Triangle &t : scene.triangles) ptrs.push_back(&t);
    scene.bvh = build_bvh(ptrs);
    
 }
 
 // Fonction pour convertir une image gKit en Mat OpenCV
 cv::Mat image_to_mat(const Image &img) {
    cv::Mat mat(img.height(), img.width(), CV_8UC3);
    for (int y = 0; y < img.height(); ++y) {
       for (int x = 0; x < img.width(); ++x) {
          Color c = img(x, y);
          cv::Vec3b color;
          color[2] =
              static_cast<unsigned char>(255 * std::clamp(c.r, 0.f, 1.f));  // R
          color[1] =
              static_cast<unsigned char>(255 * std::clamp(c.g, 0.f, 1.f));  // G
          color[0] =
              static_cast<unsigned char>(255 * std::clamp(c.b, 0.f, 1.f));  // B
          mat.at<cv::Vec3b>(img.height() - 1 - y, x) = color;
       }
    }
    return mat;
 }
 
 // fonction pour appliquer un flou gaussien
 void apply_gaussian_blur(const cv::Mat &input,
                          const std::string &output_filename) {
    cv::Mat blurred;
    cv::GaussianBlur(input, blurred, cv::Size(7, 7), 1.5);
    cv::imwrite(output_filename, blurred);
 }
 
 // fonction pour appliquer un filtre bilatéral
 void apply_bilateral_filter(const cv::Mat &input,
                             const std::string &output_filename) {
    cv::Mat filtered;
    cv::bilateralFilter(input, filtered, 9, 75, 75);
    cv::imwrite(output_filename, filtered);
 }
 
 int main() {
    auto start = high_resolution_clock::now();
 
    Scene scene;
    load_scene(scene);  // charge la scène avec les objets et le BVH
 
    std::vector<Vector> lights = {
        normalize(Vector(-4, 6, 1)),
        normalize(Vector(4, 6, -1))};  // sources lumineuses
    Color emission = Color(1.2f);      // intensité lumineuse
    Image image(512, 512);             // image de sortie
 
    cv::Mat mat = image_to_mat(image);
    apply_gaussian_blur(mat, "image_gaussian_blur.png");
    apply_bilateral_filter(mat, "image_bilateral_filter.png");
 
    delete_bvh(scene.bvh);  // libère la mémoire du BVH
 
    auto end = high_resolution_clock::now();
    std::cout << "Temps de rendu : "
              << duration_cast<seconds>(end - start).count() << "s\n";
    std::cout << "BVH Stats: nodes=" << bvh_node_count
              << ", leafs=" << bvh_leaf_count << ", depth=" << bvh_max_depth
              << std::endl;
    return 0;
 }
 