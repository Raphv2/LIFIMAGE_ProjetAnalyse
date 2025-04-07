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


float intersect_sphere( /* parametres de la sphere */ const Point& c, const float r, /* parametres du rayon */ const Point& o, const Vector& d )
{

    Vector co = Vector(c, o);
    float a= dot(d,d);
    float b= 2 * dot(d, co);
    float k= dot(co, co) - (r*r) ;

    float dd = b *b - 4* a*k;

    if ( dd < 0)
    return -1;

    float t1 = (-b + sqrt(dd)) / (2 * a);
    float t2 = (-b - sqrt(dd)) / (2 * a);
    float t = std::min(t1, t2);
    if (t < 0) return -1;
    else return t;
}

float intersect_plan( /* parametres du plan */ const Point& a, const Vector& n, /* parametres du rayon */ const Point& o, const Vector& d )
{
    // intersection avec le rayon o , d
    float t= dot (n , Vector (o , a) ) / dot (n , d);

    return t;

}




int main( )
{
    // cree l'image resultat
    Image image(1024, 512);    // par exemple...
    Point c = Point(0, 0, -3);
    float r = 2;
    Point a = Point(0,-1,0);
    Vector n = Vector(0,1,0);
    Vector l= Vector(-4, 6, 1);
    Color emission = Color(200,0,200);

    for(int py= 0; py < image.height(); py++)
    for(int px= 0; px < image.width(); px++)
    {        
        Point o= Point(0, 0, 0);    // origine
        Point e = Point(px - image.width()/2, py - image.height()/2, -100);
        Vector d= normalize(Vector(o, e));     // direction : extremite - origine
        float i = intersect_sphere( /* parametres de la sphere */ c, r, /* parametres du rayon */ o, d );
        float j = intersect_plan( /* parametres du plan */ a, n, /* parametres du rayon */ o, d );
        Point p= o + i *d;
        Vector n_sphere = normalize(p - c);
        float cos_theta = 0;
        
        if(i >= 0.0f && j < 0){
            cos_theta = std::max(dot(n_sphere, normalize(l)), 0.0f);
            image(px, py)= Red() * emission * cos_theta ;
        }
        else if(j >= 0.0f && i < 0){
            cos_theta = std::max(dot(normalize(n), normalize(l)), 0.0f);
            image(px, py) = Blue() * emission * cos_theta ;
        }
        else if(j >= 0.0f && i >= 0.0f){
            cos_theta = std::max(i < j ? dot(n_sphere, normalize(l)) : dot(normalize(n), normalize(l)), 0.0f);
            image(px, py) = i < j ? Red()  * emission * cos_theta : Blue() * emission * cos_theta  ;
        }

        
    }

    write_image(image, "image.png");
    return 0;
}