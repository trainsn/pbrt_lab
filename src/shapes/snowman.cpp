#include "shapes/sphere.h"
#include "shapes/snowman.h"
#include "sampling.h"
#include "paramset.h"
#include "efloat.h"
#include "stats.h"
// #include <iostream> // for print

namespace pbrt {

// SnowMan Method Definitions
SnowMan::SnowMan(const Transform *ObjectToWorld, const Transform *WorldToObject,
           bool reverseOrientation,
           Float radiusH, Float radiusB, 
           Point3f &posH, Point3f &posB, 
           Float phiMax): Shape(ObjectToWorld, WorldToObject, reverseOrientation),
            radiusH(radiusH),
            radiusB(radiusB),
            posH(posH),
            posB(posB),
            phiMax(Radians(Clamp(phiMax, 0, 360))),
            object2worldHead((*ObjectToWorld) * Translate(Vector3f(posH.x, posH.y, posH.z))),
            world2objectHead((*WorldToObject) * Inverse(Translate(Vector3f(posH.x, posH.y, posH.z)))),
            object2worldBody((*ObjectToWorld) * Translate(Vector3f(posB.x, posB.y, posB.z))),
            world2objectBody((*WorldToObject) * Inverse(Translate(Vector3f(posB.x, posB.y, posB.z)))),
            sphereHead(&object2worldHead, &world2objectHead, false, radiusH, -radiusH, radiusH, phiMax),  // 360.f
            sphereBody(&object2worldBody, &world2objectBody, false, radiusB, -radiusB, radiusB, phiMax){}
            

Bounds3f SnowMan::ObjectBound() const {
    // wrong: return Union(sphereHead.ObjectBound(), sphereBody.ObjectBound());
    Point3f pmin, pmax;
    pmin.x = std::min(posH.x - radiusH, posB.x - radiusB);
    pmax.x = std::max(posH.x + radiusH, posB.x + radiusB);
    pmin.y = std::min(posH.y - radiusH, posB.y - radiusB);
    pmax.y = std::max(posH.y + radiusH, posB.y + radiusB);
    pmin.z = std::min(posH.z - radiusH, posB.z - radiusB);
    pmax.z = std::max(posH.z + radiusH, posB.z + radiusB);

    return Bounds3f(pmin, pmax);
}


bool SnowMan::Intersect(const Ray &r, Float *tHit, SurfaceInteraction *isect,
                       bool testAlphaTexture) const {
    ProfilePhase p(Prof::ShapeIntersect);
    
    // get intersections information from two spheres
    SurfaceInteraction isect1, isect2;
    Float tHit1, tHit2;
    bool hit1 = sphereHead.Intersect(r, &tHit1, &isect1, testAlphaTexture);
    bool hit2 = sphereBody.Intersect(r, &tHit2, &isect2, testAlphaTexture);
    
    // union
    if (hit1) {
        *isect = isect1;
        *tHit = tHit1;
        if (hit2){
            if (tHit2 < tHit1){
                *isect = isect2;
                *tHit = tHit2;
            }
        }
        return true;
    }
    else if (hit2) {
        *isect = isect2;
        *tHit = tHit2;
        return true;
    }
    return false;
}


bool SnowMan::IntersectP(const Ray &r, bool testAlphaTexture) const {
    ProfilePhase p(Prof::ShapeIntersectP);
    
    if ( sphereHead.IntersectP(r, testAlphaTexture) ||
         sphereBody.IntersectP(r, testAlphaTexture)
        ) 
        return true;

    return false;
}


Float SnowMan::Area() const {
    return sphereHead.Area() + sphereBody.Area();
}

Interaction SnowMan::Sample(const Point2f &u, Float *pdf) const {
    LOG(FATAL) << "SnowMan::Sample not implemented.";
    return Interaction();
}


std::shared_ptr<Shape> CreateSnowManShape(const Transform *o2w,
                                         const Transform *w2o,
                                         bool reverseOrientation,
                                         const ParamSet &params) {
    Float radiusH = params.FindOneFloat("radiusH", 1.f);
    Float radiusB = params.FindOneFloat("radiusB", .5f);
    Point3f posH = params.FindOnePoint3f("posH", Point3f(0, 0, 0));
    Point3f posB = params.FindOnePoint3f("posB", Point3f(0, 0, 1.1));
    Float phimax = params.FindOneFloat("phimax", 360.f);
    return std::make_shared<SnowMan>(o2w, w2o, reverseOrientation, 
                                    radiusH, radiusB,
                                    posH, posB, phimax);
}

}  // namespace pbrt
