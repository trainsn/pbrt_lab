#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_SHAPES_SNOWMAN_H
#define PBRT_SHAPES_SNOWMAN_H

#include "shape.h"
#include "sphere.h"
#include "../core/transform.h"

namespace pbrt {

// SnowMan Declarations
class SnowMan : public Shape {
  public:
    // SnowMan Public Methods
    SnowMan(const Transform *ObjectToWorld, const Transform *WorldToObject, bool reverseOrientation,
           Float radiusH, Float radiusB, Point3f &posH, Point3f &posB, Float phiMax);
    Bounds3f ObjectBound() const;
    bool Intersect(const Ray &ray, Float *tHit, SurfaceInteraction *isect,
                   bool testAlphaTexture) const;
    bool IntersectP(const Ray &ray, bool testAlphaTexture) const;
    Interaction Sample(const Point2f &u, Float *pdf) const;
    Float Area() const;
    
  private:
    // SnowMan Private Data
    const Float radiusH, radiusB;
    const Float phiMax;
    const Point3f posH, posB;

    const Transform identity;
    const Transform object2worldHead, object2worldBody;
    const Transform world2objectHead, world2objectBody;
    Sphere sphereHead, sphereBody;
};


std::shared_ptr<Shape> CreateSnowManShape(const Transform *o2w,
                                         const Transform *w2o,
                                         bool reverseOrientation,
                                         const ParamSet &params);

}  // namespace pbrt

#endif  // PBRT_SHAPES_SNOWMAN_H
