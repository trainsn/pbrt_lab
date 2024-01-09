#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_SHAPES_MICKY_H
#define PBRT_SHAPES_MICKY_H

// shapes/micky.h*
#include "shape.h"
#include "sphere.h"
#include "../core/transform.h"

namespace pbrt {

// Micky Declarations
class Micky : public Shape {
  public:
    // Micky Public Methods
    Micky(const Transform *ObjectToWorld, const Transform *WorldToObject, bool reverseOrientation,
           Float radiusHead, Float radiusL, Float radiusR, 
           Point3f &leftPos, Point3f &rightPos, Float phiMax);
    Bounds3f ObjectBound() const;
    bool Intersect(const Ray &ray, Float *tHit, SurfaceInteraction *isect,
                   bool testAlphaTexture) const;
    bool IntersectP(const Ray &ray, bool testAlphaTexture) const;
    Float Area() const;
    Interaction Sample(const Point2f &u, Float *pdf) const;
    // Interaction Sample(const Interaction &ref, const Point2f &u, Float *pdf) const;
    // Float Pdf(const Interaction &ref, const Vector3f &wi) const;
    // Float SolidAngle(const Point3f &p, int nSamples) const;

  private:
    // Micky Private Data
    const Float radiusHead, radiusL, radiusR;
    const Float phiMax;
    const Point3f leftPos, rightPos;

    const Transform identity;
    const Transform object2worldL, object2worldR;
    const Transform world2objectL, world2objectR;
    Sphere sphereHead, sphereLeft, sphereRight;
};


std::shared_ptr<Shape> CreateMickyShape(const Transform *o2w,
                                         const Transform *w2o,
                                         bool reverseOrientation,
                                         const ParamSet &params);



}  // namespace pbrt

#endif  // PBRT_SHAPES_MICKY_H
