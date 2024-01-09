// shapes/sphere.cpp*
#include "shapes/sphere.h"
#include "shapes/micky.h"
#include "sampling.h"
#include "paramset.h"
#include "efloat.h"
#include "stats.h"
#include <cstdio> // for print

namespace pbrt {

// Micky Method Definitions
// explicitly initialize the Sphere members in the member initializer list of the Micky constructor. 
// You cannot use assignment inside the constructor body for members without default constructors.
Micky::Micky(const Transform *ObjectToWorld, const Transform *WorldToObject,
           bool reverseOrientation,
           Float radiusHead, Float radiusL, Float radiusR, 
           Point3f &leftPos, Point3f &rightPos,
           Float phiMax): Shape(ObjectToWorld, WorldToObject, reverseOrientation),
            radiusHead(radiusHead),
            radiusL(radiusL),
            radiusR(radiusR),
            leftPos(Point3f(-1.5f, 0.f, 1.f)),
            rightPos(Point3f(1.5f, 0.f, 1.f)),
            phiMax(Radians(Clamp(phiMax, 0, 360))),
            identity(Transform()),
            object2worldL(Translate(Vector3f(leftPos.x, leftPos.y, leftPos.z))),
            world2objectL(Inverse(object2worldL)),
            object2worldR(Translate(Vector3f(rightPos.x, rightPos.y, rightPos.z))),
            world2objectR(Inverse(object2worldR)),
            sphereHead(&identity, &identity, false, radiusHead, -radiusHead, radiusHead, phiMax),  // 360.f
            sphereLeft(&object2worldL, &world2objectL, false, radiusL, -radiusL, radiusL, phiMax),
            sphereRight(&object2worldR, &world2objectR, false, radiusR, -radiusR, radiusR, phiMax){}
                

Bounds3f Micky::ObjectBound() const {
    return Bounds3f(Point3f(std::min({-radiusHead, leftPos.x-radiusL}), std::min({-radiusHead, leftPos.y-radiusL, rightPos.y-radiusR}), std::min({-radiusHead, leftPos.z-radiusL, rightPos.z-radiusR})),
                    Point3f(std::max({radiusHead, rightPos.x+radiusR}), std::max({ radiusHead, leftPos.y+radiusL, rightPos.y+radiusR}), std::max({ radiusHead, leftPos.z+radiusL, rightPos.z+radiusR})));
}

bool Micky::Intersect(const Ray &r, Float *tHit, SurfaceInteraction *isect,
                       bool testAlphaTexture) const {
    // Jingyi: Start profiling a specific phase, e.g., ShapeIntersect
    ProfilePhase p(Prof::ShapeIntersect);
    // Jingyi: The destructor of ProfilePhase automatically records the time spent in the phase
    
    // get intersections information from three spheres
    SurfaceInteraction isectHead, isectL, isectR;
    Float tHitHead, tHitL, tHitR;
    bool hitHead = sphereHead.Intersect(r, &tHitHead, &isectHead, testAlphaTexture);
    bool hitLeft = sphereLeft.Intersect(r, &tHitL, &isectL, testAlphaTexture);
    bool hitRight = sphereRight.Intersect(r, &tHitR, &isectR, testAlphaTexture);
    // std::printf("tHitHead: %f, tHitL: %f, tHitR: %f\n", tHitHead, tHitL, tHitR);
    std::printf("hitHead: %d, hitLeft: %d, hitRight: %d\n", hitHead, hitLeft, hitRight);

    if(hitHead){
        std::printf("hitHead!!! %f %f %f \n", tHitHead, isectHead.uv.x, isectHead.uv.y);
        tHit = &tHitHead;
        isect = &isectHead;
        if (hitLeft){
            std::printf("hitLeft!!! %f %f %f \n", tHitL, isectL.uv.x, isectL.uv.y);
            if(tHitL < *tHit) {
                tHit = &tHitL;
                isect = &isectL;
            }
        }
        if (hitRight){
            std::printf("hitRight!!! %f %f %f \n", tHitR, isectR.uv.x, isectR.uv.y);
            if(tHitR < *tHit){
                tHit = &tHitR;
                isect = &isectR;
            }
        }
        return true;
    }
    else if (hitLeft){
        std::printf("hitLeft2!!! %f %f %f \n", tHitL, isectL.uv.x, isectL.uv.y);
        tHit = &tHitL;
        isect = &isectL;
        if (hitRight){
            std::printf("hitRight2!!! %f %f %f \n", tHitR, isectR.uv.x, isectR.uv.y);
            if(tHitR < *tHit){
                tHit = &tHitR;
                isect = &isectR;
            }
        }
        return true;
    }
    else if (hitRight){
        std::printf("hitRight3!!! %f %f %f \n", tHitR, isectR.uv.x, isectR.uv.y);
        tHit = &tHitR;
        isect = &isectR;
        return true;
    }
    
    return false;
    

    // // Ignore any alpha textures used for trimming the shape when performing
    // // this intersection. Hack for the "San Miguel" scene, where this is used
    // // to make an invisible area light.
    // if (!Intersect(ray, &tHit, &isectLight, false)) return 0;
}

bool Micky::IntersectP(const Ray &r, bool testAlphaTexture) const {
    ProfilePhase p(Prof::ShapeIntersectP);
    
    if ( sphereHead.IntersectP(r, testAlphaTexture) ||
         sphereLeft.IntersectP(r, testAlphaTexture) ||
         sphereRight.IntersectP(r, testAlphaTexture)
        ) 
        return true;

    return false;

    // Float phi;
    // Point3f pHit;
    // // Transform _Ray_ to object space
    // Vector3f oErr, dErr;
    // Ray ray = (*WorldToObject)(r, &oErr, &dErr);

    // // Compute quadratic sphere coefficients

    // // Initialize _EFloat_ ray coordinate values
    // EFloat ox(ray.o.x, oErr.x), oy(ray.o.y, oErr.y), oz(ray.o.z, oErr.z);
    // EFloat dx(ray.d.x, dErr.x), dy(ray.d.y, dErr.y), dz(ray.d.z, dErr.z);
    // EFloat a = dx * dx + dy * dy + dz * dz;
    // EFloat b = 2 * (dx * ox + dy * oy + dz * oz);
    // EFloat c = ox * ox + oy * oy + oz * oz - EFloat(radius) * EFloat(radius);

    // // Solve quadratic equation for _t_ values
    // EFloat t0, t1;
    // if (!Quadratic(a, b, c, &t0, &t1)) return false;

    // // Check quadric shape _t0_ and _t1_ for nearest intersection
    // if (t0.UpperBound() > ray.tMax || t1.LowerBound() <= 0) return false;
    // EFloat tShapeHit = t0;
    // if (tShapeHit.LowerBound() <= 0) {
    //     tShapeHit = t1;
    //     if (tShapeHit.UpperBound() > ray.tMax) return false;
    // }

    // // Compute sphere hit position and $\phi$
    // pHit = ray((Float)tShapeHit);

    // // Refine sphere intersection point
    // pHit *= radius / Distance(pHit, Point3f(0, 0, 0));
    // if (pHit.x == 0 && pHit.y == 0) pHit.x = 1e-5f * radius;
    // phi = std::atan2(pHit.y, pHit.x);
    // if (phi < 0) phi += 2 * Pi;

    // // Test sphere intersection against clipping parameters
    // if ((zMin > -radius && pHit.z < zMin) || (zMax < radius && pHit.z > zMax) ||
    //     phi > phiMax) {
    //     if (tShapeHit == t1) return false;
    //     if (t1.UpperBound() > ray.tMax) return false;
    //     tShapeHit = t1;
    //     // Compute sphere hit position and $\phi$
    //     pHit = ray((Float)tShapeHit);

    //     // Refine sphere intersection point
    //     pHit *= radius / Distance(pHit, Point3f(0, 0, 0));
    //     if (pHit.x == 0 && pHit.y == 0) pHit.x = 1e-5f * radius;
    //     phi = std::atan2(pHit.y, pHit.x);
    //     if (phi < 0) phi += 2 * Pi;
    //     if ((zMin > -radius && pHit.z < zMin) ||
    //         (zMax < radius && pHit.z > zMax) || phi > phiMax)
    //         return false;
    // }
    // return true;
}

Float Micky::Area() const {
    return sphereHead.Area() + sphereLeft.Area() + sphereRight.Area();
}

Interaction Micky::Sample(const Point2f &u, Float *pdf) const {
    // Point3f pObj = Point3f(0, 0, 0) + radiusHead * UniformSampleSphere(u);
    // Interaction it;
    // it.n = Normalize((*ObjectToWorld)(Normal3f(pObj.x, pObj.y, pObj.z)));
    // if (reverseOrientation) it.n *= -1;
    // // Reproject _pObj_ to sphere surface and compute _pObjError_
    // pObj *= radiusHead / Distance(pObj, Point3f(0, 0, 0));
    // Vector3f pObjError = gamma(5) * Abs((Vector3f)pObj);
    // it.p = (*ObjectToWorld)(pObj, pObjError, &it.pError);
    // *pdf = 1 / Area();
    // return it;
    LOG(FATAL) << "Cone::Sample not implemented.";
    return Interaction();
}

// Interaction Micky::Sample(const Interaction &ref, const Point2f &u,
//                            Float *pdf) const {
//     Point3f pCenter = (*ObjectToWorld)(Point3f(0, 0, 0));

//     // Sample uniformly on sphere if $\pt{}$ is inside it
//     Point3f pOrigin =
//         OffsetRayOrigin(ref.p, ref.pError, ref.n, pCenter - ref.p);
//     if (DistanceSquared(pOrigin, pCenter) <= radiusHead * radiusHead) {
//         Interaction intr = Sample(u, pdf);
//         Vector3f wi = intr.p - ref.p;
//         if (wi.LengthSquared() == 0)
//             *pdf = 0;
//         else {
//             // Convert from area measure returned by Sample() call above to
//             // solid angle measure.
//             wi = Normalize(wi);
//             *pdf *= DistanceSquared(ref.p, intr.p) / AbsDot(intr.n, -wi);
//         }
//         if (std::isinf(*pdf)) *pdf = 0.f;
//         return intr;
//     }

//     // Sample sphere uniformly inside subtended cone

//     // Compute coordinate system for sphere sampling
//     Float dc = Distance(ref.p, pCenter);
//     Float invDc = 1 / dc;
//     Vector3f wc = (pCenter - ref.p) * invDc;
//     Vector3f wcX, wcY;
//     CoordinateSystem(wc, &wcX, &wcY);

//     // Compute $\theta$ and $\phi$ values for sample in cone
//     Float sinThetaMax = radiusHead * invDc;
//     Float sinThetaMax2 = sinThetaMax * sinThetaMax;
//     Float invSinThetaMax = 1 / sinThetaMax;
//     Float cosThetaMax = std::sqrt(std::max((Float)0.f, 1 - sinThetaMax2));

//     Float cosTheta  = (cosThetaMax - 1) * u[0] + 1;
//     Float sinTheta2 = 1 - cosTheta * cosTheta;

//     if (sinThetaMax2 < 0.00068523f /* sin^2(1.5 deg) */) {
//         /* Fall back to a Taylor series expansion for small angles, where
//            the standard approach suffers from severe cancellation errors */
//         sinTheta2 = sinThetaMax2 * u[0];
//         cosTheta = std::sqrt(1 - sinTheta2);
//     }

//     // Compute angle $\alpha$ from center of sphere to sampled point on surface
//     Float cosAlpha = sinTheta2 * invSinThetaMax +
//         cosTheta * std::sqrt(std::max((Float)0.f, 1.f - sinTheta2 * invSinThetaMax * invSinThetaMax));
//     Float sinAlpha = std::sqrt(std::max((Float)0.f, 1.f - cosAlpha*cosAlpha));
//     Float phi = u[1] * 2 * Pi;

//     // Compute surface normal and sampled point on sphere
//     Vector3f nWorld =
//         SphericalDirection(sinAlpha, cosAlpha, phi, -wcX, -wcY, -wc);
//     Point3f pWorld = pCenter + radiusHead * Point3f(nWorld.x, nWorld.y, nWorld.z);

//     // Return _Interaction_ for sampled point on sphere
//     Interaction it;
//     it.p = pWorld;
//     it.pError = gamma(5) * Abs((Vector3f)pWorld);
//     it.n = Normal3f(nWorld);
//     if (reverseOrientation) it.n *= -1;

//     // Uniform cone PDF.
//     *pdf = 1 / (2 * Pi * (1 - cosThetaMax));

//     return it;
// }

// Float Micky::Pdf(const Interaction &ref, const Vector3f &wi) const {
//     Point3f pCenter = (*ObjectToWorld)(Point3f(0, 0, 0));
//     // Return uniform PDF if point is inside sphere
//     Point3f pOrigin =
//         OffsetRayOrigin(ref.p, ref.pError, ref.n, pCenter - ref.p);
//     if (DistanceSquared(pOrigin, pCenter) <= radiusHead * radiusHead)
//         return Shape::Pdf(ref, wi);

//     // Compute general sphere PDF
//     Float sinThetaMax2 = radiusHead * radiusHead / DistanceSquared(ref.p, pCenter);
//     Float cosThetaMax = std::sqrt(std::max((Float)0, 1 - sinThetaMax2));
//     return UniformConePdf(cosThetaMax);
// }

// Float Micky::SolidAngle(const Point3f &p, int nSamples) const {
//     Point3f pCenter = (*ObjectToWorld)(Point3f(0, 0, 0));
//     if (DistanceSquared(p, pCenter) <= radius * radius)
//         return 4 * Pi;
//     Float sinTheta2 = radius * radius / DistanceSquared(p, pCenter);
//     Float cosTheta = std::sqrt(std::max((Float)0, 1 - sinTheta2));
//     return (2 * Pi * (1 - cosTheta));
// }

std::shared_ptr<Shape> CreateMickyShape(const Transform *o2w,
                                         const Transform *w2o,
                                         bool reverseOrientation,
                                         const ParamSet &params) {
    Float radiusHead = params.FindOneFloat("hradius", 1.f);
    Float radiusL = params.FindOneFloat("rradius", .5f);
    Float radiusR = params.FindOneFloat("lradius", .5f);
    Point3f leftPos = params.FindOnePoint3f("pleft", Point3f(-1.5f, 0.f, 1.f));
    Point3f rightPos = params.FindOnePoint3f("pright", Point3f(1.5f, 0.f, 1.f));
    Float phimax = params.FindOneFloat("phimax", 360.f);
    return std::make_shared<Micky>(o2w, w2o, reverseOrientation, 
                                    radiusHead, radiusL, radiusR,
                                    leftPos, rightPos, phimax);
}

}  // namespace pbrt
