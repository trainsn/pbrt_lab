/*  lab 1: a basic ray tracer
    more advanced features like reflections, refractions, and global illumination 
    will be added later for a complete physically based rendering system.

    object(): intersect
    material(): scatter
    render(scene, width, height)
        for i , for j:
            ray = Ray(origin, direction)
            image[i, j] = trace_ray(ray, scene, depth=0)
        return image

    intersect_scene(ray, scene):
        for obj in scene:
            t = obj.intersect(ray)
            hit_point = ray.origin + t * ray.direction
            normal = (hit_point - obj.center) / obj.radius
            material = Lambertian(albedo=np.array([0.8, 0.8, 0.8]))  # Example Lambertian material
        return hit_point, normal, material

    trace_ray(ray, scene, depth=0):
        hit_point, normal, material = intersect_scene(ray, scene)
        if hit_point is not None:
            scattered_ray, attenuation = material.scatter(ray, hit_point, normal)
            return attenuation * trace_ray(scattered_ray, scene, depth + 1)
        else:
            return np.array([0.5, 0.7, 1.0]) # Background color
    
*/

#include <iostream>

#include "core/api.h"
#include "core/pbrt.h"

#include "core/paramset.h"
#include "core/parser.cpp"
#include "core/film.h"
#include "core/light.h"
#include "core/transform.h"
#include "core/medium.h"
#include "core/spectrum.h" // SpectrumType
#include "core/scene.h"

#include "filters/gaussian.h"
#include "filters/box.h"
#include "lights/spot.h"
#include "lights/diffuse.h"
#include "lights/point.h"
#include "cameras/perspective.h"
#include "shapes/sphere.h"
#include "shapes/cylinder.h"
#include "shapes/triangle.h"
#include "materials/matte.h"
#include "materials/plastic.h"
#include "textures/constant.h"
#include "textures/checkerboard.h"
#include "samplers/halton.h"
#include "integrators/path.h"
#include "accelerators/bvh.h"

using namespace std;
using namespace pbrt;

int main () {
    printf("lab 1 \n\n");

    ParallelInit();  // Threads must be launched before the profiler is initialized.
    InitProfiler();

    /* ---------------------------- create image ---------------------------- */
    ParamSet filterParams;
    // std::unique_ptr<Filter> filter(CreateGaussianFilter(filterParams));
    std::unique_ptr<Filter> filter(new GaussianFilter(Vector2f(1., 1.), 1.0));

    // std::unique_ptr<std::string[]> strFileName(new std::string[1]);
    // strFileName[0] = std::string("test_lab1.png");
    // std::unique_ptr<int[]> xres(new int[1]{300});
    // std::unique_ptr<int[]> yres(new int[1]{300});

    // ParamSet filmParams;
    // filmParams.AddInt("xresolution", std::move(xres), 1);
    // filmParams.AddInt("yresolution", std::move(yres), 1);
    // filmParams.AddString("filename", std::move(strFileName), 1);
    // Film *film = CreateFilm(filmParams, std::move(filter));
    Point2i resolution(300, 300);
    std::string filename("test_lab1.png");
    Film *film = new Film(resolution, Bounds2f(Point2f(0, 0), Point2f(1, 1)), // crop, Point2f(1, 1) --> resolution
                        std::move(filter), 1., filename, 1.);
    


    /* ---------------------------- create camera ---------------------------- */
    // Medium *medium1; // TBD: define medium
    float transformStart = 0;
    float transformEnd = 1;
    Transform transform = Transform(); 
    Transform lookAt = LookAt(Point3f(5.f, 0.f, 0.f), Point3f(0.f, 0.f, 0.f), Vector3f(0.f, 0.f, 1.f)); // LookAt(pos, look, up)
    transform = transform * lookAt; // worldToCamera
    transform = Inverse(transform); // cameraToWorld
    AnimatedTransform cam2World(&transform, transformStart, &transform, transformEnd);
    AnimatedTransform identity(new Transform, 0, new Transform, 1);
    
    // std::unique_ptr<Float[]> camFov(new Float[1]{50.f});
    // ParamSet cameraParams;
    // cameraParams.AddFloat("fov", std::move(camFov), 1);
    // PerspectiveCamera *cam = CreatePerspectiveCamera(cameraParams, identity, film, medium1);
    // std::shared_ptr<Camera> camera = std::make_shared<PerspectiveCamera>(*cam);
    std::shared_ptr<Camera> camera = std::make_shared<PerspectiveCamera>(
                    cam2World, Bounds2f(Point2f(-1, -1), Point2f(1, 1)),  // screenWindow
                    0., 1.,  // shutterOpen, shutterClose
                    0., 10., 50, film, nullptr); 
                    // lensRadius, focalDistance, fov, medium
                    // lensRadius = 0 # default to pinhole camera



    /* ---------------------------- create light ---------------------------- */
    Medium *medium2; 
    Float colorI[3] = {50.f, 50.f, 50.f},
            from[3] = {1.5f, 1.0f, 4.1f},
              to[3] = {0.f, 0.f, 0.f};
    Transform light2world = Transform();
    // Spectrum spectrumI = Spectrum::FromRGB(colorI);
    std::unique_ptr<Float[]> spectrumI(new Float[3]);
    std::copy(std::begin(colorI), std::end(colorI), spectrumI.get());
    std::unique_ptr<Point3f[]> fromP(new Point3f[1]{Point3f(from[0], from[1], from[2])});
    std::unique_ptr<Point3f[]> toP(new Point3f[1]{Point3f(to[0], to[1], to[2])});
    std::unique_ptr<Float[]> coneangle(new Float[1]{60.f});
    
    ParamSet lightParams;
    lightParams.AddRGBSpectrum("I", std::move(spectrumI), 3);
    lightParams.AddPoint3f("from", std::move(fromP), 1);
    lightParams.AddPoint3f("to", std::move(toP), 1);
    lightParams.AddFloat("coneangle", std::move(coneangle), 1);
    std::vector<std::shared_ptr<Light>> lights;
    std::shared_ptr<SpotLight> light = CreateSpotLight(light2world, medium2, lightParams);
    lights.push_back(light);



    /* ---------------------------- create checkboard plane ---------------------------- */
    // std::map<std::string, std::shared_ptr<Texture<Float>>> floatTextures;
    // std::map<std::string, std::shared_ptr<Texture<Spectrum>>> spectrumTextures;
    
    // ParamSet geomParams, materialParams;
    // std::unique_ptr<Float[]> su(new Float[1]{8.f});
    // std::unique_ptr<Float[]> sv(new Float[1]{8.f});
    // geomParams.AddFloat("uscale", std::move(su), 1);
    // geomParams.AddFloat("vscale", std::move(sv), 1);

    // std::unique_ptr<Float[]> tex1(new Float[3]{0.1f, 0.1f, 0.1f});
    // std::unique_ptr<Float[]> tex2(new Float[3]{0.8f, 0.8f, 0.8f});
    // geomParams.AddRGBSpectrum("tex1", std::move(tex1), 3);
    // geomParams.AddRGBSpectrum("tex2", std::move(tex2), 3);
    
    std::shared_ptr<Texture<Spectrum>> tex1Spect = std::make_shared<ConstantTexture<Spectrum>>(Spectrum(0.1));
    std::shared_ptr<Texture<Spectrum>> tex2Spect = std::make_shared<ConstantTexture<Spectrum>>(Spectrum(0.8));
    std::shared_ptr<Texture<Float>> sigma  = std::make_shared<ConstantTexture<Float>>(0.);
    

    // Transform tex2world = Translate(Vector3f(0.f, 0.f, -1.f));
    // std::unique_ptr<TextureMapping3D> map(new IdentityMapping3D(tex2world)); // 3D
    std::unique_ptr<TextureMapping2D> map(new UVMapping2D(8.f, 8.f, 0.f, 0.f));
    // TextureParams texParams(geomParams, materialParams, floatTextures, spectrumTextures);
    // Texture<Spectrum> *checkerboardTexture = CreateCheckerboardSpectrumTexture(tex2world, texParams);
    
    std::shared_ptr<Texture<Spectrum>> texturePtr = std::make_shared<Checkerboard2DTexture<Spectrum>>(std::move(map), tex1Spect, tex2Spect, AAMethod::ClosedForm);
    // Checkerboard2DTexture<Spectrum> *checkerboardTexture = new Checkerboard2DTexture<Spectrum>(std::move(map), tex1Spect, tex2Spect, AAMethod::ClosedForm);
    // std::shared_ptr<Texture<Spectrum>> checkerboardTexturePtr(checkerboardTexture);
    std::shared_ptr<Material> checkerMatteMaterial = std::make_shared<MatteMaterial>(texturePtr, sigma, nullptr);
    
    // create trianglemesh:
    ParamSet trimeshParams;
    std::unique_ptr<int[]> indices(new int[6]{0, 1, 2, 0, 2, 3});
    std::unique_ptr<Point3f[]> P(new Point3f[4]{
                    Point3f(-20, -20, 0), Point3f(20, -20, 0), 
                    Point3f(20, 20, 0), Point3f(-20, 20, 0)});
    std::unique_ptr<Point2f[]> st(new Point2f[4]{
                    Point2f(0, 0), Point2f(1, 0), 
                    Point2f(1, 1), Point2f(0, 1)});
    trimeshParams.AddInt("indices", std::move(indices), 6);
    trimeshParams.AddPoint3f("P", std::move(P), 4);
    trimeshParams.AddPoint2f("st", std::move(st), 4);
    std::vector<std::shared_ptr<Shape>> trianglemesh;
    Transform tex2world = Translate(Vector3f(0.f, 0.f, -1.f));
    Transform world2tex = Inverse(tex2world);
    trianglemesh = CreateTriangleMeshShape(&tex2world, &world2tex, false, trimeshParams, nullptr);

    /* ---------------------------- create shapes ---------------------------- */
    // ParamSet sphere1Params;
    // ParamSet sphere2Params;
    // ParamSet cylinder1Params;
    // std::unique_ptr<Float[]> sphere1Rad(new Float[1]{0.5});
    // sphere1Params.AddFloat("radius", std::move(sphere1Rad), 1);
    // std::unique_ptr<Float[]> sphere2Rad(new Float[1]{0.25});
    // sphere2Params.AddFloat("radius", std::move(sphere2Rad), 1);
    // std::unique_ptr<Float[]> cylinder1Rad(new Float[1]{0.2});
    // cylinder1Params.AddFloat("radius", std::move(cylinder1Rad), 1);
    // std::unique_ptr<Float[]> cylinder1zmin(new Float[1]{-0.5});
    // cylinder1Params.AddFloat("zmin", std::move(cylinder1zmin), 1);
    // std::unique_ptr<Float[]> cylinder1zmax(new Float[1]{0.5});
    // cylinder1Params.AddFloat("zmax", std::move(cylinder1zmax), 1);

    bool reverseOrientation = false;
    Transform object2world1 = Transform() * Translate(Vector3f(1.f, 0.f, -1.f));
    Transform object2world2 = object2world1 * Translate(Vector3f(0.f, -1.f, 1.f));
    Transform object2world3 = object2world2 * Translate(Vector3f(0.f, 0.5f, -0.5f));
    Transform world2object1 = Inverse(object2world1);
    Transform world2object2 = Inverse(object2world2);
    Transform world2object3 = Inverse(object2world3);
    // std::cout << "object2world1: " << object2world1 << std::endl;
    // std::cout << "object2world2: " << object2world2 << std::endl;
    // std::cout << "object2world3: " << object2world3 << std::endl;

    std::vector<std::shared_ptr<Shape>> shapes;
    // std::shared_ptr<Shape> s;
    // s = CreateSphereShape(&object2world1, &world2object1, reverseOrientation, sphere1Params);
    // shapes.push_back(s);
    // s = CreateSphereShape(&object2world2, &world2object2, reverseOrientation, sphere2Params);
    // shapes.push_back(s);
    // s = CreateCylinderShape(&object2world3, &world2object3, reverseOrientation, cylinder1Params);
    // shapes.push_back(s);
    
    std::shared_ptr<Shape> sphere1 = std::make_shared<Sphere>(&object2world1, &world2object1, reverseOrientation, 
                            0.5, -1, 1, 360);  /* radius, zMin, zMax, phiMax*/
    std::shared_ptr<Shape> sphere2 = std::make_shared<Sphere>(&object2world2, &world2object2, reverseOrientation, 
                            0.25, -1, 1, 360);  /* radius, zMin, zMax, phiMax*/
    std::shared_ptr<Shape> cylinder1 = std::make_shared<Cylinder>(&object2world3, &world2object3, reverseOrientation, 
                            0.2, -0.5, 0.5, 360);  /* radius, zMin, zMax, phiMax*/
    shapes.push_back(sphere1);
    shapes.push_back(sphere2);
    shapes.push_back(cylinder1);

    
    Float Kd1f[3] = {0.4f, 0.5f, 0.4f},
          Kd2f[3] = {0.9f, 0.8f, 0.8f},
          Kd3f[3] = {0.4f, 0.2f, 0.2f};

    // std::unique_ptr<Float[]> Kd1(new Float[3]{0.4f, 0.5f, 0.4f});
    // std::unique_ptr<Float[]> Ks1(new Float[3]{0.3f, 0.3f, 0.3f});
    // std::unique_ptr<Float[]> Kd2(new Float[3]{0.9f, 0.8f, 0.8f});
    // std::unique_ptr<Float[]> Kd3(new Float[3]{0.4f, 0.2f, 0.2f});
    // std::unique_ptr<Float[]> Ks3(new Float[3]{0.5f, 0.5f, 0.5f});
    // std::unique_ptr<Float[]> roughness(new Float[1]{0.2});

    // ParamSet GeoParams, shapeMatParams1, shapeMatParams2, shapeMatParams3;
    // shapeMatParams1.AddRGBSpectrum("Kd", std::move(Kd1), 3);
    // shapeMatParams1.AddRGBSpectrum("Ks", std::move(Ks1), 3);
    // shapeMatParams1.AddFloat("roughness", std::move(roughness), 1);
    // shapeMatParams2.AddRGBSpectrum("Kd", std::move(Kd2), 3);
    // shapeMatParams3.AddRGBSpectrum("Kd", std::move(Kd3), 3);
    // shapeMatParams3.AddRGBSpectrum("Ks", std::move(Ks3), 3);
    // shapeMatParams3.AddFloat("roughness", std::move(roughness), 1);

    // TextureParams shapeTexParams1(GeoParams, shapeMatParams1, floatTextures, spectrumTextures);
    // TextureParams shapeTexParams2(GeoParams, shapeMatParams2, floatTextures, spectrumTextures);
    // TextureParams shapeTexParams3(GeoParams, shapeMatParams3, floatTextures, spectrumTextures);
    
    std::shared_ptr<Texture<Spectrum>> Kd1 = std::make_shared<ConstantTexture<Spectrum>>(Spectrum::FromRGB(Kd1f));
    std::shared_ptr<Texture<Spectrum>> Ks1 = std::make_shared<ConstantTexture<Spectrum>>(Spectrum(0.3));
    std::shared_ptr<Texture<Spectrum>> Kd2 = std::make_shared<ConstantTexture<Spectrum>>(Spectrum::FromRGB(Kd2f));
    std::shared_ptr<Texture<Spectrum>> Kd3 = std::make_shared<ConstantTexture<Spectrum>>(Spectrum::FromRGB(Kd2f));
    std::shared_ptr<Texture<Spectrum>> Ks3 = std::make_shared<ConstantTexture<Spectrum>>(Spectrum(0.5));
    std::shared_ptr<Texture<Float>> roughness1 = std::make_shared<ConstantTexture<Float>>(0.2);
    std::shared_ptr<Texture<Float>> roughness3 = std::make_shared<ConstantTexture<Float>>(0.025);
    std::vector<std::shared_ptr<Material>> materials;
    // std::shared_ptr<Material> mat;
    // MatteMaterial *mattePtr;
    // PlasticMaterial *plasticPtr1, *plasticPtr2;
    // plasticPtr1 = CreatePlasticMaterial(shapeTexParams1);
    // mattePtr = CreateMatteMaterial(shapeTexParams2);
    // plasticPtr2 = CreatePlasticMaterial(shapeTexParams3);
    // mat = std::make_shared<PlasticMaterial>(*plasticPtr1);
    // materials.push_back(mat);
    // mat = std::make_shared<MatteMaterial>(*mattePtr);
    // materials.push_back(mat);
    // mat = std::make_shared<PlasticMaterial>(*plasticPtr2);
    // materials.push_back(mat);
    // MatteMaterial* matPtr = CreateMatteMaterial(shapeTexParams1);
    // std::shared_ptr<MatteMaterial> materialPtr = std::make_shared<MatteMaterial>(*matPtr);

    std::shared_ptr<Material> material1 = std::make_shared<PlasticMaterial>(Kd1, Ks1, roughness1, nullptr, false);
    std::shared_ptr<Material> material2 = std::make_shared<MatteMaterial>(Kd2, sigma, nullptr);
    std::shared_ptr<Material> material3 = std::make_shared<PlasticMaterial>(Kd3, Ks3, roughness3, nullptr, false);
    materials.push_back(material1);
    materials.push_back(material2);
    materials.push_back(material3);

    // shape to Primitive
    MediumInterface mediumInterface;
    std::vector<std::shared_ptr<Primitive>> prims;
    // prims.reserve(shapes.size());
    // std::cout << shapes.size() << " ,,,, " << trianglemesh.size() << std::endl;
    // prims.reserve(shapes.size() + trianglemesh.size());
    
    // for (auto s : shapes) {
    for (size_t i = 0; i < shapes.size(); ++i){
        // std::shared_ptr<AreaLight> areal;
        std::shared_ptr<Shape> s = shapes[i];
        std::shared_ptr<Material> m = materials[i];
        prims.push_back(
            std::make_shared<GeometricPrimitive>(s, m, nullptr, mediumInterface)
        );
    }

    for (size_t i = 0; i < trianglemesh.size(); ++i){
        std::shared_ptr<Shape> s = trianglemesh[i];
        std::cout << i << " success " << std::endl;
        prims.push_back(
            std::make_shared<GeometricPrimitive>(s, checkerMatteMaterial, nullptr, mediumInterface)
        );
    }

    std::shared_ptr<BVHAccel> bvhPtr = std::make_shared<BVHAccel>(prims);


    // Sampler "halton" "integer pixelsamples" [4] 
    int nsamp = 16;
    std::shared_ptr<HaltonSampler> samplerHalton = std::make_shared<HaltonSampler>(nsamp, film->GetSampleBounds());
    // params.FindOneInt("pixelsamples", 16); // 256
    // CreateHaltonSampler(paramSet, film->GetSampleBounds());

    // Integrator "path" "integer maxdepth" [1]
    int maxDepth = 1;
    std::unique_ptr<Integrator> integrator(new PathIntegrator(maxDepth, camera, samplerHalton, film->croppedPixelBounds));
    // integrator = CreatePathIntegrator(IntegratorParams, sampler, camera);

    std::unique_ptr<Scene> scene(new Scene(bvhPtr, lights));
    
    if (scene && integrator) {
        integrator->Render(*scene);
    }

    ParallelCleanup();
    CleanupProfiler();
    
}

    // SpectrumType spectrumType = SpectrumType::Reflectance; // enum class SpectrumType { Reflectance, Illuminant };

    // Transform t1 = Transform(); //Translate(Vector3f(0.2f, 0.f, 0.f));
    // Transform t2 = Inverse(t1);
    // std::shared_ptr<Shape> spheretmp = std::make_shared<Sphere>(
    //     &t1, &t2, true, 0.5, -1, 1, 360);

    // std::shared_ptr<Texture<Spectrum>> Kdtmp =
    //     std::make_shared<ConstantTexture<Spectrum>>(Spectrum(0.8));
    // std::shared_ptr<Texture<Float>> sigmatmp =
    //     std::make_shared<ConstantTexture<Float>>(0.);
    // std::shared_ptr<Material> materialtmp =
    //     std::make_shared<MatteMaterial>(Kdtmp, sigmatmp, nullptr);

    // // std::shared_ptr<AreaLight> areaLight = std::make_shared<DiffuseAreaLight>(t1, nullptr, Spectrum(0.8), 1, spheretmp);

    // prims.push_back(std::make_shared<GeometricPrimitive>(
    //     spheretmp, materialtmp, nullptr, mediumInterface));

    
    // // // std::shared_ptr<AreaLight> areaLight = std::make_shared<DiffuseAreaLight>(t1, nullptr, Spectrum(0.8), 1, spheretmp);

    // std::vector<std::shared_ptr<Primitive>> primstmp;
    // primstmp.push_back(std::make_shared<GeometricPrimitive>(
    //     spheretmp, materialtmp, nullptr, mediumInterface));

    // std::shared_ptr<BVHAccel> bvhtmp = std::make_shared<BVHAccel>(primstmp);

    // std::unique_ptr<Scene> scenetmp(new Scene(bvhtmp, lights));
    

    // Integrator *integratortmp = new PathIntegrator(8, camera, samplerHalton, film->croppedPixelBounds);

    // if (scenetmp && integratortmp) {
    //     integratortmp->Render(*scenetmp);
    // }