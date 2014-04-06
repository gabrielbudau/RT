#include <math.h>   
#include <stdlib.h> 
#include <stdio.h> 
#include <fstream>
#include <omp.h>

double erand48(unsigned short xsubi[3]){
  return (double) rand() / (double) RAND_MAX;
}

struct Vec { // Usage: time ./smallpt 5000 && xv image.ppm
  double x, y, z; // position, also color (r,g,b)

  Vec(double x_ = 0, double y_ = 0, double z_ = 0) {
    x = x_;
    y = y_;
    z = z_;
  }

  Vec operator+(const Vec &b) const {
    return Vec(x + b.x, y + b.y, z + b.z);
  }

  Vec operator-(const Vec &b) const {
    return Vec(x - b.x, y - b.y, z - b.z);
  }

  Vec operator*(double b) const {
    return Vec(x*b, y*b, z * b);
  }

  Vec mult(const Vec &b) const {
    return Vec(x * b.x, y * b.y, z * b.z);
  }

  Vec& norm() {
    return *this = *this * (1 / sqrt(x * x + y * y + z * z));
  }

  double dot(const Vec &b) const {
    return x * b.x + y * b.y + z * b.z;
  } // cross:

  Vec operator%(Vec&b) {
    return Vec(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x);
  }
};

struct Ray {
  Vec o, d;

  Ray(Vec o_, Vec d_) : o(o_), d(d_) {
  }
};

enum Refl_t {
  DIFF, SPEC, REFR
}; // material types, used in radiance()

struct Sphere {
  double rad; // radius
  Vec p, e, c; // position, emission, color
  Refl_t refl; // reflection type (DIFFuse, SPECular, REFRactive)

  Sphere(double rad_, Vec p_, Vec e_, Vec c_, Refl_t refl_) :
    rad(rad_), p(p_), e(e_), c(c_), refl(refl_) {
  }

  double intersect(const Ray &r) const { // returns distance, 0 if nohit
    Vec op = p - r.o; // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
    double 
      t, 
      eps = 1e-4,
      b = op.dot(r.d), 
      det = b * b - op.dot(op) + rad*rad;

    if (det < 0) 
      return 0;
    else 
      det = sqrt(det);

    return (t = b - det)>eps ? t : ((t = b + det) > eps ? t : 0);
  }
};
Sphere spheres [] = {//Scene: radius, position, emission, color, material
  Sphere(1e5, Vec(1e5 + 1, 40.8, 81.6), Vec(), Vec(.75, .25, .25), DIFF), //Left
  Sphere(1e5, Vec(-1e5 + 99, 40.8, 81.6), Vec(), Vec(.25, .25, .75), DIFF), //Rght
  Sphere(1e5, Vec(50, 40.8, 1e5), Vec(), Vec(.75, .75, .75), DIFF), //Back
  Sphere(1e5, Vec(50, 40.8, -1e5 + 170), Vec(), Vec(), DIFF), //Frnt
  Sphere(1e5, Vec(50, 1e5, 81.6), Vec(), Vec(.75, .75, .75), DIFF), //Botm
  Sphere(1e5, Vec(50, -1e5 + 81.6, 81.6), Vec(), Vec(.75, .75, .75), DIFF), //Top
  Sphere(16.5, Vec(27, 16.5, 47), Vec(), Vec(1, 1, 1)*.999, SPEC), //Mirr
  Sphere(16.5, Vec(73, 16.5, 78), Vec(), Vec(1, 1, 1)*.999, REFR), //Glas
  Sphere(600, Vec(50, 681.6 - .27, 81.6), Vec(12, 12, 12), Vec(), DIFF) //Lite
};

double clamp(double x) {
  return x < 0 ? 0 : x > 1 ? 1 : x;
}

int toInt(double x) {
  return int(pow(clamp(x), 1 / 2.2) * 255 + .5);
}

bool intersect(const Ray &r, double &t, int &id) {
  double 
    n = sizeof (spheres) / sizeof (Sphere), 
    d, 
    inf = t = 1e20;
  for (int i = int(n); i--;) if ((d = spheres[i].intersect(r)) && d < t) {
    t = d;
    id = i;
  }
  return t<inf;
}

Vec radiance(const Ray &r, int depth, unsigned short *Xi) {
  double t; // distance to intersection
  int id = 0; // id of intersected object
  if (!intersect(r, t, id)) 
    return Vec(); // if miss, return black

  const Sphere &obj = spheres[id]; // the hit object
  Vec 
    x = r.o + r.d*t, 
    n = (x - obj.p).norm(), 
    nl = (n.dot(r.d) < 0) ? n : n*-1, 
    f = obj.c;

  double p = f.x > f.y && f.x > f.z ? f.x : f.y > f.z ? f.y : f.z; // max refl
  if (depth > 100) 
    return obj.e;
  if (++depth > 5) 
    if (erand48(Xi) < p) 
      f = f * (1 / p);
  else 
    return obj.e; //R.R.

  if (obj.refl == DIFF) { // Ideal DIFFUSE reflection
    double 
      r1 = 2 * M_PI * erand48(Xi), 
      r2 = erand48(Xi), 
      r2s = sqrt(r2);

    Vec 
      w = nl, 
      u = ((fabs(w.x) > .1 ? Vec(0, 1) : Vec(1)) % w).norm(), 
      v = w % u;

    Vec d = (u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1 - r2)).norm();

    return obj.e + f.mult(radiance(Ray(x, d), depth, Xi));
  }
  else if (obj.refl == SPEC) // Ideal SPECULAR reflection
    return obj.e + f.mult(radiance(Ray(x, r.d - n * 2 * n.dot(r.d)), depth, Xi));

  Ray reflRay(x, r.d - n * 2 * n.dot(r.d)); // Ideal dielectric REFRACTION
  bool into = n.dot(nl) > 0; // Ray from outside going in?
  double 
    nc = 1, 
    nt = 1.5, 
    nnt = into ? nc / nt : nt / nc, 
    ddn = r.d.dot(nl), 
    cos2t;

  if ((cos2t = 1 - nnt * nnt * (1 - ddn * ddn)) < 0) // Total internal reflection
    return obj.e + f.mult(radiance(reflRay, depth, Xi));

  Vec tdir = (r.d * nnt - n * ((into ? 1 : -1)*(ddn * nnt + sqrt(cos2t)))).norm();

  double 
    a = nt - nc, 
    b = nt + nc, 
    R0 = a * a / (b * b), 
    c = 1 - (into ? -ddn : tdir.dot(n));

  double 
    Re = R0 + (1 - R0) * c * c * c * c*c, 
    Tr = 1 - Re, 
    P = .25 + .5 * Re, 
    RP = Re / P, 
    TP = Tr / (1 - P);

  return obj.e + f.mult(depth > 2 ? (erand48(Xi) < P ? // Russian roulette
    radiance(reflRay, depth, Xi) * RP : radiance(Ray(x, tdir), depth, Xi) * TP) :
    radiance(reflRay, depth, Xi) * Re + radiance(Ray(x, tdir), depth, Xi) * Tr);
}

void savebmp(const char *filename, int w, int h, int dpi, Vec *data);

int main(int argc, char *argv []) {
  clock_t time1, time2;
  time1 = clock();

  int 
    w = argc == 4 ? atoi(argv[2]) : 512,
    h = argc == 4 ? atoi(argv[3]) : 384,
    samps = argc == 2 ? atoi(argv[1]) / 4 : 1; // # samples

  Ray cam(Vec(50, 52, 295.6), Vec(0, -0.042612, -1).norm()); // cam pos, dir
  Vec cx = Vec(w * .5135 / h), cy = (cx % cam.d).norm()*.5135, r, *c;
  c = new Vec[w * h];
  //#pragma omp parallel for schedule(dynamic, 1) private(r)       // OpenMP
  
  for (int y = 0; y < h; y++) { // Loop over image rows
    fprintf(stderr, "\rRendering (%d spp) %5.2f%%", samps * 4, 100. * y / (h - 1));
    for (unsigned short x = 0, Xi[3] = { 0, 0, y * y * y }; x < w; x++) // Loop cols
    for (int sy = 0, i = (h - y - 1) * w + x; sy < 2; sy++) // 2x2 subpixel rows
    for (int sx = 0; sx < 2; sx++, r = Vec()) { // 2x2 subpixel cols
      for (int s = 0; s < samps; s++) {
        double r1 = 2 * erand48(Xi), dx = r1 < 1 ? sqrt(r1) - 1 : 1 - sqrt(2 - r1);
        double r2 = 2 * erand48(Xi), dy = r2 < 1 ? sqrt(r2) - 1 : 1 - sqrt(2 - r2);
        Vec d = cx * (((sx + .5 + dx) / 2 + x) / w - .5) +
          cy * (((sy + .5 + dy) / 2 + y) / h - .5) + cam.d;
        r = r + radiance(Ray(cam.o + d * 140, d.norm()), 0, Xi)*(1. / samps);
      } // Camera rays are pushed ^^^^^ forward to start in interior
      c[i] = c[i] + Vec(clamp(r.x), clamp(r.y), clamp(r.z))*.25;
    }
  }
  char fileTitle[50];
  sprintf(fileTitle,"samples_%d.bmp", samps*4);
  savebmp(fileTitle, w, h, 300, c);
  
  time2 = clock();
  float diff = ((float) time2 - (float) time1) / 1000000;

  printf("\n Total: %g sec \n", diff);
}

void savebmp(const char *filename, int w, int h, int dpi, Vec *data) {
  FILE *f;
  int k = w*h;
  int s = 4 * k;
  int filesize = 54 + s;

  double factor = 39.375;
  int m = static_cast<int>(factor);

  int ppm = dpi*m;

  unsigned char bmpfileheader[14] = { 'B', 'M', 0, 0, 0, 0, 0, 0, 0, 0, 54, 0, 0, 0 };
  unsigned char bmpinfoheader[40] = { 40, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 24, 0 };

  bmpfileheader[2] = (unsigned char) (filesize);
  bmpfileheader[3] = (unsigned char) (filesize >> 8);
  bmpfileheader[4] = (unsigned char) (filesize >> 16);
  bmpfileheader[5] = (unsigned char) (filesize >> 24);

  bmpinfoheader[4] = (unsigned char) (w);
  bmpinfoheader[5] = (unsigned char) (w >> 8);
  bmpinfoheader[6] = (unsigned char) (w >> 16);
  bmpinfoheader[7] = (unsigned char) (w >> 24);

  bmpinfoheader[8] = (unsigned char) (h);
  bmpinfoheader[9] = (unsigned char) (h >> 8);
  bmpinfoheader[10] = (unsigned char) (h >> 16);
  bmpinfoheader[11] = (unsigned char) (h >> 24);

  bmpinfoheader[21] = (unsigned char) (s);
  bmpinfoheader[22] = (unsigned char) (s >> 8);
  bmpinfoheader[23] = (unsigned char) (s >> 16);
  bmpinfoheader[24] = (unsigned char) (s >> 24);

  bmpinfoheader[25] = (unsigned char) (ppm);
  bmpinfoheader[26] = (unsigned char) (ppm >> 8);
  bmpinfoheader[27] = (unsigned char) (ppm >> 16);
  bmpinfoheader[28] = (unsigned char) (ppm >> 24);

  bmpinfoheader[29] = (unsigned char) (ppm);
  bmpinfoheader[30] = (unsigned char) (ppm >> 8);
  bmpinfoheader[31] = (unsigned char) (ppm >> 16);
  bmpinfoheader[32] = (unsigned char) (ppm >> 24);

  f = fopen(filename, "wb");

  fwrite(bmpfileheader, 1, 14, f);
  fwrite(bmpinfoheader, 1, 40, f);

  for (int i = k - 1; i >= 0; i--) {
    Vec rgb(toInt(data[i].x), toInt(data[i].y), toInt(data[i].z));

    double red = (data[i].x);
    double green = (data[i].y);
    double blue = (data[i].z);

    unsigned char color[3] = { rgb.x, rgb.y, rgb.z };

    fwrite(color, 1, 3, f);
  }

  fclose(f);
}