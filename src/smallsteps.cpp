#include <cstdlib>

#if defined(WIN32)
#pragma warning(disable:4996)
#include <GL/glut.h>
#ifdef NDEBUG
#pragma comment(linker, "/subsystem:\"windows\" /entry:\"mainCRTStartup\"")
#endif // NDEBUG
#elif defined(__APPLE__) || defined(MACOSX)
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <GLUT/glut.h>
#else // MACOSX
#include <GL/glut.h>
#endif // unix

#include "glm/glm.hpp"
#include "glm/gtc/quaternion.hpp"
#include "glm/gtx/quaternion.hpp"
#include "glm/ext/quaternion_geometric.hpp"
#include "glm/gtc/type_ptr.hpp"
#include "glm/gtx/norm.hpp"
//#include "glm/ext.hpp"
#include "glm/gtx/string_cast.hpp"
#include "imgui.h"
#include "imgui_impl_glut.h"
#include "imgui_impl_opengl2.h"
#include <vector>
#include <iostream>
#include <cstdint>
#include <algorithm>
#include <cfloat>

#define USE_DOUBLE     (1)
#define USE_BOT        (0)
#define USE_CAPTURE    (0)
#define USE_SUB_STEP   (1)

#if !(USE_DOUBLE)
typedef float     Float;
typedef glm::vec3 Vec3;
typedef glm::quat Quat;
typedef glm::mat3 Mat3;
#else
typedef double       Float;
typedef glm::f64vec3 Vec3;
typedef glm::dquat   Quat;
typedef glm::dmat3   Mat3;
#endif

namespace {
  Float FIXED_DT      = (Float)(1.0 / 60.0);
  Float GRAVITY       = (Float)9.8;
};

struct Plane {
  glm::vec4 v;
};

struct Light {
  glm::vec4 v;
};

struct Shadow {
  float v[4][4];
};

struct Material {
  GLfloat ambient[4];
  GLfloat diffuse[4];
  GLfloat specular[4];
  GLfloat shininess;
};

// brass
Material mat_brass = {
  { 0.329412f,  0.223529f, 0.027451f, 1.0f },
  { 0.780392f,  0.568627f, 0.113725f, 1.0f },
  { 0.992157f,  0.941176f, 0.807843f, 1.0f },
  27.89743616f,
};

// white rubber
Material mat_white_rubber = {
  { 0.05f, 0.05f, 0.05f, 1.0f },
  { 0.5f,  0.5f,  0.5f,  1.0f },
  { 0.7f,  0.7f,  0.7f,  1.0f },
  10.0f,
};

void render_pipe(GLfloat width, GLfloat length, int slice) {
  GLUquadricObj* q;
  q = gluNewQuadric();
  gluQuadricDrawStyle(q, GLU_FILL);
  gluQuadricNormals(q, GLU_SMOOTH);
  gluCylinder(q, width, width, length, slice, 1); // quad, base, top, height, slice, stacks
  glPushMatrix();
    glTranslatef(0.0f, 0.0f, length);
    gluDisk(q, 0.0f, width, slice, 1); // quad, inner, outer, slices, loops(top)
  glPopMatrix();
  gluQuadricOrientation(q, GLU_INSIDE);
  gluDisk(q, 0.0f, width, slice, 1); // quad, inner, outer, slices, loops(bottom)
  gluDeleteQuadric(q); 
}

void render_pipe(glm::vec3& start, glm::vec3& end, GLfloat width, int slice) {
  glm::vec3 vec    = end - start;
  GLfloat   length = glm::length(vec);
  GLfloat   ax;
  if (fabs(vec.z) < FLT_EPSILON) {
    ax = 57.2957795f * acos( vec.x / length ); // rotation angle in x-y plane
    if (vec.y <= 0.0f)
      ax = -ax;
  } else {
    ax = 57.2957795f * acos( vec.z / length ); // rotation angle
    if (vec.z <= 0.0f)
      ax = -ax;
  }
  GLfloat rx = -vec.y * vec.z;
  GLfloat ry =  vec.x * vec.z;
  glPushMatrix();
    glTranslatef(start.x, start.y, start.z);
    if (fabs(vec.z) < FLT_EPSILON) {
      glRotatef(90.0f,  0.0f, 1.0f, 0.0f); // Rotate & align with x axis
      glRotatef(   ax, -1.0f, 0.0f, 0.0f); // Rotate to point 2 in x-y plane
    } else {
      glRotatef(   ax,   rx,    ry, 0.0f); // Rotate about rotation vector
    }
    render_pipe(width, length, slice);
  glPopMatrix();
}

void set_material(Material& mat) {
  glMaterialfv(GL_FRONT, GL_AMBIENT,   mat.ambient);
  glMaterialfv(GL_FRONT, GL_DIFFUSE,   mat.diffuse);
  glMaterialfv(GL_FRONT, GL_SPECULAR,  mat.specular);
  glMaterialfv(GL_FRONT, GL_SHININESS, &mat.shininess);
}

void set_material(Material& in_mat, float alpha) {
  Material mat = in_mat;
  mat.ambient[3]   = alpha;
  mat.diffuse[3]   = alpha;
  mat.specular[3]  = alpha;
  set_material(mat);
}

class Point{
public:
  Float inv_mass;
  Vec3  position;
  Vec3  prev_position;
  Vec3  velocity;
  Float force;
  Float elongation;
  Float damping_coeff;
  Point(Float inv_m, Vec3& pos, Vec3& vel) : inv_mass(inv_m), position(pos), prev_position(pos), velocity(vel), force(0.0), elongation(0.0), damping_coeff(0.0f) {
  };
  Float KineticEnery(Float dt) {
    if (inv_mass < FLT_EPSILON) {
      return (Float)0.0; // fixed
    }
    Float m = 1.0f / (Float)inv_mass;
    return (Float)(1.0 / 2.0) * m * glm::length(velocity) * glm::length(velocity); // 1/2 mv^2
  }
  Float PotentialEnergy(Float origin) {
    if (inv_mass < FLT_EPSILON) {
      return (Float)0.0; // fixed
    }
    Float m = (Float)1.0 / (Float)inv_mass;
    Float h = position.y - origin; // h
    return m * GRAVITY * h; // mgh
  }
  void Predict(Float sdt) {
    if (inv_mass < FLT_EPSILON) {
      return;
    }
    velocity.y -= GRAVITY * sdt;
    prev_position = position;
    position += velocity * sdt;
  }
  void SolveVelocity(Float sdt) {
    Vec3  n  = velocity;
    Float v  = glm::length(n);
    n       /= v; // normalize()
    Float dv = -v * (std::min)(Float(1.0), damping_coeff * sdt * inv_mass);
    velocity += n * dv;
  }
  void Render(GLdouble radius, Material& mat, float alpha) {
    glPushMatrix();
      glTranslatef((GLfloat)position.x, (GLfloat)position.y, (GLfloat)position.z);
      set_material(mat, alpha);
      glFrontFace(GL_CW);
      glutSolidSphere(radius, 16, 16);
      glFrontFace(GL_CCW);
    glPopMatrix();
  }
};

class DistanceConstraint {
public:
  Point*               point0;
  Point*               point1;
  Float                compliance;
  Float                rest_length;
  Float                damping_coeff;
  bool                 unilateral;
  DistanceConstraint(Point* p0, Point* p1, GLfloat in_compliance) : point0(p0), point1(p1), compliance(in_compliance), rest_length(0.0), damping_coeff(0.0), unilateral(false) {
    Vec3 p0_to_p1 = point1->position - point0->position;
    rest_length = glm::length(p0_to_p1);
  }
  void SolvePosition(Float dt){
    Float w = point0->inv_mass + point1->inv_mass;
    if (w < FLT_EPSILON) {
      return;
    }
    Vec3 grad = point1->position - point0->position;
    Float   d    = glm::length(grad);
    grad          /= d; // normalize
    w += compliance / dt / dt;
    Float lambda = (d - rest_length) / w;
    if ((lambda < 0.0f) && unilateral) {
      return;
    }
    point1->force = lambda / dt / dt;
    point1->elongation = glm::length(d) - rest_length;
    point0->position += (grad * Float(+1.0) * point0->inv_mass * lambda);
    point1->position += (grad * Float(-1.0) * point1->inv_mass * lambda);
  }
  void SolveVelocity(Float dt) {
    Vec3 n = point1->position - point0->position;
    n = glm::normalize(n);
    Float v0  = glm::dot(n, point0->velocity);
    Float v1  = glm::dot(n, point1->velocity);
    Float dv0 = (v1 - v0) * std::min(Float(0.5), damping_coeff * dt * point0->inv_mass);
    Float dv1 = (v0 - v1) * std::min(Float(0.5), damping_coeff * dt * point1->inv_mass);
    point0->velocity += n * dv0;
    point1->velocity += n * dv1;
  }
  Float Stretch() {
    return point1->elongation;
  }
  Float Force() {
    // TODO
  }
};

struct Energy {
  Float k;
  Float u;
  Energy(Float in_k, Float in_u) : k(in_k), u(in_u) {}
};

struct DebugInfo {
  bool                show_depth;
  GLfloat             dof;
  GLfloat             focus;
  std::vector<Energy> energies;
  std::vector<Float>  stretches;
  DebugInfo() : show_depth(false), dof(0.1f), focus(0.0f), energies(), stretches() {}
};

struct Context;

class Scene {
public:
  enum {
    ePendulum,
    eCloth,
    eNum,
  };
  virtual void Update(Context& context, Float dt) = 0;
  virtual void Render(float alpha = 1.0f) = 0;
  virtual ~Scene() {}
};

struct Context {
  std::uint32_t       frame;
  Float               time;
  DebugInfo           debug_info;
  Plane               floor;
  Light               light;
  Shadow              floor_shadow;
  Scene*              scene;
  int                 num_links;
  int                 num_sub_steps;
  int                 scene_num;
  bool                use_collision;
  Context() : frame(0), time(0.0f), debug_info(), floor(), light(), floor_shadow(), scene(nullptr), num_links(3), num_sub_steps(10), scene_num(Scene::ePendulum), use_collision(false) {}
};

Context g_Context;

class ScenePendulum : public Scene {
public:
  std::vector<Point>              points;
  std::vector<DistanceConstraint> constraints;
  ScenePendulum(int num_links) : points(), constraints() {
    int particle_num = num_links + 1;
    points.reserve(particle_num);
    for(int i = 0; i < particle_num; i++) {
      Float inv_mass = (i == 0) ? (Float)0.0 : Float(1.0 / 1.0); // m= 1.0(kg)
      Float spacing  = Float(1.0) / (Float)num_links;
      Vec3 position(Float(0.0) + spacing * i, Float(1.0), Float(0.0));
      Vec3 velocity(Float(0.0), Float(0.0), Float(0.0));
      Vec3 gravity(Float(0.0), -GRAVITY, Float(0.0));
      points.push_back(Point(inv_mass, position, velocity));
    }
    Point* prev = &(*points.begin());
    for(auto itr = points.begin(); itr != points.end(); itr++) {
      if (itr != points.begin()) {
        constraints.push_back(DistanceConstraint(prev, &(*itr), 0.0f));
        prev = &(*itr);
      }
    }
  }
  ~ScenePendulum() {
    points.clear();
    points.shrink_to_fit();
    constraints.clear();
    constraints.shrink_to_fit();
  }
  void CalcTotalEnergies(Context& ctx) {
    Float k = 0;
    Float u = 0;
    for (auto& p : points) {
      k += p.KineticEnery(FIXED_DT);
      u += p.PotentialEnergy(0.0f);
    }
    ctx.debug_info.energies.push_back(Energy(k, u));
  }
  void CalcTotalStretches(Context& ctx) {
    Float stretch = (Float)0.0;
    for (auto& c : constraints) {
      stretch += c.Stretch();
    }
    ctx.debug_info.stretches.push_back(stretch);
  }
  virtual void Update(Context& ctx, Float dt) {
#if USE_SUB_STEP
    Float   sdt  = dt / ctx.num_sub_steps;
    for(int step = 0; step < ctx.num_sub_steps; step++) {
      for (auto& p : points) {
        p.Predict(sdt);
      }
      for(auto& c : constraints) {
        c.SolvePosition(sdt);
      }
      if (ctx.use_collision) {
        Float min_x = (Float)0.0; // center
        auto& p = points[ctx.num_links];
        if (p.position.x < min_x) {
          p.position.x = min_x;
          if (p.velocity.x < (Float)0.0) {
            p.prev_position.x = p.position.x + p.velocity.x * sdt;
          }
        }
      }
      for (auto& p : points) {
        p.velocity = p.position - p.prev_position;
        p.velocity *= 1.0f / sdt;
        p.SolveVelocity(sdt);
      }
      for(auto& c : constraints) {
        if (c.compliance > 0.0f) {
          c.SolveVelocity(sdt);
        }
      }
    }
#else
    for (auto& p : points) {
      p.Predict(dt);                                                     // use dt and only one time
    }
    for(int iteration = 0; iteration < ctx.num_sub_steps; iteration++) { // iteration = sub steps num
      for(auto& c : constraints) {
        c.SolvePosition(dt);                                             // use 'dt' instead of 'sdt'
      }
    }
    for (auto& p : points) {                                             // use dt and only one time
      p.velocity = p.position - p.prev_position;
      p.velocity *= 1.0f / dt;
      p.SolveVelocity(dt);
    }
    for(auto& c : constraints) {                                         // use dt and only one time
      if (c.compliance > 0.0f) {
        c.SolveVelocity(dt);
      }
    }
#endif
    CalcTotalEnergies(ctx);
    CalcTotalStretches(ctx);
  }
  virtual void Render(float alpha = 1.0f) {
    GLfloat radius   = 0.03f;
    for(auto itr = points.begin(); itr != points.end(); itr++) {
      (*itr).Render(radius, mat_brass, alpha); // particle
    }
    for(auto itr = constraints.begin(); itr != constraints.end(); itr++) {
      glm::vec3 p0((*itr).point0->position);
      glm::vec3 p1((*itr).point1->position);
      glFrontFace(GL_CW);
      set_material(mat_white_rubber, alpha);
      render_pipe(p0, p1, 0.01f, 12); // constraint
      glFrontFace(GL_CCW);
    }
  }
};

class SceneCloth : public Scene {
public:
  std::vector<Point>              points;
  std::vector<DistanceConstraint> constraints;
  SceneCloth(int w, int h) : points(), constraints() {
    // TODO
  }
  ~SceneCloth() {
    // TODO
  }
  virtual void Update(Context& ctx, Float dt) {
    // TODO
  }
  virtual void Render(float alpha = 1.0f) {
    // TODO
  }
};

void write_ppm(GLubyte* buff, GLenum format) {
  int w = glutGet(GLUT_WINDOW_WIDTH);
  int h = glutGet(GLUT_WINDOW_HEIGHT);
  int  pix_sz = (format == GL_RGBA) ? 4 : 1;
  char suffix[256];
  sprintf(suffix, (format == GL_RGBA) ? "screen.ppm" : "depth.ppm");
  char filename[1024];
  sprintf(filename, "%08d_%s", g_Context.frame, suffix);
  FILE *fp = fopen(filename, "wb");
  if (fp) {
    fprintf(fp, "P%d\n", (format == GL_RGBA) ? 6 : 5); // 5:Portable graymap(Binary), 6:Portable pixmap(Binary)
    fprintf(fp, "%u %u\n", w, h);
    fprintf(fp, "255\n");
    for(int y = 0; y < h; y++) {
      for(int x = 0; x < w; x++) {
        int index = (h - y - 1) * w * pix_sz + (x * pix_sz);
        if (format == GL_RGBA) {
          int r = buff[index];
          int g = buff[index + 1];
          int b = buff[index + 2];
          int a = buff[index + 3]; // not use here
          putc(r, fp); // binary
          putc(g, fp);
          putc(b, fp);
        } else {
          putc(buff[index], fp);
        }
      }
    }
    fclose(fp);
  }
}

void write_image(GLenum format) {
  GLsizei w = glutGet(GLUT_WINDOW_WIDTH);
  GLsizei h = glutGet(GLUT_WINDOW_HEIGHT);
  GLubyte* buff = (GLubyte*)malloc((size_t)w * (size_t)h * 4); // w x h * RGBA
  glReadBuffer(GL_BACK);
  glReadPixels(0, 0, w, h, format, GL_UNSIGNED_BYTE, buff);
  write_ppm(buff, format);
  free(buff);
}

void dump_params(const char* title) {
  char tmpstr[1024];
  sprintf(tmpstr, "%s_out.txt", title);
  FILE* fp = fopen(tmpstr, "w");
  sprintf(tmpstr, "# Time, %s\n", title);
  fprintf(fp, tmpstr);
  Float t = (Float)0.0;
  auto& ctx      = g_Context;
  if (strcmp(title, "Energy") == 0) {
    auto& energies = ctx.debug_info.energies; // energies
    for (const auto& e : energies) {
      fprintf(fp, "%f, %f\n", (float)t, (float)e.k + (float)e.u);
      t += FIXED_DT;
    }
  } else if (strcmp(title, "Stretch") == 0) {
    auto& stretches = ctx.debug_info.stretches; // stretches
    for (const auto& s : stretches) {
      fprintf(fp, "%f, %f\n", (float)t, (float)s);
      t += FIXED_DT;
    }
  }
  fprintf(fp, "\n");
  fclose(fp);
}

void pipe_gnuplot(const char* title) {
  dump_params(title);
#if defined(WIN32)
 FILE* gp = _popen("gnuplot", "w");
#else
 FILE* gp = popen("gnuplot", "w");
#endif
  fprintf(gp, "unset key\n");
  char tmpstr[1024];
  sprintf(tmpstr, "set title'%s'\n", title);
  fprintf(gp, tmpstr);
  fprintf(gp, "set xlabel 'Time'\n");
  sprintf(tmpstr, "set ylabel '%s'\n", title);
  fprintf(gp, tmpstr);
  sprintf(tmpstr, "plot '%s_out.txt' with lines\n", title);
  fprintf(gp, tmpstr);
  sprintf(tmpstr, "save 'smallsteps_%s.plt'\n", title);
  fprintf(gp, tmpstr);
  fprintf(gp, "exit\n");
  fflush(gp);
#if defined(WIN32)
  _pclose(gp);
#else
  pclose(gp);
#endif
}

void render_floor(GLfloat w, GLfloat d, int num_w, int num_d) {
  static const GLfloat color[][4] = { { 0.8f, 0.8f, 0.8f, 1.0f },   // white
                                      { 0.3f, 0.3f, 0.3f, 1.0f } }; // gray
  GLfloat center_w = (w * num_w) / 2.0f;
  GLfloat center_d = (d * num_d) / 2.0f;
  glPushMatrix();
    glNormal3f(0.0f, 1.0f, 0.0f); // up vector
    glBegin(GL_QUADS);
    for (int j = 0; j < num_d; ++j) {
      GLfloat dj  = d  * j;
      GLfloat djd = dj + d;
      for (int i = 0; i < num_w; ++i) {
        GLfloat wi  = w  * i;
        GLfloat wiw = wi + w;
        glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, color[(i + j) & 1]);
        glVertex3f(wi  - center_w,  0.0, dj  - center_d);
        glVertex3f(wi  - center_w,  0.0, djd - center_d);
        glVertex3f(wiw - center_w,  0.0, djd - center_d);
        glVertex3f(wiw - center_w,  0.0, dj  - center_d);
      }
    }
    glEnd();
  glPopMatrix();
}

void init_imgui() {
  IMGUI_CHECKVERSION();
  ImGui::CreateContext();
  ImGui::StyleColorsDark();
  ImGui_ImplGLUT_Init();
  ImGui_ImplGLUT_InstallFuncs();
  ImGui_ImplOpenGL2_Init();
}

void finalize_imgui() {
  ImGui_ImplOpenGL2_Shutdown();
  ImGui_ImplGLUT_Shutdown();
  ImGui::DestroyContext();
}

void finalize(void) {
  finalize_imgui();
  return;
}

void find_plane(Plane* plane, glm::vec3& v0, glm::vec3& v1, glm::vec3& v2) {
  glm::vec3 vec0 = v1 - v0;
  glm::vec3 vec1 = v2 - v0;
  plane->v[0] =   vec0.y * vec1.z - vec0.z * vec1.y;
  plane->v[1] = -(vec0.x * vec1.z - vec0.z * vec1.x);
  plane->v[2] =   vec0.x * vec1.y - vec0.y * vec1.x;
  plane->v[3] = -(plane->v[0] * v0.x + plane->v[1] * v0.y + plane->v[2] * v0.z);
}

void calc_shadow_matrix(Shadow* shadow, Plane& plane, Light& light) {
  GLfloat dot = glm::dot(plane.v, light.v);
  shadow->v[0][0] = dot - light.v[0] * plane.v[0];
  shadow->v[1][0] = 0.f - light.v[0] * plane.v[1];
  shadow->v[2][0] = 0.f - light.v[0] * plane.v[2];
  shadow->v[3][0] = 0.f - light.v[0] * plane.v[3];

  shadow->v[0][1] = 0.f - light.v[1] * plane.v[0];
  shadow->v[1][1] = dot - light.v[1] * plane.v[1];
  shadow->v[2][1] = 0.f - light.v[1] * plane.v[2];
  shadow->v[3][1] = 0.f - light.v[1] * plane.v[3];

  shadow->v[0][2] = 0.f - light.v[2] * plane.v[0];
  shadow->v[1][2] = 0.f - light.v[2] * plane.v[1];
  shadow->v[2][2] = dot - light.v[2] * plane.v[2];
  shadow->v[3][2] = 0.f - light.v[2] * plane.v[3];

  shadow->v[0][3] = 0.f - light.v[3] * plane.v[0];
  shadow->v[1][3] = 0.f - light.v[3] * plane.v[1];
  shadow->v[2][3] = 0.f - light.v[3] * plane.v[2];
  shadow->v[3][3] = dot - light.v[3] * plane.v[3];
}

void initialize(int argc, char* argv[]) {
  glClearColor(0.5f, 0.5f, 0.8f, 1.0f);
  glClearAccum(0.0f, 0.0f, 0.0f, 0.0f); 
  glClearDepth(1.0);
  glDepthFunc(GL_LESS);
  glEnable(GL_CULL_FACE);
  glCullFace(GL_BACK);
  glEnable(GL_DEPTH_TEST);

  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  //glEnable(GL_LIGHT1);

  GLfloat ambient[] = { 0.5, 0.5, 0.5, 1.0 };
  GLfloat diffuse[] = { 1.0, 1.0, 1.0, 1.0 };
  GLfloat specular[] = { 1.0, 1.0, 1.0, 1.0 };
  GLfloat lmodel_ambient[] = { 0.2f, 0.2f, 0.2f, 1.0f };
  GLfloat local_view[] = { 0.0 };

  glLightfv(GL_LIGHT1, GL_AMBIENT, ambient);
  glLightfv(GL_LIGHT1, GL_DIFFUSE, diffuse);
  glLightfv(GL_LIGHT1, GL_DIFFUSE, specular);

  glLightModelfv(GL_LIGHT_MODEL_AMBIENT, lmodel_ambient);
  glLightModelfv(GL_LIGHT_MODEL_LOCAL_VIEWER, local_view);

  GLfloat light0pos[] = { 5.0, 5.0, 5.0, 0.0 }; // w:0 directional
  glLightfv(GL_LIGHT0, GL_POSITION, light0pos);
  g_Context.light.v[0] = light0pos[0];
  g_Context.light.v[1] = light0pos[1];
  g_Context.light.v[2] = light0pos[2];
  g_Context.light.v[3] = light0pos[3];

  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glEnable(GL_BLEND);
  glEnable(GL_POINT_SMOOTH);
  glEnable(GL_LINE_SMOOTH);
  //glEnable(GL_POLYGON_SMOOTH);
  glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
  glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
  //glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);
  glEnable(GL_AUTO_NORMAL);
  glEnable(GL_NORMALIZE);

  init_imgui();
  atexit(finalize);

  glm::vec3 v0(-1.0f, 0.0f,  0.0f);
  glm::vec3 v1(+1.0f, 0.0f,  0.0f);
  glm::vec3 v2(+1.0f, 0.0f, -1.0f);
  find_plane(&g_Context.floor, v0, v1, v2);

  if (g_Context.scene == nullptr) {
    switch(g_Context.scene_num) {
    case Scene::ePendulum: g_Context.scene = new ScenePendulum(g_Context.num_links); break;
    case Scene::eCloth:    g_Context.scene = new SceneCloth(8, 8);                   break;
    }
  }
}

void restart() {
  if (g_Context.scene) {
    delete g_Context.scene;
    g_Context.scene = nullptr;
  }
  auto& energies = g_Context.debug_info.energies;
  energies.clear();
  energies.shrink_to_fit();
  auto& stretches = g_Context.debug_info.stretches;
  stretches.clear();
  stretches.shrink_to_fit();
}

void display_imgui() {
  ImGui_ImplOpenGL2_NewFrame();
  ImGui_ImplGLUT_NewFrame();
  {
    ImGui::SetNextWindowPos(ImVec2(  10,  10), ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(270, 130), ImGuiCond_FirstUseEver);

    ImGui::Begin("Debug");
    ImGui::Checkbox("Show Depth", &g_Context.debug_info.show_depth);
    ImGui::SliderFloat("DoF",     &g_Context.debug_info.dof,     0.0f,  0.2f);
    ImGui::SliderFloat("focus",   &g_Context.debug_info.focus,  -1.5f,  1.5f);
    ImGui::End();

    ImGui::Begin("Params");
    if (ImGui::Button("Restart")) {
      restart();
    }
    if (ImGui::SliderInt("Links",    &g_Context.num_links, 1, 16)) {
      restart();
    }
    ImGui::SliderInt("Substeps", &g_Context.num_sub_steps, 1, 1000);
    ImGui::Checkbox("Collision Handling", &g_Context.use_collision);
    ImGui::End();

    std::vector<Energy>& energies = g_Context.debug_info.energies;
    int sz = (int)energies.size();
    if (sz > 0) {
      ImGui::Begin("Energy");
      std::vector<float> K, U, E;
      for (const auto& e : energies) {
        K.push_back((float)e.k);
        U.push_back((float)e.u);
        E.push_back((float)(e.k + e.u));
      }
      ImGui::PlotLines("K / Time", &K[0], sz);
      ImGui::PlotLines("U / Time", &U[0], sz);
      ImGui::PlotLines("E / Time", &E[0], sz);
      ImGui::End();
    }
  }
  ImGui::Render();
  ImGui_ImplOpenGL2_RenderDrawData(ImGui::GetDrawData());
}

void display_depth() {
  if (g_Context.debug_info.show_depth == false) { return; }
  GLint view[4];
  GLubyte *buffer;
  glGetIntegerv(GL_VIEWPORT, view);
  buffer = (GLubyte *)malloc(size_t(view[2]) * size_t(view[3]));
  glReadPixels(view[0], view[1], view[2], view[3], GL_DEPTH_COMPONENT, GL_UNSIGNED_BYTE, buffer);
  glDisable(GL_DEPTH_TEST);
    glDrawPixels(view[2], view[3], GL_LUMINANCE, GL_UNSIGNED_BYTE, buffer);
  glEnable(GL_DEPTH_TEST);
  free(buffer);
}

void display_actor(float alpha = 1.0f) {
  if (g_Context.scene)
    g_Context.scene->Render(alpha);
}

void display(void){
  glClear(GL_ACCUM_BUFFER_BIT);
  int   num_accum = 8;
  struct jitter_point{ GLfloat x, y; };
  jitter_point j8[] = {
    {-0.334818f,  0.435331f},
    { 0.286438f, -0.393495f},
    { 0.459462f,  0.141540f},
    {-0.414498f, -0.192829f},
    {-0.183790f,  0.082102f},
    {-0.079263f, -0.317383f},
    { 0.102254f,  0.299133f},
    { 0.164216f, -0.054399f}
  };
  GLint viewport[4];
  glGetIntegerv (GL_VIEWPORT, viewport);
  for(int i = 0 ; i < num_accum; i++) {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glm::vec3 pos(0.0f, 1.7f, 4.0f);
    float eye_jitter = (pos.z - g_Context.debug_info.focus) / pos.z;
    eye_jitter = (eye_jitter < 0.1f) ? 0.1f : eye_jitter;
    pos.x += g_Context.debug_info.dof * j8[i].x * eye_jitter;
    pos.y += g_Context.debug_info.dof * j8[i].y * eye_jitter;
    glm::vec3 tgt(0.0f, 0.2f, 0.0f);
    glm::vec3 vec = tgt - pos;
    tgt.y = pos.y + vec.y * ((pos.z - g_Context.debug_info.focus) / pos.z);
    tgt.z = g_Context.debug_info.focus;
    gluLookAt(pos.x, pos.y, pos.z, tgt.x, tgt.y, tgt.z, 0.0, 1.0, 0.0); // pos, tgt, up

    render_floor(0.5f, 0.5f, 16, 8);

    glDisable(GL_DEPTH_TEST);
    glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);
    glEnable(GL_STENCIL_TEST);
    glStencilOp(GL_REPLACE, GL_REPLACE, GL_REPLACE);
    glStencilFunc(GL_ALWAYS, 1, 0xffffffff);
    render_floor(1.0f, 1.0f, 16, 24);       // floor pixels just get their stencil set to 1. 
    glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
    glEnable(GL_DEPTH_TEST);

    glStencilFunc(GL_EQUAL, 1, 0xffffffff); // draw if ==1
    glStencilOp(GL_KEEP, GL_KEEP, GL_KEEP);

    glDisable(GL_DEPTH_TEST);
    glCullFace(GL_BACK);
    glPushMatrix();
      glScalef(1.0f, -1.0f, 1.0f); // for reflection on plane(y=0.0f)
      glLightfv(GL_LIGHT0, GL_POSITION, &g_Context.light.v[0]);
      display_actor(0.1f);
    glPopMatrix();
    glCullFace(GL_FRONT);
    glLightfv(GL_LIGHT0, GL_POSITION, &g_Context.light.v[0]);

    glEnable(GL_POLYGON_OFFSET_FILL);
      calc_shadow_matrix(&g_Context.floor_shadow, g_Context.floor, g_Context.light);
      glDisable(GL_LIGHTING);        // force the 50% black
      glColor4f(0.0, 0.0, 0.0, 0.5f);
      glPushMatrix();
        glMultMatrixf((GLfloat*)g_Context.floor_shadow.v);
        glCullFace(GL_FRONT);
        display_actor();
        glCullFace(GL_BACK);
      glPopMatrix();
      glEnable(GL_LIGHTING);
    glDisable(GL_POLYGON_OFFSET_FILL);
    glDisable(GL_STENCIL_TEST);

    glEnable(GL_DEPTH_TEST);
    glCullFace(GL_FRONT);
    display_actor();               // actual draw
    glCullFace(GL_BACK);

    glAccum(GL_ACCUM, 1.0f / num_accum);
  }
  glAccum(GL_RETURN, 1.0f);

  display_depth();
  display_imgui();

  glutSwapBuffers();
  glutPostRedisplay();
}

void reshape_imgui(int width, int height) {
  ImGuiIO& io = ImGui::GetIO();
  io.DisplaySize.x = (float)width;
  io.DisplaySize.y = (float)height;
}

void reshape(int width, int height){
  glShadeModel(GL_SMOOTH);

  reshape_imgui(width, height);
  glViewport(0, 0, width, height);

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(30.0, (double)width / (double)height, 1.0f, 100.0f);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
}

void keyboard(unsigned char key, int x , int y){
  switch(key) {
  case 's': write_image(GL_RGBA); break;
  case 'd': write_image(GL_DEPTH_COMPONENT); break;
  case 'p': pipe_gnuplot("Stretch"); pipe_gnuplot("Energy"); break;
  case 27: exit(0); break; // esc
  }
}

void idle(void){
  GLfloat time      = (float)glutGet(GLUT_ELAPSED_TIME) / 1000.0f;
  GLfloat dt        = (GLfloat)FIXED_DT;
  auto&   ctx       = g_Context;
  if (ctx.scene == nullptr) {
    switch(g_Context.scene_num) {
    case Scene::ePendulum: ctx.scene = new ScenePendulum(ctx.num_links); break;
    case Scene::eCloth:    ctx.scene = new SceneCloth(8, 8);             break;
    }
  }
  ctx.scene->Update(g_Context, dt);

#if USE_BOT
  if (ctx.frame == 240) {
    keyboard('p', 0, 0); // gnuplot
    keyboard(27, 0, 0); // exit
  }
#endif
#if USE_CAPTURE
  keyboard('s', 0, 0); // screenshot
#endif
  while(1) {
    if (((float)glutGet(GLUT_ELAPSED_TIME) / 1000.0f - time) > FIXED_DT) {
      break; // keep 60fps
    }
  }
  ctx.frame++;
}

void mouse( int button, int state, int x, int y ){
  if (button == GLUT_LEFT_BUTTON) {
    switch(state){
    case GLUT_DOWN:
      break;
    case GLUT_UP:
      break;
    }
  }
}

void motion(int x, int y){

}

int main(int argc, char* argv[]) {
  glutInit(&argc, argv);
#ifdef __FREEGLUT_EXT_H__
  glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE, GLUT_ACTION_GLUTMAINLOOP_RETURNS);
#endif
  glutInitDisplayMode(GLUT_RGBA | GLUT_DEPTH | GLUT_DOUBLE | GLUT_ACCUM | GLUT_STENCIL);
  glutInitWindowSize(640, 480);
  glutCreateWindow("Small Steps in Physics Simulation");

  initialize(argc, argv);

  glutDisplayFunc(display);
  glutReshapeFunc(reshape);
  glutIdleFunc(idle);

  //glutMouseFunc(mouse);     // ImGui_ImplGLUT_MouseFunc
  //glutMotionFunc(motion);   // ImGui_ImplGLUT_MotionFunc
  //glutKeyboardFunc(keyboard); // ImGui_ImplGLUT_KeyboardFunc

  glutMainLoop();
  return 0;
}
